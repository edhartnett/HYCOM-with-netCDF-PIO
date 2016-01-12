      subroutine hybgen(m,n)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- ---------------------
c --- hybrid grid generator
c --- ---------------------
c
      logical    lpipe_hybgen
      parameter (lpipe_hybgen=.false.)
c
      integer   i,j,k,l
      character text*12
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f9.3,f9.2))
cdiag if (itest.gt.0 .and. jtest.gt.0) then
cdiag   write (lp,103) nstep,itest+i0,jtest+j0,
cdiag.  '  entering hybgen:  temp    saln    dens     thkns     dpth',
cdiag.  (k,temp(itest,jtest,k,n),saln(itest,jtest,k,n),
cdiag.  th3d(itest,jtest,k,n)+thbase,dp(itest,jtest,k,n)*qonem,
cdiag.  p(itest,jtest,k+1)*qonem,k=1,kk)
cdiag endif
c
      call xctilr(dpmixl( 1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_ps)
c
      margin = 0  ! no horizontal derivatives
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call hybgenaj(m,n, j)
      enddo
!$OMP END PARALLEL DO
c
c --- vertical momentum flux across moving interfaces (the s-dot term in the
c --- momentum equation) - required to locally conserve momentum when hybgen
c --- moves vertical coordinates first, store old interface pressures in
c --- -pu-, -pv-
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            pu(i,j,1)=0.0
            pu(i,j,2)=dpu(i,j,1,n)
          enddo
          do k=2,kk
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
              pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,k,n)
            enddo
          enddo
        enddo
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            pv(i,j,1)=0.0
            pv(i,j,2)=dpv(i,j,1,n)
          enddo
          do k=2,kk
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k,n)
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
c --- update layer thickness at -u,v- points.
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n),
     &            dpv(1-nbdy,1-nbdy,1,n),
     &            p,depthu,depthv, margin)  ! p's halo extended by dpudpv
c
      if     (lpipe .and. lpipe_hybgen) then
c ---   compare two model runs.
        do k= 1,kk+1
          write (text,'(a9,i3)') 'p      k=',k
          call pipe_compare(p( 1-nbdy,1-nbdy,k),ip,text)
          write (text,'(a9,i3)') 'pu     k=',k
          call pipe_compare(pu(1-nbdy,1-nbdy,k),iu,text)
          write (text,'(a9,i3)') 'pv     k=',k
          call pipe_compare(pv(1-nbdy,1-nbdy,k),iv,text)
        enddo
        do k= 1,kk
          write (text,'(a9,i3)') 'dp     k=',k
          call pipe_compare(dp( 1-nbdy,1-nbdy,k,n),ip,text)
          write (text,'(a9,i3)') 'dpu    k=',k
          call pipe_compare(dpu(1-nbdy,1-nbdy,k,n),iu,text)
          write (text,'(a9,i3)') 'dpv    k=',k
          call pipe_compare(dpv(1-nbdy,1-nbdy,k,n),iv,text)
        enddo
      endif
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call hybgenbj(m,n, j)
      enddo
!$OMP END PARALLEL DO
c
      return
      end
      subroutine hybgenaj(m,n,j )
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c --- --------------------------------------------
c --- hybrid grid generator, single j-row (part A).
c --- --------------------------------------------
c
      include 'common_blocks.h'
c
      integer m,n, j
c
      double precision tsum,ssum,thsum,trsum(mxtrcr),
     &                 psum,q2sum,q2lsum,rpsum
      real ttem(kdm+1),ssal(kdm+1),tthe(kdm+1),ttrc(kdm+1,mxtrcr),
     .     qq2(kdm+1),qq2l(kdm+1),
     .     pres(kdm+2),dp0ij(kdm),dp0cum(kdm+1)
      real pwidth,p_hat,p_hat0,p_hat2,p_hat3,
     .     delt,deltm,dels,delsm,thnew,q,qtr,qts,thkbop
      integer i,k,k1,ka,l,ksubl,kp,ktr,iter
ccc   real colint,colins,colout,colous
cdiag      character*12 cinfo
c
      double precision dsmll
      parameter       (dsmll=1.d-8)
c
c --- c u s h i o n   function (from Bleck & Benjamin, 1992):
c --- if delp >= qqmx*dp0 >>  dp0, -cushn- returns -delp-
c --- if delp <= qqmn*dp0 << -dp0, -cushn- returns  -dp0-
c
      real       qqmn,qqmx,cusha,cushb
      parameter (qqmn=-2.0, qqmx=4.0)  ! traditional range
*     parameter (qqmn=-4.0, qqmx=6.0)  ! somewhat wider range
      parameter (cusha=qqmn**2*(qqmx-1.0)/(qqmx-qqmn)**2)
      parameter (cushb=1.0/qqmn)
c
      real qq,cushn,delp,dp0
      include 'stmt_fns.h'
      qq(   delp,dp0)=max(qqmn,min(qqmx,delp/dp0))
      cushn(delp,dp0)=dp0*
     .                (1.0+cusha*(1.0-cushb*qq(delp,dp0))**2)*
     .                max(1.0,delp/(dp0*qqmx))
c
      thkbop = thkbot*onem
c
      do 1 l=1,isp(j)
c
      do 2 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
      dp0cum(1)=0.
      do k=1,kk
        dp0ij( k)  =min(dp0k(k),max(ds0k(k),dssk(k)*depths(i,j)))
        dp0cum(k+1)=dp0cum(k)+dp0ij(k)
        p(i,j, k+1)=p(i,j,k)+dp(i,j,k,n)
      enddo
c
c --- mark the layer containing the mixed layer base by -ksubl- if the
c --- kraus-turner mixed layer model is selected. otherwise set ksubl to -1.
      if     (mxlkta .and. thermo) then
        ksubl=-1
        do k=1,kk
          if (p(i,j,k  ).lt.dpmixl(i,j,n) .and. 
     .        p(i,j,k+1).ge.dpmixl(i,j,n)+onemm) then
            ksubl=k
          endif
        enddo
        klist(i,j)=ksubl
      else
        ksubl=-1
      endif
c
c --- does layer touch sea floor?
      do 10 k=3,kk
      if (p(i,j,k+1).gt.p(i,j,kk+1)-0.1*onemm) then
        if (dp(i,j,k,n).le.0.1*onemm) then
          if (k.le.nhybrd) then
c ---       fill massless layers on sea floor with fluid from above
            th3d(i,j,k,n)=th3d(i,j,k-1,n)
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=temp(i,j,k-1,n)
          elseif (hybflg.ne.2) then
c ---       fill massless layers on sea floor with saln from above
            th3d(i,j,k,n)=max(theta(k),th3d(i,j,k-1,n))
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          else
c ---       fill massless layers on sea floor with temp from above
            th3d(i,j,k,n)=max(theta(k),th3d(i,j,k-1,n))
            temp(i,j,k,n)=temp(i,j,k-1,n)
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
          do ktr= 1,ntracr
            tracer(i,j,k,n,ktr)=tracer(i,j,k-1,n,ktr)
          enddo
          if (mxlmy) then
            q2( i,j,k,n)=q2( i,j,k-1,n)
            q2l(i,j,k,n)=q2l(i,j,k-1,n)
          endif
*       else if (k.eq.-99) then  ! always .false. at run time
        else if (k.le.nhybrd                              .and.
     &           (theta(k)     -th3d(i,j,k,  n)).gt.0.002 .and.
     &           (th3d(i,j,k,n)-th3d(i,j,k-1,n)).gt.0.002      ) then
c
c ---     water in deepest inflated layer is too light.
c ---     split layer into 2 sublayers, one near the desired density
c ---     and one exactly matching the T&S properties of layer k-1.
c ---     To prevent "runaway" T or S, the result satisfies either
c ---       abs(T.k - T.k-1) <= abs(T.k-2 - T.k-1) or
c ---       abs(S.k - S.k-1) <= abs(S.k-2 - S.k-1) 
c ---     It is also limited to a 50% change in layer thickness.
c
          deltm=abs(temp(i,j,k-2,n)-temp(i,j,k-1,n))
          delt =abs(temp(i,j,k-1,n)-temp(i,j,k,  n))
          delsm=abs(saln(i,j,k-2,n)-saln(i,j,k-1,n))
          dels =abs(saln(i,j,k-1,n)-saln(i,j,k,  n))
          qts=0.0
          if     (delt.gt.epsil) then
            qts=max(qts,(min(deltm,2.0*delt)-delt)/delt)  ! qts<=1.0
          endif
          if     (dels.gt.epsil) then
            qts=max(qts,(min(delsm,2.0*dels)-dels)/dels)  ! qts<=1.0
          endif
          q=(theta(k)-th3d(i,j,k,n))/(theta(k)-th3d(i,j,k-1,n))
          q=min(q,qts/(1.0+qts))  ! upper sublayer <= 50% of total
          p_hat=q*(p(i,j,k+1)-p(i,j,k))
          p(i,j,k)=p(i,j,k)+p_hat
c
          if     (hybflg.eq.0) then  !T&S
            temp(i,j,k,n)=temp(i,j,k,n)+(q/(1.0-q))*(temp(i,j,k,n)  -
     .                                             temp(i,j,k-1,n) )
            saln(i,j,k,n)=saln(i,j,k,n)+(q/(1.0-q))*(saln(i,j,k,n)  -
     .                                             saln(i,j,k-1,n) )
            th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
          elseif (hybflg.eq.1) then  !th&S
            th3d(i,j,k,n)=th3d(i,j,k,n)+(q/(1.0-q))*(th3d(i,j,k,n)  -
     .                                             th3d(i,j,k-1,n) )
            saln(i,j,k,n)=saln(i,j,k,n)+(q/(1.0-q))*(saln(i,j,k,n)  -
     .                                             saln(i,j,k-1,n) )
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
          elseif (hybflg.eq.2) then  !th&T
            th3d(i,j,k,n)=th3d(i,j,k,n)+(q/(1.0-q))*(th3d(i,j,k,n)  -
     .                                             th3d(i,j,k-1,n) )
            temp(i,j,k,n)=temp(i,j,k,n)+(q/(1.0-q))*(temp(i,j,k,n)  -
     .                                             temp(i,j,k-1,n) )
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
          if     (ntracr.gt.0.0) then
            qtr=p_hat/(p(i,j,k+1)-p(i,j,k))  !ok because <1.0 and >0.0
            do ktr= 1,ntracr
              if     (trcflg(ktr).eq.2) then !temperature tracer
                tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)+
     &                             (q/(1.0-q))*(tracer(i,j,k,  n,ktr)-
     &                                          tracer(i,j,k-1,n,ktr))
              else !standard tracer - not split into two sub-layers
                tracer(i,j,k-1,n,ktr)=tracer(i,j,k-1,n,ktr)+
     &                                     qtr*(tracer(i,j,k,  n,ktr)-
     &                                          tracer(i,j,k-1,n,ktr))
              endif
            enddo !ktr
           endif !tracers
          if (mxlmy) then
            q2( i,j,k,n)=q2( i,j,k,n)+
     .            (q/(1.0-q))*(q2( i,j,k,n)-q2( i,j,k-1,n))
            q2l(i,j,k,n)=q2l(i,j,k,n)+
     .            (q/(1.0-q))*(q2l(i,j,k,n)-q2l(i,j,k-1,n))
          endif
cdiag          if (i.eq.itest .and. j.eq.jtest) then
cdiag            write(lp,'(a,i3,f6.3,5f8.3)')
cdiag     .        'hybgen, 10(+):',
cdiag     .        k,q,temp(i,j,k,n),saln(i,j,k,n),
cdiag     .            th3d(i,j,k,n)+thbase,theta(k)+thbase
cdiag            call flush(lp)
cdiag          endif
cdiag        else
cdiag          if (i.eq.itest .and. j.eq.jtest) then
cdiag            write(lp,'(a,i3,f6.3,5f8.3)')
cdiag     .        'hybgen, 10(-):',
cdiag     .        k,0.0,temp(i,j,k,n),saln(i,j,k,n),
cdiag     .            th3d(i,j,k,n)+thbase,theta(k)+thbase
cdiag            call flush(lp)
cdiag          endif
        endif
      endif
 10   continue
c
c --- store one-dimensional arrays of -temp-, -saln-, and -p- for the 'old'
c --- vertical grid before attempting to restore isopycnic conditions
      pres(1)=p(i,j,1)
      kp=0
      do 3 k=1,kk
      k1=k+kp
      if (k.ne.ksubl) then
        tthe(k1)=th3d(i,j,k,n)
        ttem(k1)=temp(i,j,k,n)
        ssal(k1)=saln(i,j,k,n)
        do ktr= 1,ntracr
          ttrc(k1,ktr)=tracer(i,j,k,n,ktr)
        enddo
        if (mxlmy) then
          qq2( k1)=q2( i,j,k,n)
          qq2l(k1)=q2l(i,j,k,n)
        endif
        pres(k1+1)=p(i,j,k+1)
      else				!  k = ksubl
c --- expand layer into two sublayers, one above and one below mixed layer base
        tthe(k)=tthe(k-1)
        ttem(k)=ttem(k-1)
        ssal(k)=ssal(k-1)
        pres(k+1)=dpmixl(i,j,n)
        q=(p(i,j,k)-dpmixl(i,j,n))/(p(i,j,k+1)-dpmixl(i,j,n))
        tthe(k+1)=th3d(i,j,k,n)+q*(tthe(k)-th3d(i,j,k,n))
        ttem(k+1)=temp(i,j,k,n)+q*(ttem(k)-temp(i,j,k,n))
        ssal(k+1)=saln(i,j,k,n)+q*(ssal(k)-saln(i,j,k,n))
        pres(k+2)=p(i,j,k+1)
        do ktr= 1,ntracr
          ttrc(k,  ktr)=ttrc(k-1,ktr)
          ttrc(k+1,ktr)=tracer(i,j,k,n,ktr)+
     &                    q*(ttrc(k,ktr)-tracer(i,j,k,n,ktr))
        enddo
        kp=1
      end if
 3    continue
c
ccc   colint=temp(i,j,1,n)*(p(i,j,2)-p(i,j,1))
ccc   colins=saln(i,j,1,n)*(p(i,j,2)-p(i,j,1))
c
c --- try to restore isopycnic conditions by moving layer interfaces
c
      do 88 k=2,nhybrd
cdiag            if (i.eq.itest .and. j.eq.jtest) then
cdiag              write(cinfo,'(a9,i2.2,1x)') '  do 88 k=',k
cdiag 109          format (i9,2i5,a,a/(i9,8x,a,a,i3,f9.2,f8.2,f9.2,f8.2))
cdiag              write(lp,109) nstep,itest+i0,jtest+j0,
cdiag     .         cinfo,':    othkns  odpth    nthkns  ndpth',
cdiag     .        (nstep,cinfo,':',k1,
cdiag     .         (pres(k1+1)-
cdiag     .          pres(k1)   )*qonem,
cdiag     .          pres(k1+1)  *qonem,
cdiag     .         (p(itest,jtest,k1+1)-
cdiag     .          p(itest,jtest,k1)   )*qonem,
cdiag     .          p(itest,jtest,k1+1)  *qonem,k1=1,kk)
cdiag              call flush(lp)
cdiag            endif
c
ccc   colint=colint+temp(i,j,k,n)*(p(i,j,k+1)-p(i,j,k))
ccc   colins=colins+saln(i,j,k,n)*(p(i,j,k+1)-p(i,j,k))
c
c --- maintain constant thickness in layer 1
      if (k.eq.2) then
        p_hat=-999.0*dp0k(1)
        go to 9
      end if
c
c --- are we dealing with a massless layer on the sea floor?
      if (p(i,j,k).gt.p(i,j,kk+1)-0.1*onemm) then
        if (k.le.nhybrd) then
c ---     fill massless layers on sea floor with fluid from above
          th3d(i,j,k,n)=th3d(i,j,k-1,n)
          saln(i,j,k,n)=saln(i,j,k-1,n)
          temp(i,j,k,n)=temp(i,j,k-1,n)
        elseif (th3d(i,j,k,n).ne.theta(k)) then
          if (hybflg.ne.2) then
c ---       fill with saln from above
            th3d(i,j,k,n)=max(theta(k),th3d(i,j,k-1,n))
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          else
c ---       fill with temp from above
            th3d(i,j,k,n)=max(theta(k),th3d(i,j,k-1,n))
            temp(i,j,k,n)=temp(i,j,k-1,n)
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
        endif
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr)=tracer(i,j,k-1,n,ktr)
        enddo
        if (mxlmy) then
          q2( i,j,k,n)=q2( i,j,k-1,n)
          q2l(i,j,k,n)=q2l(i,j,k-1,n)
        endif
      endif
c
c --- is density noticeably different from isopycnic reference value?
      if (abs(th3d(i,j,k,n)-theta(k)).lt.0.002) go to 8
c
      if (th3d(i,j,k,n).le.theta(k)) go to 7		!  layer too light
c
c --- water in layer k is too dense. try to dilute with water from layer k-1
c
c --- if layer k-1 is too light, perform this adjustment every other pair
c --- of time steps beneath the z-coordinate domain. otherwise, layers
c --- k-1 and k may not be able to reach their target densities
      if (mod(nstep-1,4).le.1 .and.
     .    p(i,j,k).gt.dp0cum(k)+onem .and.
     .    th3d(i,j,k-1,n).lt.theta(k-1)    ) go to 8
c
cdiag      if (i.eq.itest .and. j.eq.jtest) then
cdiag        write(lp,'(a,f8.4)') 'hybgen, too dense:',th3d(i,j,k,n)-theta(k)
cdiag        call flush(lp)
cdiag      endif
c
      if     ((theta(k)-th3d(i,j,k-1,n)).le.epsil) then
        p_hat=min(p(i,j,k-1),
     .            p(i,j,k)-999.0*dp0k(k-1))  ! take entire layer k-1
      else
        q=(theta(k)-th3d(i,j,k,n))/(theta(k)-th3d(i,j,k-1,n))  ! -ve
        p_hat=p(i,j,k)+q*(p(i,j,k+1)-p(i,j,k))
c
c ---   correct for cabbeling by performing an iterative procedure to make
c ---   sure that the target density is achieved to a tolerance of 0.0001.
c ---   if the water becomes more dense, do not move the interface
c
        if     (hybflg.eq.0 .and.
     .          p_hat.gt.p(i,j,k-1) .and. 
     .          p(i,j,k).gt.dp0cum(k)+onem) then
          do 58 iter=1,5
            q=(p(i,j,k)-p_hat)/max(p(i,j,k+1)-p_hat,onemm)
            delt=q*(temp(i,j,k-1,n)-temp(i,j,k,n))
            dels=q*(saln(i,j,k-1,n)-saln(i,j,k,n))
            thnew=sig(temp(i,j,k,n)+delt,saln(i,j,k,n)+dels)-thbase
            if (thnew.gt.th3d(i,j,k,n)) then  ! more dense
              p_hat=p(i,j,k)
              go to 158
            endif
            if (abs(thnew-theta(k)).le.0.0001) then
              go to 158
            endif
            q=(theta(k)-thnew)/(theta(k)-th3d(i,j,k-1,n))
            p_hat=min(p(i,j,k),p_hat+q*(p(i,j,k+1)-p_hat))
  58      continue
 158      continue
        end if
      end if
c
c --- maintain minimum layer thickess of layer k-1.
c
 9    continue
      p_hat0=p_hat
      p_hat=max(p(i,j,k-1)+dp0ij(k-1),
     &      min(p(i,j,k+1)-dp0ij(k),
     &          p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1)) ))
cdiag      if (i.eq.itest .and. j.eq.jtest) then
cdiag        write(lp,'(a,2f8.2)') 'hybgen, 9: ',
cdiag     &   (p_hat0-p(i,j,k-1))*qonem,
cdiag     &   cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))*qonem
cdiag        call flush(lp)
cdiag      endif
c
c --- if isopycnic conditions cannot be achieved because of a blocking
c --- layer in the interior ocean, move interface k-1 (and k-2 if
c --- necessary) upward
c
      if     (k.le.2) then
c ---   do nothing.
      else if (p_hat.ge.p(i,j,k) .and. 
     &         p(i,j,k-1).gt.dp0cum(k-1)+onem .and.
     &         (p(i,j,kk+1)-p(i,j,k-1).lt.thkbop .or.
     &          p(i,j,k-1) -p(i,j,k-2).gt.qqmx*dp0ij(k-2))) then  ! k.gt.2
        p_hat2=p(i,j,k-2)+cushn(p(i,j,k-1)-p_hat+p_hat0-p(i,j,k-2),
     &                          dp0ij(k-2))
        if (p_hat2.lt.p(i,j,k-1)-onemm) then
          p(i,j,k-1)=max(p_hat2,2.0*p(i,j,k-1)-p_hat)
cdiag          if (i.eq.itest .and. j.eq.jtest) then
cdiag            write(lp,'(a,i3.2,f8.2)') 'hybgen,  1blocking :',
cdiag     .                              k-1,p(i,j,k-1)*qonem
cdiag            call flush(lp)
cdiag          endif
          p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
        else if (k.le.3) then
c ---     do nothing.
        else if (p(i,j,k-2).gt.dp0cum(k-2)+onem .and.
     &           (p(i,j,kk+1)-p(i,j,k-2).lt.thkbop .or.
     &            p(i,j,k-2) -p(i,j,k-3).gt.qqmx*dp0ij(k-3))) then  ! k.gt.3
          p_hat3=p(i,j,k-3)+cushn(p(i,j,k-2)-p_hat+p_hat0-p(i,j,k-3),
     &                            dp0ij(k-3))
          if (p_hat3.lt.p(i,j,k-2)-onemm) then
            p(i,j,k-2)=max(p_hat3,2.0*p(i,j,k-2)-p(i,j,k-1))
cdiag            if (i.eq.itest .and. j.eq.jtest) then
cdiag              write(lp,'(a,i3.2,f8.2)') 'hybgen,  2blocking :',
cdiag     .                                k-2,p(i,j,k-2)*qonem
cdiag              call flush(lp)
cdiag            endif
            p_hat2=p(i,j,k-2)+cushn(p(i,j,k-1)-p_hat+p_hat0-p(i,j,k-2),
     &                              dp0ij(k-2))
            if (p_hat2.lt.p(i,j,k-1)-onemm) then
              p(i,j,k-1)=max(p_hat2,2.0*p(i,j,k-1)-p_hat)
cdiag              if (i.eq.itest .and. j.eq.jtest) then
cdiag                write(lp,'(a,i3.2,f8.2)') 'hybgen,  3blocking :',
cdiag     .                                  k-1,p(i,j,k-1)*qonem
cdiag                call flush(lp)
cdiag              endif
              p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
            end if
          end if
        end if
      end if
c
      if (p_hat.le.p(i,j,k)) then
c
c --- upper intfc moves up. entrain layer k-1 water into layer k
c
        p(i,j,k)=p_hat
cdiag        if (i.eq.itest .and. j.eq.jtest) then
cdiag          write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k) :',
cdiag     .                            k,p(i,j,k)*qonem
cdiag          call flush(lp)
cdiag        endif
c
      else				!  p_hat > p(i,j,k)
c
c --- move upper interface down and entrain layer k water into layer k-1
c
        p_hat=min(p_hat,p(i,j,k+1))
        p(i,j,k)=p_hat
cdiag        if (i.eq.itest .and. j.eq.jtest) then
cdiag          write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k-):',
cdiag     .                            k,p(i,j,k)*qonem
cdiag          call flush(lp)
cdiag        endif
c
c --- do we need to inflate layer k by lowering  l o w e r  interface?
c
        if (k.lt.kk) then
          p_hat=p(i,j,k+1)
          go to 6
        end if
      end if
      go to 8
c
c --- water in layer k is too light. try to dilute with water from layer k+1
c
 7    continue
c
c --- is this layer touching the sea floor?
      if (p(i,j,k+1).gt.p(i,j,kk+1)-0.1*onemm .and. k.gt.2) go to 8
c
c --- are we below any KT mixed layer?
      if (mxlkta .and. k.le.klist(i,j)) go to 8
c
c --- if layer k+1 is too dense, perform this adjustment every other pair
c --- of time steps beneath the z-coordinate domain. otherwise, layers
c --- k and k+1 may not be able to reach their target densities
      if (mod(nstep-1,4).ge.2 .and.
     .    p(i,j,k+1).gt.dp0cum(k+1)+onem .and.
     .    th3d(i,j,k+1,n).gt.theta(k+1)        ) go to 8
c
cdiag      if (i.eq.itest .and. j.eq.jtest) then
cdiag        write(lp,'(a,f8.4)') 'hybgen, too light:',theta(k)-th3d(i,j,k,n)
cdiag        call flush(lp)
cdiag      endif
c
      if     ((th3d(i,j,k+1,n)-theta(k)).le.epsil) then
        p_hat=max(p(i,j,k+2),
     .            p(i,j,k+1)+999.0*dp0k(k))  ! take entire layer k+1
      else
        q=(th3d(i,j,k,n)-theta(k))/(th3d(i,j,k+1,n)-theta(k))  ! +ve
        p_hat=p(i,j,k+1)+q*(p(i,j,k)-p(i,j,k+1))
c
c ---   correct for cabbeling by performing an iterative procedure to make
c ---   sure that the target density is achieved to a tolerance of 0.0001.
c ---   if the water becomes less dense, do not move the interface
c
        if     (hybflg.eq.0 .and.
     .          p_hat.lt.p(i,j,k+2) .and. 
     .          p(i,j,k+1).gt.dp0cum(k+1)) then
          do 59 iter=1,5
            q=(p_hat-p(i,j,k+1))/max(p_hat-p(i,j,k),onemm)
            delt=q*(temp(i,j,k+1,n)-temp(i,j,k,n))
            dels=q*(saln(i,j,k+1,n)-saln(i,j,k,n))
            thnew=sig(temp(i,j,k,n)+delt,saln(i,j,k,n)+dels)-thbase
            if (thnew.lt.th3d(i,j,k,n)) then  ! less dense
              p_hat=p(i,j,k+1)
              go to 159
            endif
            if (abs(thnew-theta(k)).le.0.0001) then
              go to 159
            endif
            q=(thnew-theta(k))/(th3d(i,j,k+1,n)-theta(k))
            p_hat=max(p(i,j,k+1),p_hat+q*(p(i,j,k)-p_hat))
  59      continue
 159      continue
        endif
      endif
c
c --- maintain minimum thickess of layers k and k+1
 6    continue
      if (mxlkta .and. k.eq.klist(i,j)) go to 8
      if     (k.lt.nhybrd) then
        p_hat=p(i,j,k+2)-cushn(p(i,j,k+2)-p_hat,dp0ij(k+1))
      endif
      p_hat=p(i,j,k)  +  max(p_hat    -p(i,j,k),dp0ij(k))
      p_hat=min(p_hat,
     &          max(0.5*(p(i,j,k+1)+p(i,j,k+2)),
     &              p(i,j,k+2)-dp0ij(k+1)))
c
      if (p_hat.gt.p(i,j,k+1)) then
c
c --- entrain layer k+1 water into layer k
c
        p(i,j,k+1)=p_hat
cdiag        if (i.eq.itest .and. j.eq.jtest) then
cdiag          write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k+):',
cdiag     .                            k+1,p(i,j,k+1)*qonem
cdiag          call flush(lp)
cdiag        endif
      end if
 8    continue
c
c --- if layer above is too thin, move upper interface down
      p_hat0=p(i,j,k-1)+dp0ij(k-1)
      if (p_hat0.gt.p(i,j,k)   .and.
     &    p_hat0.lt.p(i,j,kk+1)     ) then
        p(i,j,k)=p_hat0
      endif
c
c --- inforce interface order (is this necessary?), usually inexpensive.
      do ka= k+1,kk
        if     (p(i,j,ka).ge.p(i,j,k)) then
          exit  ! usually get here quickly
        else
          p(i,j,ka) = p(i,j,k)
        endif
      enddo
 88   continue
c
ccc   write (lp,'(a,1p,2e16.8,e9.1)') 'temp. column integral:',
ccc  .  colint,colout,(colout-colint)/colint
ccc   write (lp,'(a,1p,2e16.8,e9.1)') 'saln. column integral:',
ccc  .  colins,colous,(colous-colins)/colins
c
c --- remap -temp-, -saln- (and tracer) profiles from the 'old' vertical
c --- grid onto the 'new' vertical grid
c --- update tracer every other time step
      do 22 k=1,kk
      dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
      if     (pres(k)  .eq.p(i,j,k  ) .and.
     .        pres(k+1).eq.p(i,j,k+1)      ) then
        go to 22  ! layer unchanged
      endif
      psum =0.0d0
      thsum=0.0d0
      tsum =0.0d0
      ssum =0.0d0
      q2sum =0.0d0
      q2lsum=0.0d0
      trsum(:)=0.0d0
      do 24 k1=1,kk+kp
        if     (pres(k1+1).le.p(i,j,k  )) then
          go to  24  ! too shallow
        elseif (pres(k1)  .ge.p(i,j,k+1)) then
          go to 124  ! too deep
        elseif (pres(k1)  .eq.p(i,j,k  ) .and.
     .          pres(k1+1).eq.p(i,j,k+1)      ) then
          if     (tthe(k1)  .eq.temp(i,j,k,n) .and.
     .            ssal(k1)  .eq.saln(i,j,k,n)      ) then
            go to 22  ! layer unchanged
          endif
          psum =1.0d0
          thsum=tthe(k1)
          tsum =ttem(k1)
          ssum =ssal(k1)
          do ktr= 1,ntracr
            trsum(ktr)=ttrc(k1,ktr)
          enddo
          if (mxlmy) then
            q2sum =qq2( k1)
            q2lsum=qq2l(k1)
          endif
cdiag          if (i.eq.itest .and. j.eq.jtest) then
cdiag            write(lp,'(a,2i3,1p5e13.4)')
cdiag     .        'hybgen, 24:',k,k1,0.0,   psum,thsum,tsum,ssum
cdiag            call flush(lp)
cdiag          endif
          go to 124  ! exact match, but new layer (very unlikely)
        else         ! old and new layers intersect
          pwidth=min(pres(k1+1),p(i,j,k+1))-
     .           max(pres(k1  ),p(i,j,k  ))
          psum =psum +         pwidth
          thsum=thsum+tthe(k1)*pwidth
          tsum =tsum +ttem(k1)*pwidth
          ssum =ssum +ssal(k1)*pwidth
          do ktr= 1,ntracr
            trsum(ktr)=trsum(ktr)+ttrc(k1,ktr)*pwidth
          enddo
          if (mxlmy) then
            q2sum =q2sum +qq2( k1)*pwidth
            q2lsum=q2lsum+qq2l(k1)*pwidth
          endif
cdiag          if (i.eq.itest .and. j.eq.jtest) then
cdiag            write(lp,'(a,2i3,1p5e12.4)')
cdiag     .        'hybgen, 24:',k,k1,pwidth,psum,thsum,tsum,ssum
cdiag            call flush(lp)
cdiag          endif
        endif
  24  continue
 124  continue
      if(psum.gt.0.0d0) then
        rpsum=1.0d0/psum
        if     (hybflg.eq.0) then  !T&S
          temp(i,j,k,n)=tsum*rpsum
          saln(i,j,k,n)=ssum*rpsum
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
        else
          th3d(i,j,k,n)=thsum*rpsum
          if     (abs(th3d(i,j,k,n)-theta(k)).gt.0.001) then  !still T&S
            temp(i,j,k,n)=tsum*rpsum
            saln(i,j,k,n)=ssum*rpsum
            th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
          elseif (hybflg.eq.1) then  !th&S, layer is approx. isopycnal
            saln(i,j,k,n)=ssum*rpsum
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
          elseif (hybflg.eq.2) then  !th&T, layer is approx. isopycnal
            temp(i,j,k,n)=tsum*rpsum
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
        endif
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr)=trsum(ktr)*rpsum
        enddo
        if (mxlmy) then
          q2( i,j,k,n)=max(dsmll,q2sum *rpsum)
          q2l(i,j,k,n)=max(dsmll,q2lsum*rpsum)
        endif
      else
        if (k.le.nhybrd) then
c ---     fill with fluid from above
          th3d(i,j,k,n)=th3d(i,j,k-1,n)
          saln(i,j,k,n)=saln(i,j,k-1,n)
          temp(i,j,k,n)=temp(i,j,k-1,n)
        elseif (hybflg.ne.2) then
c ---     fill with saln from above
          th3d(i,j,k,n)=theta(k)
          saln(i,j,k,n)=saln(i,j,k-1,n)
          temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
          saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
        else
c ---     fill with temp from above
          th3d(i,j,k,n)=theta(k)
          temp(i,j,k,n)=temp(i,j,k-1,n)
          saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
        endif
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr)=tracer(i,j,k-1,n,ktr)
        enddo
        if (mxlmy) then
          q2( i,j,k,n)=q2( i,j,k-1,n)
          q2l(i,j,k,n)=q2l(i,j,k-1,n)
        endif
      end if
 22   continue
c
cdiag 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f9.3,f9.2))
cdiag      if (i.eq.itest .and. j.eq.jtest) then
cdiag        write (lp,103) nstep,itest+i0,jtest+j0,
cdiag     .  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
cdiag     .  (k,ttem(k),ssal(k),tthe(k)+thbase,
cdiag     .   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag     .   k,temp(i,j,k,n),saln(i,j,k,n),
cdiag     .   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
cdiag     .   p(i,j,k+1)*qonem,
cdiag     .  k=1,kk),
cdiag     .  (k,ttem(k),ssal(k),tthe(k)+thbase,
cdiag     .   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag     .   k=kk+1,kp)
cdiag        call flush(lp)
cdiag      endif
c
 2    continue
c
c --- to avoid roundoff errors in -dpudpv- after restart, make sure p=p(dp)
c
      do 1 k=1,kk
      do 1 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
 1    continue
c
      if(mxlkta) then
        do 71 l=1,isp(j)
        do 71 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
        do 71 k=1,kk
        if(dpmixl(i,j,n).gt.p(i,j,k  ) .and.
     .     dpmixl(i,j,n).le.p(i,j,k+1)) then
          t1sav(i,j,n)=tmix(i,j)
          s1sav(i,j,n)=smix(i,j)
          tmlb(i,j,n)=temp(i,j,k,n)
          smlb(i,j,n)=saln(i,j,k,n)
          nmlb(i,j,n)=k
        end if
c
 71     continue
      end if
c
      return
      end
      subroutine hybgenbj(m,n, j)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c --- --------------------------------------------
c --- hybrid grid generator, single j-row (part B).
c --- --------------------------------------------
c
      include 'common_blocks.h'
c
      integer m,n, j
c
      integer i,k,k1,kp,ksubl,l
      real    uu(kdm+1),vv(kdm+1),pres(kdm+2)
      real    q,thknss,usum,vsum,psum,pwidth
c
c --- vertical momentum flux across moving interfaces (the s-dot term in the
c --- momentum equation) - required to locally conserve momentum when hybgen
c --- moves vertical coordinates.
c
      do 412 l=1,isu(j)
      do 412 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
      thknss=max(thkmin*onem,min(depthu(i,j)-onem,
     .       .5*(dpmixl(i,j,n)+dpmixl(i-1,j,n))))
c
c --- mark the layer containing the mixed layer base by -ksubl- if the
c --- kraus-turner mixed layer model is selected. otherwise set ksubl to -1.
      if     (mxlkta .and. thermo) then
        ksubl=-1
        do k=1,kk
          if (pu(i,j,k  ).lt.thknss .and.
     .        pu(i,j,k+1).ge.thknss+onemm) then
            ksubl=k
          endif
        enddo
      else
        ksubl=-1
      endif
c
c --- store one-dimensional arrays of -u- and -p- for the 'old' vertical grid
      pres(1)=pu(i,j,1)
      kp=0
      do 414 k=1,kk
      k1=k+kp
      if (k.ne.ksubl) then
        uu(k1)=u(i,j,k,n)
        pres(k1+1)=pu(i,j,k+1)
      else                              !  k = ksubl
c --- expand layer into two sublayers, one above and one below mixed layer base
        uu(k)=uu(k-1)
        pres(k+1)=thknss
        q=(pu(i,j,k)-thknss)/(pu(i,j,k+1)-thknss)
        uu(k+1)=u(i,j,k,n)+q*(uu(k)-u(i,j,k,n))
        pres(k+2)=pu(i,j,k+1)
        kp=1
      end if
 414  continue
c
c --- remap -u- profiles from the 'old' vertical grid onto the
c --- 'new' vertical grid.
      do 415 k=1,kk
      pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,k,n)
      usum=0.
      psum=0.
      do 416 k1=1,kk+kp
        if     (pres(k1+1).le.pu(i,j,k  )) then
          go to 416  ! too shallow
        elseif (pres(k1)  .ge.pu(i,j,k+1)) then
          go to 417  ! too deep
        elseif (pres(k1)  .eq.pu(i,j,k  ) .and.
     .          pres(k1+1).eq.pu(i,j,k+1)      ) then
          psum=1.0
          usum=uu(k1)
          go to 417  ! exact match
        else         ! old and new layers intersect
          pwidth=min(pres(k1+1),pu(i,j,k+1))-
     .           max(pres(k1  ),pu(i,j,k  ))
          psum=psum+       pwidth
          usum=usum+uu(k1)*pwidth
        endif
 416  continue
 417  continue
      if(psum.gt.0.0) then
        u(i,j,k,n)=usum/psum
      endif
 415  continue
c
 104  format (i9,2i5,a/(33x,i3,f8.3,f9.3,f9.2))
cdiag if (i.eq.itest .and. j.eq.jtest) then
cdiag   write (lp,104) nstep,itest+i0,jtest+j0,
cdiag.  '   hybgen, do 412:  u       thkns     dpth',
cdiag.  (k,uu(k),
cdiag.   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag.   k,u(i,j,k,n),
cdiag.   dpu(i,j,k,n)*qonem,pu(i,j,k+1)*qonem,
cdiag.   k=1,kk),
cdiag.  (k,uu(k),
cdiag.   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag.   k=kk+1,kp)
cdiag endif
c
 412  continue
c
      do 512 l=1,isv(j)
      do 512 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
      thknss=max(thkmin*onem,min(depthv(i,j)-onem,
     .       .5*(dpmixl(i,j,n)+dpmixl(i,j-1,n))))
c
c --- mark the layer containing the mixed layer base by -ksubl- if the
c --- kraus-turner mixed layer model is selected. otherwise set ksubl to -1.
      if     (mxlkta .and. thermo) then
        ksubl=-1
        do k=1,kk
          if (pv(i,j,k  ).lt.thknss .and.
     .        pv(i,j,k+1).ge.thknss+onemm) then
            ksubl=k
          endif
        enddo
      else
        ksubl=-1
      endif
c
c --- store one-dimensional arrays of -v- and -p- for the 'old' vertical grid
      pres(1)=pv(i,j,1)
      kp=0
      do 514 k=1,kk
      k1=k+kp
      if (k.ne.ksubl) then
        vv(k1)=v(i,j,k,n)
        pres(k1+1)=pv(i,j,k+1)
      else                              !  k = ksubl
c --- expand layer into two sublayers, one above and one below mixed layer base
        vv(k)=vv(k-1)
        pres(k+1)=thknss
        q=(pv(i,j,k)-thknss)/(pv(i,j,k+1)-thknss)
        vv(k+1)=v(i,j,k,n)+q*(vv(k)-v(i,j,k,n))
        pres(k+2)=pv(i,j,k+1)
        kp=1
      end if
 514  continue
c
c --- remap -v- profiles from the 'old' vertical grid onto the
c --- 'new' vertical grid.
      do 515 k=1,kk
      pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k,n)
      vsum=0.
      psum=0.
      do 516 k1=1,kk+kp
        if     (pres(k1+1).le.pv(i,j,k  )) then
          go to 516  ! too shallow
        elseif (pres(k1)  .ge.pv(i,j,k+1)) then
          go to 517  ! too deep
        elseif (pres(k1)  .eq.pv(i,j,k  ) .and.
     .          pres(k1+1).eq.pv(i,j,k+1)      ) then
          psum=1.0
          vsum=vv(k1)
          go to 517  ! exact match
        else         ! old and new layers intersect
          pwidth=min(pres(k1+1),pv(i,j,k+1))-
     .           max(pres(k1  ),pv(i,j,k  ))
          psum=psum+       pwidth
          vsum=vsum+vv(k1)*pwidth
        endif
 516  continue
 517  continue
      if(psum.gt.0.0) then
        v(i,j,k,n)=vsum/psum
      endif
 515  continue
c
cdiag if (i.eq.itest .and. j.eq.jtest) then
cdiag   write (lp,104) nstep,itest+i0,jtest+j0,
cdiag.  '   hybgen, do 512:  v       thkns     dpth',
cdiag.  (k,vv(k),
cdiag.   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag.   k,v(i,j,k,n),
cdiag.   dpv(i,j,k,n)*qonem,pv(i,j,k+1)*qonem,
cdiag.   k=1,kk),
cdiag.  (k,vv(k),
cdiag.   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag.   k=kk+1,kp)
cdiag endif
c
 512  continue
c
      return
      end
c
c
c> Revision history:
c>
c> Feb. 2000 -- total rewrite to convert to 'newzp' approach
c> Jul. 2000 -- added hybgenj for OpenMP parallelization
c> Oct. 2000 -- added hybgenbj to simplify OpenMP logic
c> Nov. 2000 -- fill massless layers on sea floor with salinity from above
c> Nov. 2000 -- unmixing of deepest inflated layer uses th&T&S from above
c> Nov. 2000 -- ignored isopycnic variance is now 0.002
c> Nov. 2000 -- iterate to correct for cabbeling
c> Nov. 2000 -- allow for "blocking" interior layers
c> Nov. 2000 -- hybflg selects conserved fields (any two of T/S/th)
