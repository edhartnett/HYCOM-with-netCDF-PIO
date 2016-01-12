      subroutine inicon(mnth)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer mnth
c
      real    pinit,pmin(0:kdm),realat
      integer i,j,k,k1,l,m,n
cdiag character text*24
      character ptxt*12
c
      real     poflat
      external poflat
c
      include 'stmt_fns.h'
c
c --- -------------------------
c --- mass field initialization
c --- -------------------------
c
      margin = 0
c
      if     (iniflg.eq.2) then
        call rdrlax(mnth,1)
!$OMP   PARALLEL DO PRIVATE(j,l,i,k)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              do k=1,kk
                if (k.eq.1 .or. k.le.nhybrd) then
                  temp(i,j,k,1)=twall(i,j,k,1)
                  saln(i,j,k,1)=swall(i,j,k,1)
                  th3d(i,j,k,1)=sig(temp(i,j,k,1),saln(i,j,k,1))-thbase
                else  ! isopyc
                  temp(i,j,k,1)=tofsig(theta(k)+thbase,swall(i,j,k,1))
                  saln(i,j,k,1)=swall(i,j,k,1)
                  th3d(i,j,k,1)=theta(k)
                endif
c
                temp(i,j,k,2)=temp(i,j,k,1)
                saln(i,j,k,2)=saln(i,j,k,1)
                th3d(i,j,k,2)=th3d(i,j,k,1)
              enddo
            enddo
          enddo
        enddo
      else
!$OMP   PARALLEL DO PRIVATE(j,l,i,k)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              do k=1,kk
                temp(i,j,k,1)=tofsig(theta(k)+thbase,saln0)
                saln(i,j,k,1)=saln0
                th3d(i,j,k,1)=theta(k)
c
                temp(i,j,k,2)=temp(i,j,k,1)
                saln(i,j,k,2)=saln(i,j,k,1)
                th3d(i,j,k,2)=th3d(i,j,k,1)
              enddo
            enddo
          enddo
        enddo
      endif
c
      if     (mnproc.eq.1) then
      write (lp,'('' sigma(k):'',9f7.2/(15x,9f7.2))')
     &   (sigma(k),k=1,kk)
      write (lp,'('' theta(k):'',9f7.2/(15x,9f7.2))')
     &   (theta(k)+thbase,k=1,kk)
      endif !1st tile
      call xcsync(flush_lp)
      if     (iniflg.lt.0 .or. iniflg.gt.2) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in inicon - invalid iniflg value'
        write(lp,*) 'iniflg = ',iniflg
        write(lp,*)
        endif !1st tile
        call xcstop('(inicon)')
               stop '(inicon)'
      endif
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,pmin,realat,pinit)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 54 j=1-margin,jj+margin
      do 54 l=1,isp(j)
      do 54 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,1)=0.0
      pmin(0)=0.0
      do 55 k=1,kk
      if     (k.le.nhybrd) then
        pmin(k)=pmin(k-1)+min(dp0k(k),dp0sig(i,j))
      else  ! isopyc
        pmin(k)=pmin(k-1)
      endif
c
      if     (iniflg.le.1) then
c
c       initial interfaces from zonal mean climatology.
c
        tracer(i,j,k)=0.
        if (k.lt.kk) then
          if (iniflg.eq.0) then
c
c ---       initial interfaces are flat,
c ---       based on zonal mean climatology at center of the basin.
c
            if (mapflg.le.1) then
              realat=alat( jtdm/2-ypivn,gridn)
            else
              realat=alatu(jtdm/2-ypivn,gridn)
            endif
          else  ! iniflg==1
            if (mapflg.le.1) then
              realat=alat( (j+j0)-ypivn,gridn)
            else
              realat=alatu((j+j0)-ypivn,gridn)
            endif
          endif
          pinit=poflat(.5*(theta(k)+theta(k+1))+thbase,realat)
        else  ! k==kk
          pinit=huge
        endif
        p(i,j,k+1)=max(pmin(k),pinit)
        if     (k.gt.2 .and. k.le.nhybrd+1 .and.
     &          p(i,j,k  ).le.pmin(k-1)    .and.
     &          p(i,j,k+1).gt.pmin(k  )         ) then
          do k1=1,k-1
            th3d(i,j,k1,1)   =theta(k-1)-(k-k1)*sigjmp
            temp(i,j,k1,1)   =tofsig(th3d(i,j,k1,1)+thbase,saln0)
            saln(i,j,k1,1)   =saln0
c
            th3d(i,j,k1,2)=th3d(i,j,k1,1)
            temp(i,j,k1,2)=temp(i,j,k1,1)
            saln(i,j,k1,2)=saln(i,j,k1,1)
c
            if (tbaric) then
              thstar(i,j,k1)=th3d(i,j,k1,1)+kappaf(temp(i,j,k1,1),
     &                                             saln(i,j,k1,1),
     &                                             p(   i,j,k1))
            else
              thstar(i,j,k1)=th3d(i,j,k1,1)
            endif
          end do
        end if
        if (k.eq.kk) then
          do k1=1,kk
            p( i,j,k1+1)=min(p(i,j,k1+1),depths(i,j)*onem)
            dp(i,j,k1,1)=    p(i,j,k1+1)-p(i,j,k1)
            dp(i,j,k1,2)=   dp(i,j,k1,1)
            if (tbaric) then
              thstar(i,j,k1)=th3d(i,j,k1,1)+kappaf(temp(i,j,k1,1),
     &                                             saln(i,j,k1,1),
     &                                             p(   i,j,k1))
            else
              thstar(i,j,k1)=th3d(i,j,k1,1)
            endif
          enddo
        endif
      elseif (iniflg.eq.2) then
c
c       initial interfaces from relaxation fields.
c
        if     (k.lt.kk) then
          p(i,j,k+1) = pwall(i,j,k+1,1)
        else
          p(i,j,k+1) = depths(i,j)*onem
        endif
        dp(i,j,k,1) = p(i,j,k+1)-p(i,j,k)
        dp(i,j,k,2) = dp(i,j,k,1)
        if (tbaric) then
          thstar(i,j,k)=th3d(i,j,k,1)+
     &                  kappaf(temp(i,j,k,1),saln(i,j,k,1),p(i,j,k))
        else
          thstar(i,j,k)=th3d(i,j,k,1)
        endif
      endif
c
cdiag if (mod(k,3).ne.1) go to 55
cdiag write (text,'(''intf.pressure (m), k='',i3)') k+1
cdiag call prtmsk(ip,p(1-nbdy,1-nbdy,k+1),util1,idm,ii,jj,0.,1./onem,text)
c
 55   continue
c
      if     (isopyc) then
c
c ---   MICOM-like mixed layer no thinner than thkmin.
c
        p( i,j,2)  =max(p(i,j,2),min(depths(i,j),thkmin)*onem)
        dp(i,j,1,1)=p(i,j,2)-p(i,j,1)
        dp(i,j,1,2)=dp(i,j,1,1)
        do k=2,kk
          p( i,j,k+1)=max(p(i,j,k+1),p(i,j,k))
          dp(i,j,k,1)=    p(i,j,k+1)-p(i,j,k)
          dp(i,j,k,2)=dp(i,j,k,1)
        enddo
      endif
 54   continue
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 50 j=1-margin,jj+margin
      do 51 l=1,isp(j)
      do 51 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      pbavg(i,j,1)=0.
      pbavg(i,j,2)=0.
      pbavg(i,j,3)=0.
      pbot(i,j)=p(i,j,kk+1)
c
      montg(i,j,1)=0.
      do k=1,kk-1
        if(p(i,j,k+1).lt.thkmin*onem) then
          tracer(i,j,k)=1.
        endif
        montg(i,j,k+1)=montg(i,j,k)-
     &    p(i,j,k+1)*(thstar(i,j,k+1)-thstar(i,j,k))*thref**2
      enddo
c
      thkk(i,j)=thstar(i,j,kk)
      psikk(i,j)=montg(i,j,kk)
c
c --- start with a thin mixed layer
      if     (hybrid) then
        dpmixl(i,j,1)=min(depths(i,j)*onem-onem,
     &                    max(thkmin*onem,p(i,j,2)))
      else  ! isopyc
        dpmixl(i,j,1)=p(i,j,2)
      endif
      dpmixl(i,j,2)=dpmixl(i,j,1)
      dpbl(  i,j)  =dpmixl(i,j,1)
c
      temice(i,j) = temp(i,j,1,1)
      covice(i,j) = 0.0
      thkice(i,j) = 0.0
 51   continue
      do i=1-margin,ii+margin
        do k= 1,3
          ubavg(i,j,k) = 0.0
          vbavg(i,j,k) = 0.0
        enddo
        do k= 1,kk
          u(i,j,k,1) = 0.0
          u(i,j,k,2) = 0.0
          v(i,j,k,1) = 0.0
          v(i,j,k,2) = 0.0
        enddo
      enddo
 50   continue
!$OMP END PARALLEL DO
c
      if(mxlkrt) then
!$OMP   PARALLEL DO PRIVATE(j,l,i,k)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              do k=1,kk
                if(dpmixl(i,j,1).gt.p(i,j,k  ) .and.
     &             dpmixl(i,j,1).le.p(i,j,k+1)) then
                  t1sav(i,j,1)=temp(i,j,k,1)
                  s1sav(i,j,1)=saln(i,j,k,1)
                  tmlb( i,j,1)=temp(i,j,k,1)
                  smlb( i,j,1)=saln(i,j,k,1)
                  nmlb( i,j,1)=k
                  t1sav(i,j,2)=t1sav(i,j,1)
                  s1sav(i,j,2)=s1sav(i,j,1)
                  tmlb( i,j,2)=tmlb(i,j,1)
                  smlb( i,j,2)=smlb(i,j,1)
                  nmlb( i,j,2)=k
                end if
              enddo
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
      end if
c
      if (hybrid) then
        m=2
        n=1
            call pipe_comparall(m,n, 'inicon, step')
        call xctilr(p(      1-nbdy,1-nbdy,2  ),1,kk, 1,1, halo_ps)
        call dpudpv(p,depthu,depthv,dpu(1-nbdy,1-nbdy,1,n),
     &                              dpv(1-nbdy,1-nbdy,1,n))
            if     (lpipe) then
              do k= 1,kk
                write (ptxt,'(a9,i3)') 'dpu    k=',k
                call pipe_compare(dpu(1-nbdy,1-nbdy,1,n),iu,ptxt)
                write (ptxt,'(a9,i3)') 'dpv    k=',k
                call pipe_compare(dpv(1-nbdy,1-nbdy,1,n),iv,ptxt)
              enddo
            endif
        call hybgen(m,n)
            call pipe_comparall(m,n, 'inicn1, step')
        m=1
        n=2
        call xctilr(p(      1-nbdy,1-nbdy,2  ),1,kk, 1,1, halo_ps)
        call dpudpv(p,depthu,depthv,dpu(1-nbdy,1-nbdy,1,n),
     &                              dpv(1-nbdy,1-nbdy,1,n))
            if     (lpipe) then
              do k= 1,kk
                write (ptxt,'(a9,i3)') 'dpu    k=',k
                call pipe_compare(dpu(1-nbdy,1-nbdy,1,n),iu,ptxt)
                write (ptxt,'(a9,i3)') 'dpv    k=',k
                call pipe_compare(dpv(1-nbdy,1-nbdy,1,n),iv,ptxt)
              enddo
            endif
        call hybgen(m,n)
            call pipe_comparall(m,n, 'inicn2, step')
      endif
c
      if     (itest.gt.0 .and. jtest.gt.0) then
         write (lp,103) nstep,i0+itest,j0+jtest,
     &   '  istate:  temp    saln  thstar   thkns    dpth   montg',
     &   dpmixl(itest,jtest,1)/onem,
     &   (k,temp(itest,jtest,k,1),saln(itest,jtest,k,1),
     &   thstar(itest,jtest,k)+thbase,dp(itest,jtest,k,1)/onem,
     &   p(itest,jtest,k+1)/onem,montg(itest,jtest,k)/g,k=1,kk)
         write(lp,104) depths(itest,jtest)
 103     format (i9,2i5,a/23x,'mxl',32x,     f8.1/
     &                   (23x,i3,2f8.2,f8.2,2f8.1,f8.3))
 104     format (         23x,'bot',32x,     f8.1)
      endif !test tile
      call xcsync(flush_lp)
c
      return
      end
c
c
c> Revision history:
c>
c> Nov. 1999 - added code to initialize homogeneous values of thermodynamical
c>             variables near the surface
c> May  2000 - conversion to SI units
c> Aug. 2000 - added hybrid and isopycnic vertical coordinate options
