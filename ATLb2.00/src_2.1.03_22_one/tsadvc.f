      module mod_advem
c --- module for advem only
      contains
      subroutine advem(iord,fld,u,v,scal,scali,dt,fco,fc, util1,util2)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
      implicit none
c
      integer, intent(in)    :: iord
      real,    intent(in)    :: dt
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(inout) :: fld
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: u,v,scal,scali,fco,fc
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(out)   :: util1,util2
c
c combined monotone scheme, for details see section 3.3 (eqs. 34 to 37)
c in smolarkiewicz and clark, 1986, j.comput.phys.,67,no 2, p. 396-438
c and smolarkiewicz and grabowski, 1989, j.comput.phys.
c  iord   - order of scheme (1 or 2, 1 for simple donor cell scheme)
c  u,v    - mass fluxes satisfying continuity equation
c  scal   - spatial increments (squared)
c  scali  - inverse of scal
c  dt     - temporal increment
c  fco,fc - depth of the layer at previous and new time step
c  util1  - work array
c  util2  - work array
c
c  on return, fld's valid halo will be margin-5 wide.
c
      logical    lpipe_advem
      parameter (lpipe_advem=.false.)
c
      real, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & fmx,fmn,flp,fln,flx,fly,u1,v1,flxdiv
*     real, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
*    & tmpu1,tmpu2,tmpv1,tmpv2
c
      real    amount,fc2,fco2,flxdn,flxdp,flydn,flydp,q
      real*8  clip,vlume
ccc   real*8  bfore,after
      integer i,j,l,ia,ib,ja,jb,mbdy_a
c
      real onemu
      data onemu/0.0098/
c
      integer       itest,jtest,ittest,jttest
      common/testpt/itest,jtest,ittest,jttest
      save  /testpt/
c
      if     (margin.ge.5) then
        mbdy_a = margin  ! original margin (.ge. 5)
      else
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i2 /)')
     &      'error: advem called with margin.lt.5, margin = ',margin
        endif
        call xcstop('advem')
               stop 'advem'
      endif
ccc
ccc --- optional code for checking conservation properties
ccc      do j=1,jj
ccc        do l=1,isp(j)
ccc          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
ccc            util3(i,j) = fld(i,j)*fco(i,j)*scal(i,j)
ccc          enddo
ccc        enddo
ccc      enddo
ccc      call xcsum(bfore, util3,ip)
c
c --- compute low-order and part of antidiffusive fluxes
c
c --- rhs: u, v, fld+
c --- lhs: flx, fly, fmx, fmn
c
      margin = mbdy_a - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,q,jb,ja,ib,ia)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            u1(i,j)=.5*abs(u(i,j))*(fld(i,j)-fld(i-1,j))
            if (u(i,j).ge.0.) then
              q=fld(i-1,j)
            else
              q=fld(i  ,j)
            endif
            flx(i,j)=u(i,j)*q
          enddo
        enddo
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            v1(i,j)=.5*abs(v(i,j))*(fld(i,j)-fld(i,j-1))
            if (v(i,j).ge.0.) then
              q=fld(i,j-1)
            else
              q=fld(i,j  )
            endif
            fly(i,j)=v(i,j)*q
          enddo
        enddo
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            ia=i-1
            if (ip(ia,j).eq.0) ia=i
            ib=i+1
            if (ip(ib,j).eq.0) ib=i
            ja=j-1
            if (ip(i,ja).eq.0) ja=j
            jb=j+1
            if (ip(i,jb).eq.0) jb=j
            fmx(i,j)=max(fld(i,j),fld(ia,j),fld(ib,j),
     &                            fld(i,ja),fld(i,jb))
            fmn(i,j)=min(fld(i,j),fld(ia,j),fld(ib,j),
     &                            fld(i,ja),fld(i,jb))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
      if     (lpipe .and. lpipe_advem) then
c ---   compare two model runs.
        call pipe_compare_sym2(u1,  iu,'ad:11:u1    ',
     &                         v1,  iv,'ad:11:v1    ')
        call pipe_compare_sym2(flx, iu,'ad:11:flx   ',
     &                         fly, iv,'ad:11:fly   ')
        call pipe_compare_sym1(fmx, ip,'ad:11:fmx   ')
        call pipe_compare_sym1(fmn, ip,'ad:11:fmn   ')
      endif
c
c --- rhs: u, v, fld+
c --- lhs: flx, fly, fmx, fmn
c
      margin = mbdy_a - 1
c
      do j=1-margin,jj+margin
        do l=1,isp(j)
          if     (ifp(j,l).ge. 1-margin) then
            flx(ifp(j,l)  ,j)=0.
          endif
          if     (ilp(j,l).lt.ii+margin) then
            flx(ilp(j,l)+1,j)=0.
          endif
        enddo
      enddo
c
      do i=1-margin,ii+margin
        do l=1,jsp(i)
          if     (jfp(i,l).ge. 1-margin) then
            fly(i,jfp(i,l)  )=0.
          endif
          if     (jlp(i,l).lt.jj+margin) then
            fly(i,jlp(i,l)+1)=0.
          endif
        enddo
      enddo
      if     (lpipe .and. lpipe_advem) then
c ---   compare two model runs.
        call pipe_compare_sym2(flx, iu,'ad:22:flx   ',
     &                         fly, iv,'ad:33:fly   ')
      endif
cdiag if     (itest.gt.0 .and. jtest.gt.0) then
cdiag   i=itest
cdiag   j=jtest
cdiag   write (lp,'(a,2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3,
cdiag.             1pe9.2,0pf9.3/1pe39.2/0pf39.3)')
cdiag.    'advem (1)',i+i0,j+j0,
cdiag.    fld(i-1,j),u(i,j),fld(i,j-1),v(i,j),
cdiag.    fld(i,j),v(i,j+1),fld(i,j+1),u(i+1,j),fld(i+1,j)
cdiag endif
c
c --- rhs: flx+, fly+, fco, fmn, fmx, fc
c --- lhs: flxdiv, fld, util1, util2
c
      margin = mbdy_a - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i,q,amount)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            flxdiv(i,j)=((flx(i+1,j)-flx(i,j))+
     &                   (fly(i,j+1)-fly(i,j)))*dt*scali(i,j)
            q=fld(i,j)*fco(i,j)-flxdiv(i,j)
            amount=max(fmn(i,j)*fc(i,j),
     &                 min(q,fmx(i,j)*fc(i,j)))
            fld(i,j)=(fld(i,j)*onemu+amount)/(onemu+fc(i,j))
c
            util1(i,j)=fc(i,j)   *scal(i,j)
            util2(i,j)=(q-amount)*scal(i,j)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
c ---   compare two model runs.
        call pipe_compare_sym1(fld,    ip,'ad:610:fld  ')
        call pipe_compare_sym1(flxdiv, ip,'ad:610:flxdv')
      endif
c
c --- if iord=1, scheme reduces to simple donor cell scheme.
      if (iord.le.1) go to 100
c
c --- finish computation of antidiffusive fluxes
c
c --- rhs: u1, u, flxdiv+, fco+, fc+
c --- lhs: flx, fly
c
      margin = mbdy_a - 2
c
*     tmpu1=0.0
*     tmpu2=0.0
*     tmpv1=0.0
*     tmpv2=0.0
*
!$OMP PARALLEL DO PRIVATE(j,l,i,fc2,fco2)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        fco2=0.0
        fc2 =0.0
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            fco2=fco(i,j)+fco(i-1,j)  ! inforce order on flx calc
            fc2 =fc( i,j)+fc( i-1,j)  ! inforce order on flx calc
*           tmpu1(i,j)=fco2
*           tmpu2(i,j)=fc2
            flx(i,j)=u1(i,j)-u(i,j)*(flxdiv(i,j)+flxdiv(i-1,j))
     &         /((fco2+fc2)+onemu)
          enddo
        enddo
        if (fco2*fc2.eq.1.e30) flx(1-nbdy,j)=0.0  ! prevent removal of fc*2
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            fco2=fco(i,j)+fco(i,j-1)  ! inforce order on fly calc
            fc2 =fc( i,j)+fc( i,j-1)  ! inforce order on fly calc
*           tmpv1(i,j)=fco2
*           tmpv2(i,j)=fc2
            fly(i,j)=v1(i,j)-v(i,j)*(flxdiv(i,j)+flxdiv(i,j-1))
     &         /((fco2+fc2)+onemu)
          enddo
        enddo
        if (fco2*fc2.eq.1.e30) flx(1-nbdy,j)=0.0  ! prevent removal of fc*2
      enddo
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
c ---   compare two model runs.
*       call pipe_compare_sym2(tmpu1,iu,'ad: 8:tmpu1 ',
*    &                         tmpv1,iv,'ad: 8:tmpv1 ')
*       call pipe_compare_sym2(tmpu2,iu,'ad: 8:tmpu2 ',
*    &                         tmpv2,iv,'ad: 8:tmpv2 ')
        call pipe_compare_sym2(flx,  iu,'ad: 8:flx   ',
     &                         fly,  iv,'ad: 8:fly   ')
      endif
c
c---- limit antidiffusive fluxes
c
c --- rhs: fmx, fmn, fld, fc, flx+, fly+
c --- lhs: flp, fln
c
      margin = mbdy_a - 3
c
!$OMP PARALLEL DO PRIVATE(j,l,i,flxdn,flxdp,flydn,flydp)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            flxdp=min(0.,flx(i+1,j))-max(0.,flx(i,j))
            flxdn=max(0.,flx(i+1,j))-min(0.,flx(i,j))
            flydp=min(0.,fly(i,j+1))-max(0.,fly(i,j))
            flydn=max(0.,fly(i,j+1))-min(0.,fly(i,j))
            flp(i,j)=(fmx(i,j)-fld(i,j))*(fc(i,j)*scal(i,j))/
     &       ((onemu-(flxdp+flydp))*dt)
            fln(i,j)=(fld(i,j)-fmn(i,j))*(fc(i,j)*scal(i,j))/
     &       ((onemu+(flxdn+flydn))*dt)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
c ---   compare two model runs.
        call pipe_compare_sym1(flp, ip,'ad:16:flp   ')
        call pipe_compare_sym1(fln, ip,'ad:16:fln   ')
      endif
c
c --- rhs: flx, fly, flp+, fln+
c --- lhs: flx, fly
c
      margin = mbdy_a - 4
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            flx(i,j)=max(0.,flx(i,j))*min(1.,flp(i,j),fln(i-1,j))
     &              +min(0.,flx(i,j))*min(1.,flp(i-1,j),fln(i,j))
          enddo
        enddo
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            fly(i,j)=max(0.,fly(i,j))*min(1.,flp(i,j),fln(i,j-1))
     &              +min(0.,fly(i,j))*min(1.,flp(i,j-1),fln(i,j))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
c ---   compare two model runs.
        call pipe_compare_sym2(flx, iu,'ad:18:flx   ',
     &                         fly, iv,'ad:18:fly   ')
      endif
c
cdiag i=itest
cdiag j=jtest
cdiag write (lp,'(''advem (2)'',2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3,
cdiag.1pe9.2,0pf9.3/1pe39.2/0pf39.3)') i,j,fld(i-1,j),u(i,j),fld(i,j-1),
cdiag.v(i,j),fld(i,j),v(i,j+1),fld(i,j+1),u(i+1,j),fld(i+1,j)
c
c --- rhs: flx+, fly+, fld, fc, fmx, util2
c --- lhs: flxdiv, fld, util2
c
      margin = mbdy_a - 5
c
!$OMP PARALLEL DO PRIVATE(j,l,i,amount,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            flxdiv(i,j)=((flx(i+1,j)-flx(i,j))+
     &                   (fly(i,j+1)-fly(i,j)))*dt*scali(i,j)
            q=fld(i,j)*fc(i,j)-flxdiv(i,j)
            amount=max(fmn(i,j)*fc(i,j),min(q,fmx(i,j)*fc(i,j)))
            fld(i,j)=(fld(i,j)*onemu+amount)/(onemu+fc(i,j))
c
            util2(i,j)=util2(i,j)+(q-amount)*scal(i,j)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
c ---   compare two model runs.
        call pipe_compare_sym1(flxdiv, ip,'ad:620:flxdv')
        call pipe_compare_sym1(fld,    ip,'ad:1620:fld ')
      endif
c
  100 continue
c
c --- return 'clipped' amount to field.
      call xcsum(vlume, util1,ip)
      call xcsum(clip,  util2,ip)
c
*     if (vlume.ne.0.0 .and. clip.ne.0.0) then
      if (.false.) then
        q=clip/vlume
cdiag   write (lp,'(a,1pe11.3)') 'tracer drift in advem:',-q
c
c ---   rhs: fld
c ---   lhs: glf
c
        margin = mbdy_a - 5
c
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              fld(i,j)=fld(i,j)+q
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
        if     (lpipe .and. lpipe_advem) then
c ---     compare two model runs.
          call pipe_compare_sym1(fld, ip,'ad:13:fld   ')
        endif
      endif
ccc
ccc --- optional code for checking conservation properties
ccc      do j=1,jj
ccc        do l=1,isp(j)
ccc          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
ccc            util3(i,j) = fld(i,j)*fc(i,j)*scal(i,j)
ccc          enddo
ccc        enddo
ccc      enddo
ccc      call xcsum(after, util3,ip)
ccc      if (mnproc.eq.1) then
ccc      write (lp,'(a,1p,3e14.6,e11.1)') 'advem conservation:',
ccc     .  bfore,after,after-bfore,(after-bfore)/bfore
ccc      call flush(lp)
ccc      endif
c
c --- restore margin, but fld's valid halo is now only mbdy_a-5 wide.
c
      margin = mbdy_a
      return
      end subroutine advem
c
c> Revision history:
c>
c> Mar. 2000 - removed 'cushn' and added logic to assure global conservation
c> Apr. 2000 - conversion to SI units
c> Apr. 2000 - changed i/j loop nesting to j/i
c> Feb. 2001 - placed advem in a module
c
      end module mod_advem

      subroutine tsadvc(m,n)
      use mod_xc     ! HYCOM communication interface
      use mod_pipe   ! HYCOM debugging interface
      use mod_advem  ! defined above
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- ---------------------------------------------------
c --- thermodynamic variable(s): advection and diffusion.
c --- ---------------------------------------------------
c
      logical    lpipe_tsadvc
      parameter (lpipe_tsadvc=.false.)
c
      real, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & sold,told,q2old,q2lold
      real, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,mxtrcr) ::
     & trold
c
      logical latemp,lath3d
      integer i,isave,j,jsave,k,ktr,l,ia,ib,ja,jb,mbdy
      real sminn,smaxx,tminn,tmaxx,flxdiv
     &    ,offset,factor,pold,pmid,pnew,wts2dp
      real xmin(2*kdm),xmax(2*kdm)
c
      character*12 text,textu,textv
c
c --- select posdef:
c ---   1. as a power of 2
c ---   2. so that the ratio of the standard deviation to the mean of
c ---      each field is approximately the same:
c ---        0 for -saln-,   256 for -temp-, 32 for -th3d-,
c ---        0 for -tracer-,   1 for -q2-
      integer    pdtemp,pdth3d,pdq2
      parameter (pdtemp=256.0, pdth3d=32.0, pdq2=1.0)
c
      real harmon,a,b
      include 'stmt_fns.h'
c
c --- harmonic mean
      harmon(a,b)=2.*a*b/(a+b)
c
      mbdy = 6
c
      call xctilr(dp(     1-nbdy,1-nbdy,1,n),1,  kk, 6,6, halo_ps)
      call xctilr(dpold(  1-nbdy,1-nbdy,1  ),1,  kk, 6,6, halo_ps)
      call xctilr(saln(   1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_ps)
      call xctilr(temp(   1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_ps)
      call xctilr(th3d(   1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_ps)
      call xctilr(uflx(   1-nbdy,1-nbdy,1  ),1,  kk, 6,6, halo_uv)
      call xctilr(vflx(   1-nbdy,1-nbdy,1  ),1,  kk, 6,6, halo_vv)
      do ktr= 1,ntracr
        call xctilr(tracer( 1-nbdy,1-nbdy,1,1,ktr),1,2*kk, 6,6, halo_ps)
      enddo
      if (mxlmy) then
        call xctilr(q2(   1-nbdy,1-nbdy,0,1),1,2*kk+2, 6,6, halo_ps)
        call xctilr(q2l(  1-nbdy,1-nbdy,0,1),1,2*kk+2, 6,6, halo_ps)
      endif
c
      do 81 k=1,kk
c
c --- ---------------------------------------------------
c --- advection of thermodynamic variable(s) (and tracer)
c --- ---------------------------------------------------
c
c --- for isopycnic vertical coordinates:
c ---   advect -th3d- and -saln- in the mixed layer (layer 1),
c ---   advect            -saln- only in all other layers
c --- for hybrid vertical coordinates:
c ---   advect -temp- and -saln- in all layers if advflg==0,
c ---   advect -th3d- and -saln- in all layers if advflg==1
c
      latemp =  k.le.nhybrd .and. advflg.eq.0       ! advect temp
      lath3d = (k.le.nhybrd .and. advflg.eq.1) .or.
     &         (k.eq.1      .and. isopyc     )      ! advect th3d
c
c --- smooth mixed-layer mass fluxes in lateral direction
      if(isopyc .and. k.eq.1) then
c
c ---   rhs: vflx+, uflx+
c ---   lhs: vflux, uflux
c
        margin = mbdy - 1
c
        do j=1-margin,jj+margin
          do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              ia=max(i-1,ifv(j,l))
              ib=min(i+1,ilv(j,l))
              vflux(i,j)=.5*vflx(i,j,1)+.25*(vflx(ia,j,1)+vflx(ib,j,1))
            enddo
          enddo
        enddo
c
        do i=1-margin,ii+margin
          do l=1,jsu(i)
            do j=max(1-margin,jfu(i,l)),min(jj+margin,jlu(i,l))
              ja=max(j-1,jfu(i,l))
              jb=min(j+1,jlu(i,l))
              uflux(i,j)=.5*uflx(i,j,1)+.25*(uflx(i,ja,1)+uflx(i,jb,1))
            enddo
          enddo
        enddo
      endif
c
c ---   rhs: temp, saln, uflux+, vflux+, dpold
c ---   lhs: told, sold, util1, util2, temp, th3d
c
        margin = mbdy - 1  ! util[12] at mbdy-2
c
!$OMP PARALLEL DO PRIVATE(j,l,i,ktr,flxdiv,offset)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
c ---       save for time smoothing
            if     (latemp) then
              told(i,j)=temp(i,j,k,n)
            elseif (lath3d) then
              told(i,j)=th3d(i,j,k,n)
            endif
            sold(i,j)=saln(i,j,k,n)
            do ktr= 1,ntracr
              trold(i,j,ktr)=tracer(i,j,k,n,ktr)
              if     (trcflg(ktr).eq.2) then
               tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)+pdtemp
              endif
            enddo
            if (mxlmy) then
              q2old( i,j)=q2( i,j,k,n)
              q2lold(i,j)=q2l(i,j,k,n)
            endif
c
c --- before calling 'advem', make sure (a) mass fluxes are consistent
c --- with layer thickness change, and (b) all fields are positive-definite
            if(isopyc .and. k.eq.1) then
              flxdiv=((uflux(i+1,j)  -uflux(i,j)  )
     &               +(vflux(i,j+1)  -vflux(i,j)  ))*delt1*scp2i(i,j)
            else
              flxdiv=((uflx( i+1,j,k)-uflx( i,j,k))
     &               +(vflx( i,j+1,k)-vflx( i,j,k)))*delt1*scp2i(i,j)
            endif
            util2(i,j)=.5*(dpold(i,j,k)+dp(i,j,k,n)-flxdiv)
            util1(i,j)=.5*(dpold(i,j,k)+dp(i,j,k,n)+flxdiv)
            offset=min(0.,util1(i,j),util2(i,j))
            util2(i,j)=util2(i,j)-offset
            util1(i,j)=util1(i,j)-offset
c
            if     (latemp) then
              temp(i,j,k,n)=temp(i,j,k,n)+pdtemp
            elseif (lath3d) then
              th3d(i,j,k,n)=th3d(i,j,k,n)+pdth3d
            endif
            if     (mxlmy) then
              q2( i,j,k,n)=q2( i,j,k,n)+pdq2
              q2l(i,j,k,n)=q2l(i,j,k,n)+pdq2
            endif
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
      if     (lpipe .and. lpipe_tsadvc) then
c ---   compare two model runs.
        write(text,'(a10,i2)') '49:sold,k=',k
        call pipe_compare_sym1(sold,ip,text)
        write(text,'(a10,i2)') '49:told,k=',k
        call pipe_compare_sym1(told,ip,text)
        write(text,'(a10,i2)') '49:utl1,k=',k
        call pipe_compare_sym1(util1,ip,text)
        write(text,'(a10,i2)') '49:utl2,k=',k
        call pipe_compare_sym1(util2,ip,text)
        write (textu,'(a9,i3)') 'uflx   k=',k
        write (textv,'(a9,i3)') 'vflx   k=',k
        call pipe_compare_sym2(uflx(1-nbdy,1-nbdy,k),  iu,textu,
     &                         vflx(1-nbdy,1-nbdy,k),  iv,textv)
        write (text,'(a9,i3)') 'temp.n k=',k
        call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'saln.n k=',k
        call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'th3d.n k=',k
        call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,n),ip,text)
      endif
c
c --- rhs: temp.n, th3d.n, saln.n, uflx, vflx, uflux, vflux
c --- lhs: temp.n, th3d.n, saln.n
c
      margin = mbdy - 1  ! util[12] at mbdy-2
c
      if     (latemp) then
        call advem(2,temp(1-nbdy,1-nbdy,k,n),
     &               uflx(1-nbdy,1-nbdy,k),vflx(1-nbdy,1-nbdy,k),
     &               scp2,scp2i,delt1,util1,util2,util3,util4)
        call advem(2,saln(1-nbdy,1-nbdy,k,n),
     &               uflx(1-nbdy,1-nbdy,k),vflx(1-nbdy,1-nbdy,k),
     &               scp2,scp2i,delt1,util1,util2,util3,util4)
      elseif (lath3d .and. hybrid) then
        call advem(2,th3d(1-nbdy,1-nbdy,k,n),
     &               uflx(1-nbdy,1-nbdy,k),vflx(1-nbdy,1-nbdy,k),
     &               scp2,scp2i,delt1,util1,util2,util3,util4)
        call advem(2,saln(1-nbdy,1-nbdy,k,n),
     &               uflx(1-nbdy,1-nbdy,k),vflx(1-nbdy,1-nbdy,k),
     &               scp2,scp2i,delt1,util1,util2,util3,util4)
      elseif (lath3d .and. isopyc) then  ! MICOM-like upper layer
        call advem(2,th3d(1-nbdy,1-nbdy,k,n),uflux,vflux,
     &               scp2,scp2i,delt1,util1,util2,util3,util4)
        call advem(2,saln(1-nbdy,1-nbdy,k,n),uflux,vflux,
     &               scp2,scp2i,delt1,util1,util2,util3,util4)
      else   ! exactly isopycnal layer
        call advem(2,saln(1-nbdy,1-nbdy,k,n),
     &               uflx(1-nbdy,1-nbdy,k),vflx(1-nbdy,1-nbdy,k),
     &               scp2,scp2i,delt1,util1,util2,util3,util4)
      endif
      do ktr= 1,ntracr
*       call advem(2,tracer(1-nbdy,1-nbdy,k,n,ktr),  !MPDATA     fails
        call advem(1,tracer(1-nbdy,1-nbdy,k,n,ktr),  !donor cell works
     &               uflx(1-nbdy,1-nbdy,k),vflx(1-nbdy,1-nbdy,k),
     &               scp2,scp2i,delt1,util1,util2,util3,util4)
      enddo
      if (mxlmy) then
        call advem(2,q2( 1-nbdy,1-nbdy,k,n),
     &               uflx(1-nbdy,1-nbdy,k),vflx(1-nbdy,1-nbdy,k),
     &               scp2,scp2i,delt1,util1,util2,util3,util4)
        call advem(2,q2l(1-nbdy,1-nbdy,k,n),
     &               uflx(1-nbdy,1-nbdy,k),vflx(1-nbdy,1-nbdy,k),
     &               scp2,scp2i,delt1,util1,util2,util3,util4)
      endif
c
      if     (lpipe .and. lpipe_tsadvc) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'temp.n k=',k
        call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'saln.n k=',k
        call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'th3d.n k=',k
        call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,n),ip,text)
      endif
c
c --- rhs: temp.n, th3d.n, saln.n, dpold, dp.m, dp.n, sold, told
c --- lhs: temp.n, th3d.n, dp.m, saln.m, temp.m, th3d.m
c
      margin = mbdy - 6  ! after advem
c
      sminn=999.
      tminn=999.
      smaxx=-999.
      tmaxx=-999.
c
!$OMP PARALLEL DO PRIVATE(j,l,i,ktr,pold,pmid,pnew,wts2dp)
!$OMP&            REDUCTION(MIN:sminn,tminn) REDUCTION(MAX:smaxx,tmaxx)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            if     (latemp) then
              temp(i,j,k,n)=temp(i,j,k,n)-pdtemp
              if     (dp(i,j,k,n).gt.onemm) then
                sminn=min(sminn,saln(i,j,k,n))
                smaxx=max(smaxx,saln(i,j,k,n))
                tminn=min(tminn,temp(i,j,k,n))
                tmaxx=max(tmaxx,temp(i,j,k,n))
              endif
            elseif (lath3d) then
              th3d(i,j,k,n)=th3d(i,j,k,n)-pdth3d
              if     (dp(i,j,k,n).gt.onemm) then
                sminn=min(sminn,saln(i,j,k,n))
                smaxx=max(smaxx,saln(i,j,k,n))
                tminn=min(tminn,th3d(i,j,k,n))
                tmaxx=max(tmaxx,th3d(i,j,k,n))
              endif
            endif
            if     (mxlmy) then
              q2( i,j,k,n)=q2( i,j,k,n)-pdq2
              q2l(i,j,k,n)=q2l(i,j,k,n)-pdq2
            endif
            do ktr= 1,ntracr
              if     (trcflg(ktr).eq.2) then
               tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)-pdtemp
              endif
            enddo
c
c ---       time smoothing of thickness field
            pold=max(0.,dpold(i,j,k))
            pmid=max(0.,dp(i,j,k,m))
            pnew=max(0.,dp(i,j,k,n))
            dp(i,j,k,m)=pmid*wts1+(pold+pnew)*wts2
c ---       time smoothing of thermodynamic variables (and tracer)
            pmid=max(0.,dp(i,j,k,m))
            wts2dp=wts2/(pmid+onemm)
            saln(i,j,k,m)=saln(i,j,k,m)
     &                   +wts2dp*(pold*(sold(i,j)    -saln(i,j,k,m))+
     &                            pnew*(saln(i,j,k,n)-saln(i,j,k,m)) )
            if     (latemp) then
              temp(i,j,k,m)=temp(i,j,k,m)
     &                     +wts2dp*(pold*(told(i,j)    -temp(i,j,k,m))+
     &                              pnew*(temp(i,j,k,n)-temp(i,j,k,m)) )
              th3d(i,j,k,m)=sig(temp(i,j,k,m),saln(i,j,k,m))-thbase
            elseif (lath3d) then
              th3d(i,j,k,m)=th3d(i,j,k,m)
     &                     +wts2dp*(pold*(told(i,j)    -th3d(i,j,k,m))+
     &                              pnew*(th3d(i,j,k,n)-th3d(i,j,k,m)) )
              temp(i,j,k,m)=tofsig(th3d(i,j,k,m)+thbase,saln(i,j,k,m))
            else   ! exactly isopycnal layer
              th3d(i,j,k,m)=theta(k)
              temp(i,j,k,m)=tofsig(th3d(i,j,k,m)+thbase,saln(i,j,k,m))
            endif
            do ktr= 1,ntracr
              tracer(i,j,k,m,ktr)=tracer(i,j,k,m,ktr)
     &         +wts2dp*(pold*( trold(i,j,    ktr)-tracer(i,j,k,m,ktr))+
     &                  pnew*(tracer(i,j,k,n,ktr)-tracer(i,j,k,m,ktr)) )
            enddo
            if (mxlmy) then
              q2( i,j,k,m)=q2( i,j,k,m)
     &                 +wts2dp*(pold*(q2old( i,j)  -q2( i,j,k,m))+
     &                          pnew*(q2( i,j,k,n) -q2( i,j,k,m)) )
              q2l(i,j,k,m)=q2l(i,j,k,m)
     &                 +wts2dp*(pold*(q2lold(i,j)  -q2l(i,j,k,m))+
     &                          pnew*(q2l(i,j,k,n) -q2l(i,j,k,m)) )
            endif
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
      xmin(k)   =sminn
      xmin(k+kk)=tminn
      xmax(k)   =smaxx
      xmax(k+kk)=tmaxx
c
      if     (lpipe .and. lpipe_tsadvc) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'temp.m k=',k
        call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,m),ip,text)
        write (text,'(a9,i3)') 'temp.n k=',k
        call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'sold   k=',k
        call pipe_compare_sym1(sold,ip,text)
        write (text,'(a9,i3)') 'saln.m k=',k
        call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,m),ip,text)
        write (text,'(a9,i3)') 'saln.n k=',k
        call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'told   k=',k
        call pipe_compare_sym1(told,ip,text)
        write (text,'(a9,i3)') 'th3d.m k=',k
        call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,m),ip,text)
        write (text,'(a9,i3)') 'th3d.n k=',k
        call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,n),ip,text)
      endif
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag.write (lp,'(i9,2i5,i3,'' th,s,dp after advection  '',2f9.3,f8.2)')
cdiag.nstep,itest,jtest,k,temp(itest,jtest,k,n),saln(itest,jtest,k,n),
cdiag.dp(itest,jtest,k,n)*qonem
c
 81   continue  ! k=1,kk
*
      call pipe_comparall(m,n, 'advem,  step')
c
c --- check for negative scalar fields.
c
      call xcminr(xmin(1:2*kk))
      call xcmaxr(xmax(1:2*kk))
c
      do k= 1,kk
        sminn=xmin(k)
        tminn=xmin(k+kk)
        smaxx=xmax(k)
        tmaxx=xmax(k+kk)
c
 101    format (i9,' i,j,k =',2i5,i3,a,2f8.2)
c
        if     (latemp .and. tminn+pdtemp.lt.0.0) then
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (temp(i,j,k,n).eq.tminn) then
                  write (lp,101) nstep,i+i0,j+j0,k,
     &              ' neg. temp after advem call ',
     &              temp(i,j,k,n)
                endif
              enddo
            enddo
          enddo
          call xcsync(flush_lp)
        elseif (lath3d .and. tminn+pdth3d.lt.0.0) then
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (th3d(i,j,k,n).eq.tminn) then
                  write (lp,101) nstep,i+i0,j+j0,k,
     &              ' neg. th3d after advem call ',
     &              th3d(i,j,k,n)
                endif
              enddo
            enddo
          enddo
          call xcsync(flush_lp)
        endif
c
        if (sminn.lt.0.0) then
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (saln(i,j,k,n).eq.sminn) then
                  write (lp,101) nstep,i+i0,j+j0,k,
     &              ' neg. saln after advem call ',
     &              saln(i,j,k,n)
                endif
              enddo
            enddo
          enddo
          call xcsync(flush_lp)
        endif
c
        if (diagno) then
          if     (mnproc.eq.1) then
            write (lp,'(i9,i3, a,2f7.3, a,1pe9.2,a)')
     &        nstep,k,
     &        ' min/max of s after advection:',sminn,smaxx,
     &        '   (range:',smaxx-sminn,')'
            call flush(lp)
          endif
        endif
      enddo
c
      mbdy = 2
c
      call xctilr(saln(   1-nbdy,1-nbdy,1,n),1,kk, 2,2, halo_ps)
      call xctilr(temp(   1-nbdy,1-nbdy,1,n),1,kk, 2,2, halo_ps)
      call xctilr(th3d(   1-nbdy,1-nbdy,1,n),1,kk, 2,2, halo_ps)
      do ktr= 1,ntracr
        call xctilr(tracer( 1-nbdy,1-nbdy,1,n,ktr),1,kk, 2,2, halo_ps)
      enddo
      if (mxlmy) then
        call xctilr(q2(   1-nbdy,1-nbdy,0,n),1,kk+1, 2,2, halo_ps)
        call xctilr(q2l(  1-nbdy,1-nbdy,0,n),1,kk+1, 2,2, halo_ps)
      endif
c
      do 82 k=1,kk
c
c --- --------------------------------------
c --- diffusion of thermodynamic variable(s)
c --- --------------------------------------
c
c --- for isopycnic vertical coordinates:
c ---   diffuse -th3d- and -saln- in the mixed layer (layer 1),
c ---   diffuse            -saln- only in all other layers
c --- for hybrid vertical coordinates:
c ---   diffuse -temp- and -saln- in all layers if advflg==0,
c ---   diffuse -th3d- and -saln- in all layers if advflg==1
c
      latemp =  k.le.nhybrd .and. advflg.eq.0       ! diffus temp
      lath3d = (k.le.nhybrd .and. advflg.eq.1) .or.
     &         (k.eq.1      .and. isopyc     )      ! diffus th3d
c
c --- rhs: dp.n+, temp.n+, th3d.n+, saln.n+
c --- lhs: uflux, uflux2, vflux, vflux2
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,factor)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            factor=temdf2*aspux(i,j)*
     &             scuy(i,j)*harmon(max(dp(i-1,j,k,n),onemm)
     &                             ,max(dp(i  ,j,k,n),onemm))
            if     (latemp) then
              uflux (i,j)=factor*(temp(i-1,j,k,n)-temp(i,j,k,n))
            elseif (lath3d) then
              uflux (i,j)=factor*(th3d(i-1,j,k,n)-th3d(i,j,k,n))
            endif
            uflux2(i,j)=factor*(saln(i-1,j,k,n)-saln(i,j,k,n))
          enddo
        enddo
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            factor=temdf2*aspvy(i,j)*
     &             scvx(i,j)*harmon(max(dp(i,j-1,k,n),onemm)
     &                             ,max(dp(i,j  ,k,n),onemm))
            if     (latemp) then
              vflux (i,j)=factor*(temp(i,j-1,k,n)-temp(i,j,k,n))
            elseif (lath3d) then
              vflux (i,j)=factor*(th3d(i,j-1,k,n)-th3d(i,j,k,n))
            endif
            vflux2(i,j)=factor*(saln(i,j-1,k,n)-saln(i,j,k,n))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
        if     (lpipe .and. lpipe_tsadvc) then
c ---     compare two model runs.
          write (textu,'(a9,i3)') 'uflux  k=',k
          write (textv,'(a9,i3)') 'vflux  k=',k
          call pipe_compare_sym2(uflux, iu,textu,
     &                           vflux, iv,textv)
          write (textu,'(a9,i3)') 'uflux2 k=',k
          write (textv,'(a9,i3)') 'vflux2 k=',k
          call pipe_compare_sym2(uflux2,iu,textu,
     &                           vflux2,iv,textv)
        endif
c
c --- rhs: dp.n, uflux+, vflux+, uflux2+, vflux2+
c --- lhs: saln.n, temp.n, th3d.n
c
      margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i,factor)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            factor=-delt1/(scp2(i,j)*max(dp(i,j,k,n),onemm))
            util2(i,j)=((uflux2(i+1,j)-uflux2(i,j))
     &                 +(vflux2(i,j+1)-vflux2(i,j)))*factor
            saln(i,j,k,n)=saln(i,j,k,n)+util2(i,j)
            if     (latemp) then
              util1(i,j)=((uflux (i+1,j)-uflux (i,j))
     &                   +(vflux (i,j+1)-vflux (i,j)))*factor
              temp(i,j,k,n)=temp(i,j,k,n)+util1(i,j)
              th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
            elseif (lath3d) then
              util1(i,j)=((uflux (i+1,j)-uflux (i,j))
     &                   +(vflux (i,j+1)-vflux (i,j)))*factor
              th3d(i,j,k,n)=th3d(i,j,k,n)+util1(i,j)
              temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            else   ! exactly isopycnal layer
              th3d(i,j,k,n)=theta(k)
              temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            endif
c
cdiag       if (i.eq.itest.and.j.eq.jtest) then
cdiag         if (1.le.i .and. i.le.ii .and.
cdiag.            1.le.j .and. j.le.jj      ) then
cdiag.          write (lp,100) nstep,i+i0,j+j0,k,'t,s,dt,ds',
cdiag.          temp(i,j,k,n),saln(i,j,k,n),util1(i,j),util2(i,j)
cdiag.          call flush(lp)
cdiag.        endif
cdiag.      endif
c
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
      if     (lpipe .and. lpipe_tsadvc) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'util1  k=',k
        call pipe_compare_sym1(util1,ip,text)
        write (text,'(a9,i3)') 'util2  k=',k
        call pipe_compare_sym1(util2,ip,text)
        write (text,'(a9,i3)') 'temp.n k=',k
        call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'saln.n k=',k
        call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'th3d.n k=',k
        call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,n),ip,text)
      endif
c
cdiag if (itest.gt.0.and.jtest.gt.0) then
cdiag.  write (lp,'(i9,2i5,i3,a,2f9.3,f8.2)')
cdiag.    nstep,itest+i0,jtest+j0,k,
cdiag.    ' t,s,dp after isopyc.mix.',
cdiag.    temp(itest,jtest,k,n),saln(itest,jtest,k,n),
cdiag.    dp(itest,jtest,k,n)*qonem
cdiag.  call flush(lp)
cdiag.endif
c
c --- tracers.
c
      do ktr= 1,ntracr
c
c --- rhs: dp.n+, tracer.n+
c --- lhs: uflux, vflux
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,factor)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            factor=temdf2*aspux(i,j)*
     &             scuy(i,j)*harmon(max(dp(i-1,j,k,n),onemm)
     &                             ,max(dp(i  ,j,k,n),onemm))
            uflux(i,j)=factor*(tracer(i-1,j,k,n,ktr)-
     &                         tracer(i,  j,k,n,ktr))
          enddo
        enddo
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            factor=temdf2*aspvy(i,j)*
     &             scvx(i,j)*harmon(max(dp(i,j-1,k,n),onemm)
     &                             ,max(dp(i,j  ,k,n),onemm))
            vflux(i,j)=factor*(tracer(i,j-1,k,n,ktr)-
     &                         tracer(i,j,  k,n,ktr))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
        if     (lpipe .and. lpipe_tsadvc) then
c ---     compare two model runs.
          write (textu,'(a9,i3)') 'uflux  k=',k
          write (textv,'(a9,i3)') 'vflux  k=',k
          call pipe_compare_sym2(uflux, iu,textu,
     &                           vflux, iv,textv)
        endif
c
c --- rhs: dp.n, uflux+, vflux+
c --- lhs: tracer.n
c
      margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i,factor)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            factor=-delt1/(scp2(i,j)*max(dp(i,j,k,n),onemm))
            util1(i,j)=(uflux(i+1,j)-uflux(i,j)
     &                 +vflux(i,j+1)-vflux(i,j))*factor
            tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)+util1(i,j)
c
cdiag       if (i.eq.itest.and.j.eq.jtest) then
cdiag         if (1.le.i .and. i.le.ii .and.
cdiag.            1.le.j .and. j.le.jj      ) then
cdiag.          write (lp,100) nstep,i+i0,j+j0,k,'t,s,dt,ds',
cdiag.          tracer(i,j,k,n,ktr),0.0,util1(i,j),0.0
cdiag.          call flush(lp)
cdiag.        endif
cdiag.      endif
c
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
      if     (lpipe .and. lpipe_tsadvc) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'util1  k=',k
        call pipe_compare_sym1(util1,ip,text)
        write (text,'(a9,i3)') 'trcr.n k=',k
        call pipe_compare_sym1(tracer(1-nbdy,1-nbdy,k,n,ktr),ip,text)
      endif
c
cdiag if (itest.gt.0.and.jtest.gt.0) then
cdiag.  write (lp,'(i9,2i5,i3,a,f9.3,f8.2)')
cdiag.    nstep,itest+i0,jtest+j0,k,
cdiag.    ' trc,dp after isopyc.mix.',
cdiag.    tracer(itest,jtest,k,n,ktr),
cdiag.    dp(itest,jtest,k,n)*qonem
cdiag.  call flush(lp)
cdiag.endif
c
      enddo !ktr
c
c --- MYL2.5
c
      if (mxlmy) then
c
c --- rhs: dp.n+, temp.n+, th3d.n+, saln.n+
c --- lhs: uflux, uflux2, vflux, vflux2
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,factor)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            factor=temdf2*aspux(i,j)*
     &             scuy(i,j)*harmon(max(dp(i-1,j,k,n),onemm)
     &                             ,max(dp(i  ,j,k,n),onemm))
            uflux (i,j)=factor*(q2( i-1,j,k,n)-q2( i,j,k,n))
            uflux2(i,j)=factor*(q2l(i-1,j,k,n)-q2l(i,j,k,n))
          enddo
        enddo
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            factor=temdf2*aspvy(i,j)*
     &             scvx(i,j)*harmon(max(dp(i,j-1,k,n),onemm)
     &                             ,max(dp(i,j  ,k,n),onemm))
            vflux (i,j)=factor*(q2( i,j-1,k,n)-q2( i,j,k,n))
            vflux2(i,j)=factor*(q2l(i,j-1,k,n)-q2l(i,j,k,n))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
        if     (lpipe .and. lpipe_tsadvc) then
c ---     compare two model runs.
          write (textu,'(a9,i3)') 'uflux  k=',k
          write (textv,'(a9,i3)') 'vflux  k=',k
          call pipe_compare_sym2(uflux, iu,textu,
     &                           vflux, iv,textv)
          write (textu,'(a9,i3)') 'uflux2 k=',k
          write (textv,'(a9,i3)') 'vflux2 k=',k
          call pipe_compare_sym2(uflux2,iu,textu,
     &                           vflux2,iv,textv)
        endif
c
c --- rhs: dp.n, uflux+, vflux+, uflux2+, vflux2+
c --- lhs: saln.n, temp.n, th3d.n
c
      margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i,factor)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            factor=-delt1/(scp2(i,j)*max(dp(i,j,k,n),onemm))
            util1(i,j)=(uflux (i+1,j)-uflux (i,j)
     &                 +vflux (i,j+1)-vflux (i,j))*factor
            util2(i,j)=(uflux2(i+1,j)-uflux2(i,j)
     &                 +vflux2(i,j+1)-vflux2(i,j))*factor
            q2( i,j,k,n)=q2( i,j,k,n)+util1(i,j)
            q2l(i,j,k,n)=q2l(i,j,k,n)+util2(i,j)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
      endif !mxlmy
c
 82   continue  ! k=1,kk
c
      return
      end subroutine tsadvc
c
c  Revision history:
c
c> June 1995 - eliminated setting of salinity in massless layers (loop 46)
c>             (this is now done in mxlayr.f)
c> Aug. 1995 - omitted t/s/dp time smoothin, case of abrupt mxlayr.thk.change
c> Sep. 1995 - increased temdf2 if mixed layer occupies >90% of column
c> Aug. 2000 - temp advection and diffusion only for hybrid vertical coordinate
c> Nov. 2000 - nhybrd T&S advection layers, kdm-nhybrd S advection layers
c> Nov. 2000 - T&S or th&S advection/diffusion based on advflg
c> May  2002 - diffusion coefficent based on max(sc?x,sc?y)
