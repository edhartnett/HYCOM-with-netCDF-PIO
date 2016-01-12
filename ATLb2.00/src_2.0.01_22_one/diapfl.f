      subroutine diapf1(m,n)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
c --- KPP-style implicit interior diapycnal mixing
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --------------------
c --- diapycnal mixing
c --------------------
c
c --- interior diapycnal mixing due to three processes:
c ---   shear instability
c ---   double diffusion
c ---   background internal waves
c
c --- this is essentially the k-profile-parameterization (kpp) mixing model
c --- (mxkpp.f) with all surface boundary layer processes removed
c
c --- uses the same tri-diagonal matrix solution of vertical diffusion 
c --- equation as mxkpp.f
c
      integer j
c
      if     (mod(nstep,  mixfrq).ne.0 .and.
     .        mod(nstep+1,mixfrq).ne.0      ) then
        return  ! diapycnal mixing only every mixfrq,mixfrq+1 time steps
      endif
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call diapf1aj(m,n, j)
      enddo
!$OMP END PARALLEL DO
c
c --- momentum mixing
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call diapf1bj(m,n, j)
      enddo
!$OMP END PARALLEL DO
c
      return
      end
      subroutine diapf1aj(m,n, j)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, j
      integer i,l
c
      do l=1,isp(j)
        do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
          call diapf1aij(m,n, i,j)
        enddo
      enddo
c
      return
      end
      subroutine diapf1bj(m,n, j)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, j
      integer i,l
c
      do l=1,isu(j)
        do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
          call diapf1uij(m,n, i,j)
        enddo
      enddo
c
      do l=1,isv(j)
        do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
          call diapf1vij(m,n, i,j)
        enddo
      enddo
c
      return
      end
      subroutine diapf1aij(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
c --- KPP-style implicit interior diapycnal mixing
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j
c
c -----------------------------------------------
c --- diapycnal mixing, single i,j point (part A)
c -----------------------------------------------
c
c --- interior diapycnal mixing due to three processes:
c ---   shear instability
c ---   double diffusion
c ---   background internal waves
c
c --- this is essentially the k-profile-parameterization (kpp) mixing model
c --- (mxkpp.f) with all surface boundary layer processes removed
c
c --- uses the same tri-diagonal matrix solution of vertical diffusion 
c --- equation as mxkpp.f
c
c local variables for kpp mixing
      real shsq(kdm+1)         ! velocity shear squared
      real alfadt(kdm+1)       ! t contribution to density jump
      real betads(kdm+1)       ! s contribution to density jump
      real dbloc(kdm+1)        ! buoyancy jump across interface
      real hwide(kdm)          ! layer thicknesses in m
      real rrho                ! double diffusion parameter
      real diffdd              ! double diffusion diffusivity scale
      real prandtl             ! prandtl number
      real rigr                ! local richardson number
      real fri                 ! function of Rig
c
c --- local 1-d arrays for matrix inversion
      real t1do(kdm+1),t1dn(kdm+1),s1do(kdm+1),s1dn(kdm+1),tr1do(kdm+1),
     .     tr1dn(kdm+1),difft(kdm+1),diffs(kdm+1),
     .     zm(kdm+1),hm(kdm),dzb(kdm)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     +     tcc(kdm),         ! central ...     (k  ) ..
     +     tcl(kdm),         ! lower .....     (k-1) ..
     +     rhs(kdm)          ! right-hand-side terms
c
      real ratio,froglp
      integer k,ka,kmask,nlayer,mixflg
c
      include 'stmt_fns.h'
      froglp=.5*max(2,mixfrq)
c
c --- locate lowest substantial mass-containing layer. avoid near-zero
c --- thickness layers near the bottom
      klist(i,j)=0
      kmask=0
c
      do k=1,kk
        p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n) 
        if (dp(i,j,k,n).lt.onemm) kmask=1
        if (p(i,j,k).lt.p(i,j,kk+1)-onem.and.kmask.eq.0) klist(i,j)=k
      enddo
c
c --- calculate vertical grid and layer widths
      do k=1,kk
        if (k.eq.1) then
          hwide(k)=dp(i,j,k,n)/onem
          zgrid(i,j,k)=-.5*hwide(k)
        elseif (k.lt.klist(i,j)) then
          hwide(k)=dp(i,j,k,n)/onem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
        elseif (k.eq.klist(i,j)) then
          hwide(k)=dp(i,j,k,n)/onem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
          zgrid(i,j,k+1)=zgrid(i,j,k)-.5*hwide(k)
        else
          hwide(k)=0.
        endif
      enddo
c
c --- calculate interface variables required to estimate interior diffusivities
      do k=1,kk
        ka=min(k+1,kk)
        if (k.le.klist(i,j)) then
          shsq(  k+1)=(u(i,j,k, n)+u(i+1,j,k, n)-
     .                 u(i,j,ka,n)-u(i+1,j,ka,n))**2+
     .                (v(i,j,k, n)+v(i,j+1,k, n)-
     .                 v(i,j,ka,n)-v(i,j+1,ka,n))**2
          alfadt(k+1)=.5*(dsigdt(temp(i,j,k, n),saln(i,j,k, n))+
     .                    dsigdt(temp(i,j,ka,n),saln(i,j,ka,n)))*
     .                          (temp(i,j,k, n)-temp(i,j,ka,n))
          betads(k+1)=.5*(dsigds(temp(i,j,k, n),saln(i,j,k, n))+
     .                    dsigds(temp(i,j,ka,n),saln(i,j,ka,n)))*
     .                          (saln(i,j,k, n)-saln(i,j,ka,n))
          dbloc( k+1)=-g* thref*(th3d(i,j,k, n)-th3d(i,j,ka,n))
        endif
      enddo
c
c --- determine interior diffusivity profiles throughout the water column
c --- limit mixing to the stratified interior of the ocean
c
      do k=1,kk+1
        vcty(i,j,k)=0.
        dift(i,j,k)=0.
        difs(i,j,k)=0.
      enddo
c
c --- shear instability plus background internal wave contributions
      do k=2,kk+1
        if (k-1.le.klist(i,j) .and. p(i,j,k).gt.dpmixl(i,j,n)) then
          if (shinst) then
            rigr=max(0.,dbloc(k)*(zgrid(i,j,k-1)-zgrid(i,j,k))/
     .          (shsq(k)+epsil))
            ratio=min(rigr/rinfty,1.)
            fri=(1.-ratio*ratio)
            fri=fri*fri*fri
            vcty(i,j,k)=difm0*fri+difmiw
            difs(i,j,k)=difs0*fri+difsiw
          else
            vcty(i,j,k)=difmiw
            difs(i,j,k)=difsiw
          endif
          dift(i,j,k)=difs(i,j,k)
        endif
      enddo
c
c --- double-diffusion (salt fingering and diffusive convection)
      if (dbdiff) then
        do k=2,kk+1
          if (k-1.le.klist(i,j) .and. p(i,j,k).gt.dpmixl(i,j,n)) then
c
c --- salt fingering case
            if (-alfadt(k).gt.betads(k) .and. betads(k).gt.0.) then
              rrho= min(-alfadt(k)/betads(k),rrho0)
              diffdd=1.-((rrho-1.)/(rrho0-1.))**2
              diffdd=dsfmax*diffdd*diffdd*diffdd
              dift(i,j,k)=dift(i,j,k)+0.7*diffdd
              difs(i,j,k)=difs(i,j,k)+diffdd
c
c --- diffusive convection case
            elseif (alfadt(k).gt.0.0 .and. betads(k).lt.0.0 .and.
     .             -alfadt(k).gt.betads(k)) then
              rrho=-alfadt(k)/betads(k)
              diffdd=1.5e-6*9.*.101*exp(4.6*exp(-.54*(1./rrho-1.)))
              prandtl=.15*rrho
              if (rrho.gt..5) prandtl=(1.85-.85/rrho)*rrho
              dift(i,j,k)=dift(i,j,k)+diffdd
              difs(i,j,k)=difs(i,j,k)+prandtl*diffdd
            endif
          endif
        enddo
      endif
c
cdiag do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
cdiag if (i.eq.itest.and.j.eq.jtest) write (lp,101) (nstep,i,j,k,
cdiag. hwide(k),1.e4*vcty(i,j,k),1.e4*dift(i,j,k),
cdiag. 1.e4*difs(i,j,k),k=1,kk+1)
cdiag end do
c
c --- perform the vertical mixing at p points
c
      mixflg=0
      do k=1,klist(i,j)
        if (dift(i,j,k+1).gt.0. .or.
     .      difs(i,j,k+1).gt.0.) mixflg=mixflg+1
        difft(k+1)=froglp*dift(i,j,k+1)
        diffs(k+1)=froglp*difs(i,j,k+1)
        t1do(k)=temp(i,j,k,n)
        s1do(k)=saln(i,j,k,n)
        if (trcout) tr1do(k)=tracer(i,j,k)
        hm(k)=hwide(k)
        zm(k)=zgrid(i,j,k)
      enddo
c
      if (mixflg.le.1) return
      nlayer=klist(i,j)
      k=nlayer+1
      ka=min(k,kk)
      difft(k)=0.
      diffs(k)=0.
      t1do(k)=temp(i,j,ka,n)
      s1do(k)=saln(i,j,ka,n)
      if (trcout) tr1do(k)=tracer(i,j,ka)
      zm(k)=zgrid(i,j,k)
c
c --- compute factors for coefficients of tridiagonal matrix elements.
c       tri(k=1:NZ,0) : dt/hwide(k)/ dzb(k-1)=z(k-1)-z(k)=dzabove)
c       tri(k=1:NZ,1) : dt/hwide(k)/(dzb(k  )=z(k)-z(k+1)=dzbelow)
c
      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo
c
      tri(1,1)=delt1/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=delt1/(hm(k)*dzb(k))
        tri(k,0)=delt1/(hm(k)*dzb(k-1))
      enddo
c
c --- solve the diffusion equation
c
c --- t solution
      call tridcof(difft,tri,nlayer,tcu,tcc,tcl)
      do k=1,nlayer
        rhs(k)=t1do(k)
      enddo
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,t1do,t1dn,difft)
      do k=1,nlayer
        temp(i,j,k,n)=t1dn(k)
      enddo
c
c --- s solution and th3d reset
      call tridcof(diffs,tri,nlayer,tcu,tcc,tcl)
      do k=1,nlayer
        rhs(k)=s1do(k)
      enddo
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,s1do,s1dn,diffs)
      do k=1,nlayer
        saln(i,j,k,n)=s1dn(k)
        th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
      enddo
c
c --- tracer solution
      if (trcout) then
        do k=1,nlayer
          rhs(k)=tr1do(k)
        enddo
        call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,tr1do,tr1dn,diffs)
        do k=1,nlayer
          tracer(i,j,k)=tr1dn(k)
        enddo
      endif
c
cdiag if (i.eq.itest.and.j.eq.jtest) write (lp,102) (nstep,i,j,k,
cdiag.difft(k),t1do(k),t1dn(k),t1dn(k)-t1do(k),
cdiag.diffs(k),s1do(k),s1dn(k),s1dn(k)-s1do(k),k=1,nlayer)
c
      return
c
 101  format(25x,'   thick      viscty    t diff    s diff  '
     .     /(i9,3i4,2x,4f10.2))
 102  format(25x,
     .' diff t  t old   t new   t chng  diff s  s old   s new   s chng'
     .     /(i9,3i4,1x,8f8.3))
      end
      subroutine diapf1uij(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
c --- KPP-style implicit interior diapycnal mixing
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j 
c
c -----------------------------------------------------------------
c --- diapycnal mixing, single i,j point, momentum at u grid points
c -----------------------------------------------------------------
c
c --- local 1-d arrays for matrix inversion
      real u1do(kdm+1),u1dn(kdm+1),
     .     diffm(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     +     tcc(kdm),         ! central ...     (k  ) ..
     +     tcl(kdm),         ! lower .....     (k-1) ..
     +     rhs(kdm)          ! right-hand-side terms
c
      real presu,froglp
      integer k,ka,kmask(idm),nlayer,mixflg
c
      froglp=.5*max(2,mixfrq)
c
      presu=0.
      kmask(1)=0
      mixflg=0
      do k=1,kk+1
        ka=min(k,kk)
        if (dpu(i,j,ka,n).lt.tencm.or.k.eq.kk+1) kmask(1)=1
        if (presu.lt.depthu(i,j)-tencm.and.kmask(1).eq.0) then
          diffm(k+1)=.5*froglp*(vcty(i,j,k+1)+vcty(i-1,j,k+1))
          if (diffm(k+1).gt.0.) mixflg=mixflg+1
          u1do(k)=u(i,j,k,n)
          hm(k)=dpu(i,j,k,n)/onem
          if (k.eq.1) then
            zm(k)=-.5*hm(k)
          else
            zm(k)=zm(k-1)-.5*(hm(k-1)+hm(k))
          endif
          presu=presu+dpu(i,j,k,n)
          nlayer=k
        elseif (k.eq.nlayer+1) then
          diffm(k)=0.
          u1do(k)=u1do(k-1)
          zm(k)=zm(k-1)-.5*hm(k-1)
        endif
      enddo
      if (mixflg.le.1) return
c
c --- compute factors for coefficients of tridiagonal matrix elements.
      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo
c
      tri(1,1)=delt1/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=delt1/(hm(k)*dzb(k))
        tri(k,0)=delt1/(hm(k)*dzb(k-1))
      enddo
c
c --- solve the diffusion equation
      call tridcof(diffm,tri,nlayer,tcu,tcc,tcl)
      do k=1,nlayer
        rhs(k)= u1do(k)
      enddo
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,u1do,u1dn,diffm)
      do k=1,nlayer
        u(i,j,k,n)=u1dn(k)
      enddo
c
cdiag if (i.eq.itest.and.j.eq.jtest) write (lp,106) (nstep,i,j,k,
cdiag. hm(k),u1do(k),u1dn(k),k=1,nlayer)
c
      return
 106  format(23x,'   thick   u old   u new'/(i9,3i4,1x,f10.3,2f8.3))
      end
      subroutine diapf1vij(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
c --- KPP-style implicit interior diapycnal mixing
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j 
c
c -----------------------------------------------------------------
c --- diapycnal mixing, single i,j point, momentum at v grid points
c -----------------------------------------------------------------
c
c --- local 1-d arrays for matrix inversion
      real v1do(kdm+1),v1dn(kdm+1),
     .     diffm(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     +     tcc(kdm),         ! central ...     (k  ) ..
     +     tcl(kdm),         ! lower .....     (k-1) ..
     +     rhs(kdm)          ! right-hand-side terms
c
      real presv,froglp
      integer k,ka,kmask(idm),nlayer,mixflg
c
      froglp=.5*max(2,mixfrq)
c
      presv=0.
      kmask(1)=0
      mixflg=0
      do k=1,kk+1
        ka=min(k,kk)
        if (dpv(i,j,ka,n).lt.tencm.or.k.eq.kk+1) kmask(1)=1
        if (presv.lt.depthv(i,j)-tencm.and.kmask(1).eq.0) then
          diffm(k+1)=.5*froglp*(vcty(i,j,k+1)+vcty(i,j-1,k+1))
          if (diffm(k+1).gt.0.) mixflg=mixflg+1
          v1do(k)=v(i,j,k,n)
          hm(k)=dpv(i,j,k,n)/onem
          if (k.eq.1) then
            zm(k)=-.5*hm(k)
          else
            zm(k)=zm(k-1)-.5*(hm(k-1)+hm(k))
          endif
          presv=presv+dpv(i,j,k,n)
          nlayer=k
        elseif (k.eq.nlayer+1) then
          diffm(k)=0.
          v1do(k)=v1do(k-1)
          zm(k)=zm(k-1)-.5*hm(k-1)
        endif
      enddo
      if (mixflg.le.1) return
c
c --- compute factors for coefficients of tridiagonal matrix elements.
c
      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo
c
      tri(1,1)=delt1/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=delt1/(hm(k)*dzb(k))
        tri(k,0)=delt1/(hm(k)*dzb(k-1))
      enddo
c
c --- solve the diffusion equation
      call tridcof(diffm,tri,nlayer,tcu,tcc,tcl)
      do k=1,nlayer
        rhs(k)=v1do(k)
      enddo
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,v1do,v1dn,diffm)
      do k=1,nlayer
        v(i,j,k,n)=v1dn(k)
      enddo
c
cdiag if (i.eq.itest.and.j.eq.jtest) write (lp,107) (nstep,i,j,k,
cdiag. hm(k),v1do(k),v1dn(k),k=1,nlayer)
c
      return
 107  format(23x,'   thick   v old   v new'/(i9,3i4,1x,f10.3,2f8.3))
      end
c
      subroutine diapf2(m,n)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
c --- MICOM-style explict interior diapycnal mixing for hybrid coordinates
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
      integer j
c
      if (diapyc.eq.0. .or. (mod(nstep  ,mixfrq).ne.0 .and.
     .                       mod(nstep+1,mixfrq).ne.0)) return
cdiag write (lp,'(i9,3x,a)') nstep,'entering   d i a p f l'
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 31 j=1-margin,jj+margin
        call diapf2j(m,n, j)
   31 continue
c
      call dpudpv(p,depthu,depthv,dpu(1-nbdy,1-nbdy,1,n),
     &                            dpv(1-nbdy,1-nbdy,1,n))
c
cdiag write (lp,'(i9,3x,a)') nstep,'exiting    d i a p f l'
      return
      end
      subroutine diapf2j(m,n, j)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, j
c
      integer i,k,ka,kmin(idm),kmax(idm),l
      real flxu(idm,kdm),flxl(idm,kdm),pdot(idm,kdm),flngth(idm,kdm),
     .     ennsq,alfa,beta,q,qmin,qmax,amount,froglp,delp,
     .     trflxu(idm,0:kdm+1),trflxl(idm,0:kdm+1),
     .     tflxu(idm,0:kdm+1),tflxl(idm,0:kdm+1),clipt(idm),
     .     sflxu(idm,0:kdm+1),sflxl(idm,0:kdm+1),clips(idm),
     .     told(idm,2),sold(idm,2),trold(idm,2)
*     real totem(idm),tosal(idm),tndcyt,tndcys	! col.integrals (diag.use only)
c
      real       small
      parameter (small=1.e-6)
c
      include 'stmt_fns.h'
c
c --- -------------------------------
c --- diapycnal mixing, single j-row.
c --- -------------------------------
c
c --- if mixfrq > 1, apply mixing algorithm to both time levels
      froglp=max(2,mixfrq)
c
      do 31 l=1,isp(j)
c
c --- t/s conservation diagnostics (optional):
*     do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
*       totem(i)=0.
*       tosal(i)=0.
*       do k=1,kk
*         totem(i)=totem(i)+temp(i,j,k,n)*dp(i,j,k,n)
*         tosal(i)=tosal(i)+saln(i,j,k,n)*dp(i,j,k,n)
*       end do
*     end do
c
      do 33 k=1,kk
      do 33 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
   33 continue
c
      do 32 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      sold(i,1)=saln(i,j,kk,n)
      told(i,1)=temp(i,j,kk,n)
      trold(i,1)=tracer(i,j,kk)
      tflxl(i,   0)=0.
      tflxu(i,kk+1)=0.
      sflxl(i,   0)=0.
      sflxu(i,kk+1)=0.
c
cdiag if (i.eq.itest.and.j.eq.jtest)
cdiag.  write (lp,'(i9,2i5,3x,a/(i36,4f10.3))') nstep,
cdiag.  i,j,'before diapf2: thickness  salinity temperature density',
cdiag.  (k,dp(i,j,k,n)/onem,saln(i,j,k,n),
cdiag.  temp(i,j,k,n),th3d(i,j,k,n)+thbase,k=1,kk)
c
      kmin(i)=kk+1
      kmax(i)=1
   32 continue
c
      do 36 k=2,kk
c
      do 36 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
c --- locate lowest mass-containing layer and upper edge of stratified region
      if (p(i,j,k).lt.p(i,j,kk+1)-onemm)  then
        kmax(i)=k
        if (kmin(i).eq.kk+1 .and.
     .      th3d(i,j,k,n).gt.th3d(i,j,k-1,n)+sigjmp) then
          kmin(i)=k
        endif
      end if
   36 continue
c
cdiag if (j.eq.jtest.and.itest.ge.ifp(j,l).and.itest.le.ilp(j,l))
cdiag.  write (lp,'(i9,2i5,a,2i5)') nstep,itest,j,' kmin,kmax =',
cdiag.  kmin(itest),kmax(itest)
c
c --- find buoyancy frequency for each layer
c
      do 43 k=2,kk-1
c
      do 43 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      if (k.gt.kmin(i) .and. k.lt.kmax(i)) then
c --- ennsq = buoy.freq.^2 / g^2
        ennsq=max(0.,min(th3d(i,j,k+1,n)-th3d(i,j,k,  n),
     .                   th3d(i,j,k,  n)-th3d(i,j,k-1,n)))
     .    /max(p(i,j,k+1)-p(i,j,k),onem)
c --- store (exch.coeff x buoy.freq.^2 / g x time step) in -flngth-
c --- (dimensions of flngth: length in pressure units)
c -----------------------------------------------------------------------
c --- use the following if exch.coeff. = diapyc / buoyancy frequency
        flngth(i,k)=diapyc*sqrt(ennsq) * baclin*froglp * onem
c -----------------------------------------------------------------------
c --- use the following if exch.coeff. = diapyc
ccc        flngth(i,k)=diapyc*ennsq*g * baclin*froglp * onem
c -----------------------------------------------------------------------
c
      end if
   43 continue
c
c --- find t/s fluxes at the upper and lower interface of each layer
c --- (compute only the part common to t and s fluxes)
c
      do 37 k=1,kk
c
      do 37 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      flxu(i,k)=0.
      flxl(i,k)=0.
c
      if (k.gt.kmin(i) .and. k.lt.kmax(i)) then
c
        alfa=-thref*dsigdt(temp(i,j,k,n),saln(i,j,k,n))
        beta= thref*dsigds(temp(i,j,k,n),saln(i,j,k,n))
c
        flxu(i,k)=flngth(i,k)/
     .    max(beta*(saln(i,j,k,n)-saln(i,j,k-1,n))
     .       -alfa*(temp(i,j,k,n)-temp(i,j,k-1,n)),small)
        flxl(i,k)=flngth(i,k)/
     .    max(beta*(saln(i,j,k+1,n)-saln(i,j,k,n))
     .       -alfa*(temp(i,j,k+1,n)-temp(i,j,k,n)),small)
c
        q=min(1.,.5*min(p(i,j,k)-p(i,j,k-1),p(i,j,k+2)-p(i,j,k+1))/
     .    max(flxu(i,k),flxl(i,k),epsil))
c
cdiag   if (q.ne.1.) write (lp,'(i9,2i5,i3,a,1p,2e10.2,0p,2f7.2,f5.2)') 
cdiag.    nstep,i,j,k,' flxu/l,dpu/l,q=',flxu(i,k),flxl(i,k),
cdiag.    (p(i,j,k)-p(i,j,k-1))/onem,(p(i,j,k+2)-p(i,j,k+1))/onem,q
c
        flxu(i,k)=flxu(i,k)*q
        flxl(i,k)=flxl(i,k)*q
c
      end if				!  kmin < k < kmax
c
cdiag if (i.eq.itest.and.j.eq.jtest.and.k.ge.kmin(i).and.k.le.kmax(i))
cdiag.   write (lp,'(i9,2i5,i3,3x,a/22x,f9.3,2f7.3,1p,3e10.3)')nstep,
cdiag.   i,j,k,'thknss   temp   saln    flngth      flxu      flxl',
cdiag.   dp(i,j,k,n)/onem,temp(i,j,k,n),saln(i,j,k,n),flngth(i,k),
cdiag.   flxu(i,k)/onem,flxl(i,k)/onem
c
   37 continue
c
c --- determine mass flux -pdot- implied by t/s fluxes.
c
      do 38 k=1,kk
c
      do 38 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      pdot(i,k)=0.
      if (k.gt.kmin(i) .and. k.le.kmax(i))
     .    pdot(i,k)=flxu(i,k)-flxl(i,k-1)
   38 continue
c
c --- convert flxu,flxl into actual t/s (and tracer) fluxes
c --- cut tracer fluxes in half to compensate for execution at 2 consecutive
c --- time steps
c
      do 35 k=1,kk
c
      do 35 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      tflxu(i,k)=0.
      tflxl(i,k)=0.
      sflxu(i,k)=0.
      sflxl(i,k)=0.
      trflxu(i,k)=0.
      trflxl(i,k)=0.
      if (k.gt.kmin(i) .and. k.lt.kmax(i)) then
        tflxu(i,k)=flxu(i,k)*temp(i,j,k-1,n)
        sflxu(i,k)=flxu(i,k)*saln(i,j,k-1,n)
        trflxu(i,k)=flxu(i,k)*tracer(i,j,k-1) * .5
c
        tflxl(i,k)=flxl(i,k)*temp(i,j,k+1,n)
        sflxl(i,k)=flxl(i,k)*saln(i,j,k+1,n)
        trflxl(i,k)=flxl(i,k)*tracer(i,j,k+1) * .5
      end if
   35 continue
c
      do 34 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      clipt(i)=0.
      clips(i)=0.
   34 continue
c
c --- update interface pressure and layer temperature/salinity
      do 39 k=kk,1,-1
      ka=max(1,k-1)
c
      do 39 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      sold(i,2)=sold(i,1)
      sold(i,1)=saln(i,j,k,n)
      told(i,2)=told(i,1)
      told(i,1)=temp(i,j,k,n)
      trold(i,2)=trold(i,1)
      trold(i,1)=tracer(i,j,k)
c
      dpold(i,j,k)=dp(i,j,k,n)
      p(i,j,k)=p(i,j,k)-pdot(i,k)
      dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
c
      if (k.ge.kmin(i) .and. k.le.kmax(i)) then
        delp=dp(i,j,k,n)
        if (delp.gt.0.) then
        amount=temp(i,j,k,n)*dpold(i,j,k)
     .    -(tflxu(i,k+1)-tflxu(i,k)+tflxl(i,k-1)-tflxl(i,k))
        q=amount
        qmax=max(temp(i,j,ka,n),told(i,1),told(i,2))
        qmin=min(temp(i,j,ka,n),told(i,1),told(i,2))
        amount=max(qmin*delp,min(amount,qmax*delp))
        clipt(i)=clipt(i)+(q-amount)
        temp(i,j,k,n)=amount/delp
c
        amount=saln(i,j,k,n)*dpold(i,j,k)
     .    -(sflxu(i,k+1)-sflxu(i,k)+sflxl(i,k-1)-sflxl(i,k))
        q=amount
        qmax=max(saln(i,j,ka,n),sold(i,1),sold(i,2))
        qmin=min(saln(i,j,ka,n),sold(i,1),sold(i,2))
        amount=max(qmin*delp,min(amount,qmax*delp))
        clips(i)=clips(i)+(q-amount)
        saln(i,j,k,n)=amount/delp
c
        amount=tracer(i,j,k)*dpold(i,j,k)
     .    -(trflxu(i,k+1)-trflxu(i,k)+trflxl(i,k-1)-trflxl(i,k))
        q=amount
        qmax=max(tracer(i,j,ka),trold(i,1),trold(i,2))
        qmin=min(tracer(i,j,ka),trold(i,1),trold(i,2))
        amount=max(qmin*delp,min(amount,qmax*delp))
        tracer(i,j,k)=amount/delp
        end if
      end if
   39 continue
c
      do 30 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      clipt(i)=clipt(i)/pbot(i,j)
      clips(i)=clips(i)/pbot(i,j)
   30 continue
c
      do 41 k=1,kk
      do 41 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
c --- restore 'clipped' t/s amount to column
      temp(i,j,k,n)=temp(i,j,k,n)+clipt(i)
      saln(i,j,k,n)=saln(i,j,k,n)+clips(i)
      th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
c
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,k,n)-dpold(i,j,k))	! diapyc.flx.
c --- make sure p is computed from dp, not the other way around (roundoff!)
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
   41 continue
c
c --- t/s conservation diagnostics (optional):
*     do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
*       tndcyt=-totem(i)
*       tndcys=-tosal(i)
*       do k=1,kk
*         tndcyt=tndcyt+temp(i,j,k,n)*dp(i,j,k,n)
*         tndcys=tndcys+saln(i,j,k,n)*dp(i,j,k,n)
*       end do
*       if (abs(tndcyt/totem(i)).gt.1.e-11)
*    .  write (lp,100) i,j,'  diapf2 temp.col.intgl.:',totem(i),tndcyt,
*    .  clipt(i)
*       if (abs(tndcys/tosal(i)).gt.1.e-11)
*    .  write (lp,100) i,j,'  diapf2 saln.col.intgl.:',tosal(i),tndcys,
*    .  clips(i)
*100    format(2i4,a,1p,e16.8,2e13.5)
*     end do
c
cdiag do 31 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
cdiag if (i.eq.itest.and.j.eq.jtest)
cdiag.  write (lp,'(i9,2i5,3x,a/(i36,0p,4f10.3))') nstep,
cdiag.  i,j,'after  diapf2: thickness  salinity temperature density',
cdiag.  (k,dp(i,j,k,n)/onem,saln(i,j,k,n),
cdiag.  temp(i,j,k,n),th3d(i,j,k,n)+thbase,k=1,kk)
c
   31 continue
c
      return
      end
c
      subroutine diapf3(m,n)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0 (adapted from micom version 2.8)
c --- MICOM-style explict interior diapycnal mixing for isopycnal coords
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
      integer j
c
      if (diapyc.eq.0. .or. (mod(nstep  ,mixfrq).ne.0 .and.
     .                       mod(nstep+1,mixfrq).ne.0)) return
cdiag write (lp,'(i9,3x,a)') nstep,'entering   d i a p f l'
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 31 j=1-margin,jj+margin
        call diapf3j(m,n, j)
   31 continue
c
      call dpudpv(p,depthu,depthv,dpu(1-nbdy,1-nbdy,1,n),
     &                            dpv(1-nbdy,1-nbdy,1,n))
c
cdiag write (lp,'(i9,3x,a)') nstep,'exiting    d i a p f l'
      return
      end
      subroutine diapf3j(m,n, j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0 (adapted from micom version 2.8)
c --- MICOM-style explict interior diapycnal mixing for isopycnic coordinates
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n,j
c
      integer i,k,l
      integer kmin(idm),kmax(idm)
      real flxu(idm,kdm),flxl(idm,kdm),pdot(idm,kdm),flngth(idm,kdm),
     .     ennsq,alfa,beta,smax,smin,sold(idm,2),q,salt,froglp,
     .     trflx(idm,kdm+1),trold(idm,2),trmax,trmin
c
      real       small
      parameter (small=1.e-6)
c
      include 'stmt_fns.h'
c
c --- ------------------------------
c --- diapycnal mixing, single j-row.
c --- -------------------------------
c
c --- if mixfrq > 1, apply mixing algorithm to both time levels
      froglp=max(2,mixfrq)
c
ccc      salt=0.
c
      do 31 l=1,isp(j)
c
      do 33 k=1,kk
      do 33 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      trflx(i,k+1)=0.
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
 33   continue
c
      do 32 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      sold(i,1)=saln(i,j,kk,n)
      trold(i,1)=tracer(i,j,kk)
      flxl(i,1)=0.
c
cdiag if (i.eq.itest .and. j.eq.jtest)
cdiag.  write (lp,'(i9,2i5,3x,a/(i36,0p,3f10.3,3p,f10.3))') nstep,
cdiag.  i,j,'before diapfl: thickness  salinity temperature density',
cdiag.  1,dp(i,j,1,n)/onem,saln(i,j,1,n),temp(i,j,1,n),
cdiag.  th3d(i,j,1,n)+thbase,(k,dp(i,j,k,n)/onem,saln(i,j,k,n),
cdiag.  temp(i,j,k,n),th3d(i,j,k,n)+thbase,k=2,kk)
c
      kmin(i)=kk+1
 32   kmax(i)=1
c
      do 36 k=2,kk
c
      do 36 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
c --- locate lowest mass-containing layer
      if (p(i,j,k).lt.p(i,j,kk+1))  then
        kmax(i)=k
c
c --- make sure salinity is compatible with density in layer k
ccc     q=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
ccc     salt=salt+(q-saln(i,j,k,n))*dp(i,j,k,n)*scp2(i)
ccc     saln(i,j,k,n)=q
c
c --- locate uppermost isopycnic layer heavier than mixed layer
        if (kmin(i).eq.kk+1 .and.
     .    max(th3d(i,j,1,m),th3d(i,j,1,n))+sigjmp.le.th3d(i,j,k,n))
     .    kmin(i)=k
      end if
 36   continue
c
cdiag if (j.eq.jtest.and.itest.ge.ifp(j,l).and.itest.le.ilp(j,l))
cdiag.  write (lp,'(i9,2i5,a,2i5)') nstep,itest,j,' kmin,kmax =',
cdiag.  kmin(itest),kmax(itest)
c
c --- temporarily swap layers  1  and  kmin-1
c
      do 40 k=3,kk
c
      do 40 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      if (k.eq.kmin(i)) then
        q=temp(i,j,k-1,n)
        temp(i,j,k-1,n)=temp(i,j,1,n)
        temp(i,j,1,n)=q
c
        q=saln(i,j,k-1,n)
        saln(i,j,k-1,n)=saln(i,j,1,n)
        saln(i,j,1,n)=q
c
        q=th3d(i,j,k-1,n)
        th3d(i,j,k-1,n)=th3d(i,j,1,n)
        th3d(i,j,1,n)=q
c
        flxl(i,k-1)=0.
      end if
 40   continue
c
c --- find buoyancy frequency for each layer
c
      do 43 k=2,kk
c
      do 43 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      if (k.ge.kmin(i) .and. k.lt.kmax(i)) then
c --- ennsq = buoy.freq.^2 / g^2
        ennsq=max(0.,min(th3d(i,j,k+1,n)-th3d(i,j,k,  n),
     .                   th3d(i,j,k,  n)-th3d(i,j,k-1,n)))
     .        /max(p(i,j,k+1)-p(i,j,k),onem)
c --- store (exch.coeff x buoy.freq.^2 / g x time step) in -flngth-
c -----------------------------------------------------------------------
c --- use the following if exch.coeff. = diapyc / buoyancy frequency
        flngth(i,k)=(diapyc*sqrt(ennsq) * baclin*froglp * onem)
c -----------------------------------------------------------------------
c --- use the following if exch.coeff. = diapyc
ccc        flngth(i,k)=diapyc*ennsq*g * baclin*froglp * onem
c -----------------------------------------------------------------------
c
      end if
 43   continue
c
c --- find t/s fluxes at the upper and lower interface of each layer
c --- (compute only the part common to t and s fluxes)
c
      do 37 k=2,kk
c
      do 37 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      flxu(i,k)=0.
      flxl(i,k)=0.
c
      if (k.ge.kmin(i) .and. k.lt.kmax(i)) then
c
        alfa=-thref*dsigdt(temp(i,j,k,n),saln(i,j,k,n))
        beta= thref*dsigds(temp(i,j,k,n),saln(i,j,k,n))
c
        flxu(i,k)=flngth(i,k)/
     .    max(beta*(saln(i,j,k,n)-saln(i,j,k-1,n))
     .       -alfa*(temp(i,j,k,n)-temp(i,j,k-1,n)),small)
        flxl(i,k)=flngth(i,k)/
     .    max(beta*(saln(i,j,k+1,n)-saln(i,j,k,n))
     .       -alfa*(temp(i,j,k+1,n)-temp(i,j,k,n)),small)
c
        q=min(1.,min(p(i,j,k)-p(i,j,k-1),p(i,j,k+2)-p(i,j,k+1))/
     .    max(flxu(i,k),flxl(i,k),epsil))
c
cdiag   if (i.eq.itest .and. j.eq.jtest .and. q.ne.1.)
cdiag.    write (lp,'(i9,2i5,i3,a,1p,2e10.2,0p,2f7.2,f4.2)') 
cdiag.    nstep,i,j,k,' flxu/l,dpu/l,q=',flxu(i,k),flxl(i,k),
cdiag.    (p(i,j,k)-p(i,j,k-1))/onem,(p(i,j,k+2)-p(i,j,k+1))/onem,q
c
        flxu(i,k)=flxu(i,k)*q
        flxl(i,k)=flxl(i,k)*q
c
      end if				!  kmin < k < kmax
c
cdiag if (i.eq.itest.and.j.eq.jtest.and.k.ge.kmin(i).and.k.le.kmax(i))
cdiag.   write (lp,'(i9,2i5,i3,3x,a/22x,f9.3,2f7.3,1p,3e10.3)')nstep,
cdiag.   i,j,k,' thknss   temp   saln   flngth      flxu      flxl',
cdiag.   dp(i,j,k,n)/onem,temp(i,j,k,n),saln(i,j,k,n),flngth(i,k),
cdiag.   flxu(i,k)/onem,flxl(i,k)/onem
c
 37   continue
c
c --- determine mass flux -pdot- implied by t/s fluxes.
c
      do 38 k=2,kk
c
      do 38 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      if (k.ge.kmin(i) .and. k.le.kmax(i))
     .    pdot(i,k)=flxu(i,k)-flxl(i,k-1)
 38   continue
c
c --- convert flxu,flxl into actual salt (and tracer) fluxes
c
      do 35 k=2,kk
c
      do 35 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      if (k.ge.kmin(i) .and. k.lt.kmax(i)) then
        flxu(i,k)=-flxu(i,k)*(saln(i,j,k,n)-saln(i,j,k-1,n))
        flxl(i,k)=-flxl(i,k)*(saln(i,j,k+1,n)-saln(i,j,k,n))
c --- upstream fluxes for tracer.
c --- cut fluxes in half to compensate for execution at 2 consecutive time steps
        trflx(i,k+1)=
     .     .25*(tracer(i,j,k  )*(pdot(i,k+1)+abs(pdot(i,k+1)))
     .         +tracer(i,j,k+1)*(pdot(i,k+1)-abs(pdot(i,k+1))))
      end if
 35   continue
c
c --- update interface pressure and layer salinity
      do 39 k=kk,2,-1
c
      do 39 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      sold(i,2)=sold(i,1)
      sold(i,1)=saln(i,j,k,n)
      trold(i,2)=trold(i,1)
      trold(i,1)=tracer(i,j,k)
c
      if (k.ge.kmin(i) .and. k.le.kmax(i)) then
        p(i,j,k)=p(i,j,k)-pdot(i,k)
        saln(i,j,k,n)=saln(i,j,k,n)-(flxl(i,k)-flxu(i,k))
     .      /max(p(i,j,k+1)-p(i,j,k),small)
        if(trcout) tracer(i,j,k)=tracer(i,j,k)-(trflx(i,k+1)-trflx(i,k))
     .      /max(p(i,j,k+1)-p(i,j,k),small)
c
c --- avoid excessive salinity (and tracer) values in totally eroded layers
        smin=min(saln(i,j,k-1,n),sold(i,1),sold(i,2))
        smax=max(saln(i,j,k-1,n),sold(i,1),sold(i,2))
        saln(i,j,k,n)=max(smin,min(saln(i,j,k,n),smax))
c
        if(trcout) then
          trmin=min(tracer(i,j,k-1),trold(i,1),trold(i,2))
          trmax=max(tracer(i,j,k-1),trold(i,1),trold(i,2))
          tracer(i,j,k)=max(trmin,min(tracer(i,j,k),trmax))
        end if
c
        temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
      else if (kmin(i).ne.kk+1 .and. k.lt.kmin(i)) then
        p(i,j,k)=min(p(i,j,k),p(i,j,k+1))
      end if
 39   continue
c
c --- undo effect of loop 40 (i.e., restore layers  1  and  kmin-1)
c
      do 44 k=3,kk
c
      do 44 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      if (k.eq.kmin(i)) then
        q=temp(i,j,k-1,n)
        temp(i,j,k-1,n)=temp(i,j,1,n)
        temp(i,j,1,n)=q
c
        q=saln(i,j,k-1,n)
        saln(i,j,k-1,n)=saln(i,j,1,n)
        saln(i,j,1,n)=q
c
        q=th3d(i,j,k-1,n)
        th3d(i,j,k-1,n)=th3d(i,j,1,n)
        th3d(i,j,1,n)=q
      end if
 44   continue
c
      do 42 k=kk,2,-1
      do 42 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
c --- if layer 1 has been totally eroded, transfer layer -kmin- to layer 1
      if (k.eq.kmin(i) .and. p(i,j,k).lt..1*onemm) then
c
cdiag   if (i.eq.itest .and. j.eq.jtest) then
cdiag   write (lp,'(i9,2i5,3x,a,i3,a)') nstep,i,j,'diapfl -- layer',k,
cdiag.  ' erodes mixed layer'
cdiag   write (lp,'(i9,2i5,i3,3x,a/22x,f9.3,2f7.3,1p,3e10.3)') nstep,
cdiag.   i,j,k,' thknss   temp   saln   flngth      flxu      flxl',
cdiag.   (p(i,j,k+1)-p(i,j,k))/onem,temp(i,j,k,n),
cdiag.   saln(i,j,k,n),flngth(i,k),flxu(i,k),flxl(i,k)
cdiag   end if
c
        kmin(i)=kmin(i)+1
        temp(i,j,1,n)=temp(i,j,k,n)
        saln(i,j,1,n)=saln(i,j,k,n)
        th3d(i,j,1,n)=th3d(i,j,k,n)
      end if
 42   continue
c
      do 41 k=1,kk
      do 41 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,k+1)=max(p(i,j,k),p(i,j,k+1))
      dpold(i,j,k)=dp(i,j,k,n)					! diapyc.flx.
      dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,k,n)-dpold(i,j,k))	! diapyc.flx.
c --- make sure p is computed from dp, not the other way around (roundoff!)
 41   p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
c
      do 31 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      dpmixl(i,j,n)=dp(i,j,1,n)
cdiag if (i.eq.itest.and.j.eq.jtest)
cdiag.  write (lp,'(i9,2i5,3x,a/(i36,0p,3f10.3,3p,f10.3))') nstep,
cdiag.  i,j,'after  diapfl: thickness  salinity temperature density',
cdiag.  1,dp(i,j,1,n)/onem,saln(i,j,1,n),temp(i,j,1,n),
cdiag.  th3d(i,j,1,n)+thbase,(k,dp(i,j,k,n)/onem,saln(i,j,k,n),
cdiag.  temp(i,j,k,n),th3d(i,j,k,n)+thbase,k=2,kk)
c
 31   continue
c
cdiag write (lp,'(i9,7x,1p,e9.2,a)') nstep,salt*1.e-6/g,
cdiag.  ' kg salt added in diapfl'
c
      call dpudpv(p,depthu,depthv,dpu(1-nbdy,1-nbdy,1,n),
     &                            dpv(1-nbdy,1-nbdy,1,n))
c
ccc      write (lp,'(i9,3x,a)') nstep,'exiting    d i a p f l'
      return
      end
c
c> Revision history:
c>
c> Mar. 2000 - conversion to SI units
c> May  2000 - converted T/S advection equations to flux form
c> Jul. 2000 - added diapf2j for OpenMP parallelization.
c> Aug. 2000 - adapted diapf3 from micom 2.8 to run within hycom 1.0
