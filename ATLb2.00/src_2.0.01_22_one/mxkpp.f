      subroutine mxkpp(m,n)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c ----------------------------------------------------
c --- large, mc williams, doney kpp vertical diffusion
c ----------------------------------------------------
c
      real    delp,sigmlj
      integer i,j,k,l
c
      include 'stmt_fns.h'
c
      call xctilr(u(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
      call xctilr(v(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
      call xctilr(p(      1-nbdy,1-nbdy,2  ),1,kk, 1,1, halo_ps)
      call xctilr(dpmixl( 1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_ps)
c
      margin = 0  ! no horizontal derivatives
c
c --- diffisuvity/viscosity calculation
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call mxkppaj(m,n, j)
      enddo
!$OMP END PARALLEL DO
c
c --- optional spatial smoothing of viscosity and diffusivities on interior
c --- interfaces.
c
      if(difsmo) then
        call xctilr(vcty(   1-nbdy,1-nbdy,1  ),1, 1, 1,1, halo_ps)
        do k=2,kk
          call psmooth_max(vcty(1-nbdy,1-nbdy,k), 1)
          call psmooth_max(dift(1-nbdy,1-nbdy,k), 0)
          call psmooth_max(difs(1-nbdy,1-nbdy,k), 0)
        enddo
      else
        call xctilr(vcty(   1-nbdy,1-nbdy,1  ),1,kk, 1,1, halo_ps)
      endif
c
c ---   final mixing of variables at p points
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call mxkppbj(m,n, j)
      enddo
!$OMP END PARALLEL DO
c
c --- final velocity mixing at u,v points
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call mxkppcj(m,n, j)
      enddo
!$OMP END PARALLEL DO
c
c --- mixed layer diagnostics
c
      if (diagno) then
c
c --- diagnose new mixed layer depth based on density jump criterion
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,sigmlj)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
c
c --- depth of mixed layer base set to interpolated depth where
c --- the density jump is equivalent to a tmljmp temperature jump.
c --- this may not vectorize, but is used infrequently.
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              sigmlj = -tmljmp*dsigdt(temp(i,j,1,n),saln(i,j,1,n))
              dpmixl(i,j,n)=-zgrid(i,j,klist(i,j))*onem
              do k=2,klist(i,j)
                if ((th3d(i,j,k,n)-th3d(i,j,1,n)).ge.sigmlj) then
                  dpmixl(i,j,n)=max(dp(i,j,1,n),
     .                              onem*(-zgrid(i,j,k-1)+
     .               ((zgrid(i,j,k-1)-zgrid(i,j,k))*
     .               (th3d(i,j,1,n)+sigmlj-th3d(i,j,k-1,n)))/
     .               (th3d(i,j,k,n) +epsil -th3d(i,j,k-1,n))))
                  go to 1
                endif
              enddo
 1            continue
            enddo
          enddo
        enddo
c
!$OMP   END PARALLEL DO
c
c --- calculate bulk mixed layer t, s, theta
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,delp)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
c
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              dpmixl(i,j,m)=dpmixl(i,j,n)
              tmix(i,j)=temp(i,j,1,n)*dp(i,j,1,n)
              smix(i,j)=saln(i,j,1,n)*dp(i,j,1,n)
            enddo
c
            do k=2,kk
              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                delp=min(p(i,j,k+1),dpmixl(i,j,n))
     .              -min(p(i,j,k  ),dpmixl(i,j,n))
                tmix(i,j)=tmix(i,j)+delp*temp(i,j,k,n)
               smix(i,j)=smix(i,j)+delp*saln(i,j,k,n)
              enddo
            enddo
c
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              tmix(i,j)=tmix(i,j)/dpmixl(i,j,n)
              smix(i,j)=smix(i,j)/dpmixl(i,j,n)
              thmix(i,j)=sig(tmix(i,j),smix(i,j))-thbase
           enddo
c
          enddo
        enddo
c
!$OMP   END PARALLEL DO
c
c --- calculate bulk mixed layer u
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,delp)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isu(j)
c
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
              umix(i,j)=u(i,j,1,n)*2.*dpu(i,j,1,n)
            enddo
c
            do k=2,kk
              do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
                delp=
     .             (min(p(i,j,k+1)+p(i-1,j,k+1),
     .                  dpmixl(i,j,n)+dpmixl(i-1,j,n))
     .             -min(p(i,j,k  )+p(i-1,j,k  ),
     .                  dpmixl(i,j,n)+dpmixl(i-1,j,n)))
                umix(i,j)=umix(i,j)+delp*u(i,j,k,n)
              enddo
            enddo
c
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
              umix(i,j)=umix(i,j)/(dpmixl(i,j,n)+dpmixl(i-1,j,n))
            enddo
c
          enddo
        enddo
!$OMP   END PARALLEL DO
c
c --- calculate bulk mixed layer v
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,delp)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isv(j)
c
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              vmix(i,j)=v(i,j,1,n)*2.*dpv(i,j,1,n)
            enddo
c
            do k=2,kk
              do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
                delp=
     .             (min(p(i,j,k+1)+p(i,j-1,k+1),
     .                  dpmixl(i,j,n)+dpmixl(i,j-1,n))
     .             -min(p(i,j,k  )+p(i,j-1,k  ),
     .                  dpmixl(i,j,n)+dpmixl(i,j-1,n)))
                vmix(i,j)=vmix(i,j)+delp*v(i,j,k,n)
              enddo
            enddo
c
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              vmix(i,j)=vmix(i,j)/(dpmixl(i,j,n)+dpmixl(i,j-1,n))
            enddo
c
          enddo
        enddo
!$OMP   END PARALLEL DO
      endif                                           ! diagno
c
      return
      end
      subroutine mxkppaj(m,n, j)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, j
c
c --- calculate viscosity and diffusivity
c
      integer i,l
c
      do l=1,isp(j)
        do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
          call mxkppaij(m,n, i,j)
        enddo
      enddo
c
      return
      end
      subroutine mxkppbj(m,n, j)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, j
c
c --- final mixing at p points
c
      integer i,l
c
      do l=1,isp(j)
        do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
          call mxkppbij(m,n, i,j)
        enddo
      enddo
c
      return
      end
      subroutine mxkppcj(m,n, j)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, j
c
c --- final velocity mixing at u,v points
c
      integer i,l
c
      do l=1,isu(j)
        do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
          call mxkppciju(m,n, i,j)
        enddo
      enddo
c
      do l=1,isv(j)
        do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
          call mxkppcijv(m,n, i,j)
        enddo
      enddo
c
      return
      end
      subroutine mxkppaij(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j
c
c -------------------------------------------------
c --- kpp vertical diffusion, single j-row (part A)
c -------------------------------------------------
c
c local variables for kpp mixing
      real delta               ! fraction hbl lies beteen zgrid neighbors
      real zref                ! nearsurface reference z
      real uref                ! nearsurface reference u
      real vref                ! nearsurface reference v
      real bref                ! nearsurface reference buoyancy
      real swfrac(kdm+1)       ! fractional surface shortwave radiation flux
      real shsq(kdm+1)         ! velocity shear squared
      real alfadt(kdm+1)       ! t contribution to density jump
      real betads(kdm+1)       ! s contribution to density jump
      real swfrml              ! fractional surface sw rad flux at ml base
      real ritop(kdm)          ! numerator of bulk richardson number
      real dbloc(kdm+1)        ! buoyancy jump across interface 
      real dvsq(kdm)           ! squared current shear for bulk richardson no.
      real hwide(kdm)          ! layer thicknesses in m
      real case                ! 1 in case A; =0 in case B
      real hbl                 ! boundary layer depth
      real rib(2)              ! bulk richardson number
      real rrho                ! double diffusion parameter
      real diffdd              ! double diffusion diffusivity scale
      real prandtl             ! prandtl number
      real rigr                ! local richardson number
      real fri                 ! function of Rig
      real stable              ! = 1 in stable forcing; =0 in unstable
      real dkm1(3)             ! boundary layer diffusions at nbl-1 level
      real gat1(3)             ! shape functions at dnorm=1
      real dat1(3)             ! derivative of shape functions at dnorm=1
      real blmc(kdm+1,3)       ! boundary layer mixing coefficients
      real wm                  ! momentum velocity scale
      real ws                  ! scalar velocity scale
      real dnorm               ! normalized depth
      real bfsfc               ! surface buoyancy forcing
c
      integer nbl              ! layer containing boundary layer base
      integer k1,k2            ! bulk richardson number indices
c
c --- local 1-d arrays for matrix solution
      real u1do(kdm+1),u1dn(kdm+1),v1do(kdm+1),v1dn(kdm+1),t1do(kdm+1),
     .     t1dn(kdm+1),s1do(kdm+1),s1dn(kdm+1),tr1do(kdm+1),
     .     tr1dn(kdm+1),diffm(kdm+1),difft(kdm+1),diffs(kdm+1),
     .     ghat(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
c
c --- local 1-d arrays for iteration loops
      real uold(kdm+1),vold (kdm+1),told (kdm+1),
     .     sold(kdm+1),thold(kdm+1)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     +     tcc(kdm),         ! central ...     (k  ) ..
     +     tcl(kdm),         ! lower .....     (k-1) ..
     +     rhs(kdm)          ! right-hand-side terms
c
      real dtemp,dsaln,wz,del,ratio,q,ghatflux,
     .     dvdzup,dvdzdn,viscp,difsp,diftp,f1,sigg,aa1,aa2,aa3,gm,gs,gt,
     .     dkmp2,dstar,hmonob,hblmin,hblmax,sflux1,bvsq,vtsq,
     .     vctyh,difsh,difth
c
      integer k,ka,kb,nlayer,ksave,iter,kmask,jrlv
c
      include 'stmt_fns.h'
c
c --- locate lowest substantial mass-containing layer.
      klist(i,j)=0
      kmask=0
c
      do k=1,kk
        p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
        if (dp(i,j,k,n).lt.tencm) kmask=1
        if (p(i,j,k).le.p(i,j,kk+1)-tencm.and.kmask.eq.0)
     .      klist(i,j)=k
      enddo
c
c --- forcing of t,s by surface fluxes. shortwave flux penetration depends on
c --- jerlov water type, which is set in thermf.f.  flux positive into ocean
c
      jrlv=jerlov(i,j)
c
      do k=1,kk
        if (k.eq.1) then
          if (p(i,j,k+1).lt.tenm) then
            swfrac(k+1)=redfac(jrlv) *exp(-p(i,j,k+1)/betard(jrlv))+
     .              (1.-redfac(jrlv))*exp(-p(i,j,k+1)/betabl(jrlv))
          else
            swfrac(k+1)=(1.-redfac(jrlv))*
     .                  exp(-p(i,j,k+1)/betabl(jrlv))
          endif
          sflux1=surflx(i,j)-sswflx(i,j)
          dtemp=(sflux1+(1.-swfrac(k+1))*sswflx(i,j))*delt1*g/
     .          (spcifh*max(onemm,dp(i,j,k,n)))
          dsaln=salflx(i,j)                         *delt1*g/
     .                 (max(onemm,dp(i,j,k,n)))
cdiag     if (i.eq.itest.and.j.eq.jtest) then
cdiag       write (lp,101) nstep,i+i0,j+j0,k, 
cdiag.        0.,1.-swfrac(k+1),dtemp,dsaln
cdiag       call flush(lp)
cdiag     endif
        elseif (k.le.klist(i,j)) then
          if (p(i,j,k+1).lt.tenm) then
            swfrac(k+1)=redfac(jrlv) *exp(-p(i,j,k+1)/betard(jrlv))+
     .              (1.-redfac(jrlv))*exp(-p(i,j,k+1)/betabl(jrlv))
          else
            swfrac(k+1)=(1.-redfac(jrlv))*
     .                  exp(-p(i,j,k+1)/betabl(jrlv))
          endif
          if(k.ne.klist(i,j)) then
            dtemp=(swfrac(k)-swfrac(k+1))*sswflx(i,j)*delt1*g/
     .            (spcifh*max(onemm,dp(i,j,k,n)))
          else
            dtemp= swfrac(k)             *sswflx(i,j)*delt1*g/
     .            (spcifh*max(onemm,dp(i,j,k,n)))
          endif
          dsaln=0.
c
cdiag     if (i.eq.itest.and.j.eq.jtest) then
cdiag       write (lp,101) nstep,i+i0,j+j0,k,
cdiag.        1.-swfrac(k),1.-swfrac(k+1),dtemp
cdiag       call flush(lp)
cdiag     endif
        endif
c
c --- modify t and s; set old value arrays at p points for initial iteration
        if (k.le.klist(i,j)) then
          temp(i,j,k,n)=temp(i,j,k,n)+dtemp
          saln(i,j,k,n)=saln(i,j,k,n)+dsaln
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
          told (k)=temp(i,j,k,n)
          sold (k)=saln(i,j,k,n)
          thold(k)=th3d(i,j,k,n)
          uold (k)=.5*(u(i,j,k,n)+u(i+1,j  ,k,n))
          vold (k)=.5*(v(i,j,k,n)+v(i  ,j+1,k,n))
        endif
      enddo
c
      k=klist(i,j)
      ka=k+1
      kb=min(ka,kk)
      told (ka)=temp(i,j,kb,n)
      sold (ka)=saln(i,j,kb,n)
      thold(ka)=th3d(i,j,kb,n)
      uold (ka)=.5*(u(i,j,k,n)+u(i+1,j  ,k,n))
      vold (ka)=.5*(v(i,j,k,n)+v(i  ,j+1,k,n))
c
c --- calculate z at vertical grid levels - this array is the z values in m
c --- at the mid-depth of each micom layer except for index klist+1, where it
c --- is the z value of the bottom
c
c --- calculate layer thicknesses in m
      do k=1,kk
        if (k.eq.1) then
          hwide(k)=dp(i,j,k,n)/onem
          zgrid(i,j,k)=-.5*hwide(k)
        else if (k.lt.klist(i,j)) then
          hwide(k)=dp(i,j,k,n)/onem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
        else if (k.eq.klist(i,j)) then
          hwide(k)=dp(i,j,k,n)/onem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
          zgrid(i,j,k+1)=zgrid(i,j,k)-.5*hwide(k)
        else
          hwide(k)=0.
        endif
      enddo
c
c --- perform niter iterations to execute the semi-implicit solution of the
c --- diffusion equation. at least two iterations are recommended
c
      do iter=1,niter
c
c --- calculate layer variables required to estimate bulk richardson number
c
c --- calculate nearsurface reference variables. these variables are 
c --- equal to upper layer values if the upper layer is at least
c --- 7.5 m thick, and are otherwise averaged over 2*epsilon*zgrid.
        do k=1,kk
c
          if (k.le.klist(i,j)) then
            if (dp00.ge.7.5*onem) then
              zref=epsilon*zgrid(i,j,1)
              uref=uold(1)
              vref=vold(1)
              bref=-g*thref*(thold(1)+thbase)
            else
              zref=epsilon*zgrid(i,j,k)
              wz=max(zgrid(i,j,1),zref)
              uref=uold(1)*wz/zref
              vref=vold(1)*wz/zref
              bref=-g*thref*(thold(1)+thbase)*wz/zref
              do ka=1,kk-1
                if (zref.ge.zgrid(i,j,ka)) go to 10
                wz=min(zgrid(i,j,ka)-zgrid(i,j,ka+1),zgrid(i,j,ka)-zref) 
                del=0.5*wz/(zgrid(i,j,ka)-zgrid(i,j,ka+1))
                uref=uref-wz*(uold(ka)+del*(uold(ka+1)-uold(ka)))/zref
                vref=vref-wz*(vold(ka)+del*(vold(ka+1)-vold(ka)))/zref
                bref=bref+wz*g*thref*(thold(ka)+thbase+del*
     .              (thold(ka+1)-thold(ka)))/zref
              enddo
  10          continue
            endif
c
            ritop(k)=(zref-zgrid(i,j,k))*
     .               (bref+g*thref*(thold(k)+thbase))
            dvsq(k)=(uref-uold(k))**2+(vref-vold(k))**2
            if (-zgrid(i,j,k)*onem.lt.tenm) then
              swfrac(k)=redfac(jrlv)*
     .                  exp(zgrid(i,j,k)*onem/betard(jrlv))+
     .                  (1.-redfac(jrlv))*
     .                  exp(zgrid(i,j,k)*onem/betabl(jrlv))
            else
              swfrac(k)=(1.-redfac(jrlv))*
     .                  exp(zgrid(i,j,k)*onem/betabl(jrlv))
            endif
          endif
        enddo
c
c --- calculate interface variables required to estimate interior diffusivities
        do k=1,kk
          ka=k+1
          if (k.le.klist(i,j)) then
            shsq  (ka)=(uold(k)-uold(ka))**2+(vold(k)-vold(ka))**2
            alfadt(ka)=.5*(dsigdt(told(k ),sold(k ))+
     .                     dsigdt(told(ka),sold(ka)))*(told(k)-told(ka))
            betads(ka)=.5*(dsigds(told(k ),sold(k ))+
     .                     dsigds(told(ka),sold(ka)))*(sold(k)-sold(ka))
            dbloc (ka)=-g*thref*(thold(k)-thold(ka))
          endif
        enddo
c
c --- determine interior diffusivity profiles throughout the water column
c
        do k=1,kk+1
          vcty (i,j,k)  =0.0
          blmc(     k,1)=0.0
          blmc(     k,2)=0.0
          blmc(     k,3)=0.0
          ghats(i,j,k)  =0.0
          dift (i,j,k)  =0.0
          difs (i,j,k)  =0.0
        enddo
c
c --- shear instability plus background internal wave contributions
        do k=2,kk+1
          if (k-1.le.klist(i,j)) then
            if (shinst) then
              rigr=max(0.,dbloc(k)*(zgrid(i,j,k-1)-zgrid(i,j,k))/
     .            (shsq(k)+epsil))
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
            if (k-1.le.klist(i,j)) then
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
              else if ( alfadt(k).gt.0.0 .and. betads(k).lt.0.0
     .           .and. -alfadt(k).gt.betads(k)) then
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
cdiag   if (i.eq.itest.and.j.eq.jtest) then
cdiag      write (lp,102) (nstep,iter,i+i0,j+j0,k,
cdiag.     hwide(k),1.e4*vcty(i,j,k),1.e4*dift(i,j,k),1.e4*difs(i,j,k),
cdiag.       k=1,kk+1)
cdiag      call flush(lp)
cdiag   endif
c
c --- calculate boundary layer diffusivity profiles and match these to the
c --- previously-calculated interior diffusivity profiles
c
c --- diffusivities within the surface boundary layer are parameterized
c --- as a function of boundary layer thickness times a depth-dependent
c --- turbulent velocity scale (proportional to ustar) times a third-order
c --- polynomial shape function of depth. boundary layer diffusivities depend
c --- on surface forcing (the magnitude of this forcing and whether it is
c --- stabilizing or de-stabilizing) and the magnitude and gradient of interior
c --- mixing at the boundary layer base. boundary layer diffusivity profiles
c --- are smoothly matched to interior diffusivity profiles at the boundary
c --- layer base (the profiles and their first derivatives are continuous
c --- at z=-hbl). the turbulent boundary layer depth is diagnosed first, the
c --- boundary layer diffusivity profiles are calculated, then the boundary
c --- and interior diffusivity profiles are combined.
c 
c --- diagnose the new boundary layer depth as the depth where a bulk
c --- richardson number exceeds ric
c
c --- initialize hbl and nbl to bottomed out values
        k1=1
        k2=2
        rib(k1)=0.0
        nbl=klist(i,j)
        hbl=-zgrid(i,j,klist(i,j))
c
c --- diagnose hbl and nbl
        do k=2,kk
          if (k.le.nbl) then
            case=-zgrid(i,j,k)
            bfsfc=buoyfl(i,j)-swfrac(k)*buoysw(i,j)
            stable=.5+sign(.5,-bfsfc)
            dnorm=stable+(1.-stable)*epsilon
c
c --- compute turbulent velocity scales at dnorm, for
c --- hbl= case = -zgrid(i,j,k)
            call wscale(i,j,case,dnorm,bfsfc,wm,ws)
c
c --- compute the turbulent shear contribution to rib
            bvsq=.5*(dbloc(k  )/(zgrid(i,j,k-1)-zgrid(i,j,k  ))+
     .               dbloc(k+1)/(zgrid(i,j,k  )-zgrid(i,j,k+1)))
            vtsq=-zgrid(i,j,k)*ws*sqrt(abs(bvsq))*vtc
c
c --- compute bulk richardson number at new level
c --- linearly interpolate to find hbl as the depth where rib = ricr
c --- minimum hbl is 1 cm below -zgrid(i,j,1)+.5*hwide(1)
            rib(k2)=ritop(k)/(dvsq(k)+vtsq+epsil)
            if (nbl.eq.klist(i,j).and.rib(k2).ge.ricr) then
              hbl=-zgrid(i,j,k-1)+(zgrid(i,j,k-1)-zgrid(i,j,k))*
     .               (ricr-rib(k1))/(rib(k2)-rib(k1)+epsil)
              nbl=k
              hblmin=-zgrid(i,j,1)+.5*hwide(1)+0.01
              hblmax=-zgrid(i,j,klist(i,j))-0.01
              if (hbl.lt.hblmin) then
                hbl=hblmin
                nbl=2
              endif
              if (hbl.gt.hblmax) then
                hbl=hblmax
                nbl=klist(i,j)
              endif
            endif
          endif
c
          ksave=k1
          k1=k2
          k2=ksave
        enddo
c
c --- calculate swfrml, the fraction of solar radiation absorbed by depth hbl
        k=nbl
        q=(zgrid(i,j,k-1)+hbl)/(zgrid(i,j,k-1)-zgrid(i,j,k))
        swfrml=swfrac(k-1)+q*(swfrac(k)-swfrac(k-1))
c
c --- limit check on hbl for negative (stablizing) surface buoyancy forcing
        bfsfc=buoyfl(i,j)-swfrml*buoysw(i,j)
        stable=.5+sign(.5,-bfsfc)
        bfsfc=bfsfc-stable*epsil                 !insures bfsfc never=0
        if (bfsfc.lt.0.) then
          hmonob=-cmonob*ustar(i,j)*ustar(i,j)*ustar(i,j)/(vonk*bfsfc)
          hbl=max(min(hbl,hekman(i,j),hmonob,(p(i,j,kk+1)-onem)/onem),
     .           -zgrid(i,j,1)+.5*hwide(1)+0.01)
        endif
        nbl=klist(i,j)
c
c--- find new nbl and re-calculate swfrml
        do k=2,kk
          if (k.le.klist(i,j)) then
            if ((nbl.eq.klist(i,j)).and.(-zgrid(i,j,k).gt.hbl)) then
              nbl=k
              q=(zgrid(i,j,k-1)+hbl)/(zgrid(i,j,k-1)-zgrid(i,j,k))
              swfrml=swfrac(k-1)+q*(swfrac(k)-swfrac(k-1))
            endif
          endif
        enddo
c
c --- find forcing stability and buoyancy forcing for final hbl values
c --- determine case (for case=0., hbl lies between -zgrid(i,j,nbl)
c --- and the interface above. for case=1., hbl lies between 
c --- -zgrid(i,j,nbl-1) and the interface below)
c
c --- velocity scales at hbl
        bfsfc=buoyfl(i,j)-swfrml*buoysw(i,j)
        stable=.5+sign(.5,-bfsfc)
        bfsfc=bfsfc-stable*epsil                  !insures bfsfc never=0
        dnorm=stable+(1.-stable)*epsilon
        case=.5+sign(.5,-zgrid(i,j,nbl)-.5*hwide(nbl)-hbl)
c
        call wscale(i,j,hbl,dnorm,bfsfc,wm,ws)
c
c --- compute the boundary layer diffusivity profiles. first, find interior
c --- viscosities and their vertical derivatives at hbl
        ka=ifix(case+epsil)*(nbl-1)+(1-ifix(case+epsil))*nbl
        q=(hbl*onem-p(i,j,ka))/dp(i,j,ka,n)
        vctyh=vcty(i,j,ka)+q*(vcty(i,j,ka+1)-vcty(i,j,ka))
        difsh=difs(i,j,ka)+q*(difs(i,j,ka+1)-difs(i,j,ka))
        difth=dift(i,j,ka)+q*(dift(i,j,ka+1)-dift(i,j,ka))
c
        q=(hbl+zgrid(i,j,nbl-1))/(zgrid(i,j,nbl-1)-zgrid(i,j,nbl))
        dvdzup=(vcty(i,j,nbl-1)-vcty(i,j,nbl  ))/hwide(nbl-1)
        dvdzdn=(vcty(i,j,nbl  )-vcty(i,j,nbl+1))/hwide(nbl  )
        viscp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
        dvdzup=(difs(i,j,nbl-1)-difs(i,j,nbl  ))/hwide(nbl-1)
        dvdzdn=(difs(i,j,nbl  )-difs(i,j,nbl+1))/hwide(nbl  )
        difsp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
        dvdzup=(dift(i,j,nbl-1)-dift(i,j,nbl  ))/hwide(nbl-1) 
        dvdzdn=(dift(i,j,nbl  )-dift(i,j,nbl+1))/hwide(nbl  )
        diftp=.5*((1.-q)*(dvdzup+abs(dvdzup))+q*(dvdzdn+abs(dvdzdn)))
c
        f1=-stable*c11*bfsfc/(ustar(i,j)**4+epsil) 
c
        gat1(1)=vctyh/hbl/(wm+epsil)
        dat1(1)=min(0.,-viscp/(wm+epsil)+f1*vctyh)
c
        gat1(2)=difsh/hbl/(ws+epsil)
        dat1(2)=min(0.,-difsp/(ws+epsil)+f1*difsh) 
c
        gat1(3)=difth/hbl/(ws+epsil)
        dat1(3)=min(0.,-diftp/(ws+epsil)+f1*difth)
c
c --- compute turbulent velocity scales on the interfaces
        do k=2,kk+1
          if (k.le.min(nbl+1,klist(i,j))) then
            sigg=p(i,j,k)/(hbl*onem)
            dnorm=stable*sigg+(1.-stable)*min(sigg,epsilon)
c
            call wscale(i,j,hbl,dnorm,bfsfc,wm,ws)
c
c --- compute the dimensionless shape functions at the interfaces
            aa1=sigg-2.
            aa2=3.-2.*sigg
            aa3=sigg-1.
c
            gm=aa1+aa2*gat1(1)+aa3*dat1(1) 
            gs=aa1+aa2*gat1(2)+aa3*dat1(2)
            gt=aa1+aa2*gat1(3)+aa3*dat1(3)
c
c --- compute boundary layer diffusivities at the interfaces
            blmc(k,1)=hbl*wm*sigg*(1.+sigg*gm)
            blmc(k,2)=hbl*ws*sigg*(1.+sigg*gs)
            blmc(k,3)=hbl*ws*sigg*(1.+sigg*gt)
c
c --- compute nonlocal transport forcing term = ghats * <ws>o
            if (nonloc) ghats(i,j,k)=(1.-stable)*cg/(ws*hbl+epsil)
          endif
        enddo
c
c
c --- enhance diffusivities on the interface closest to hbl
c
c --- first compute diffusivities at nbl-1 grid level 
        sigg=-zgrid(i,j,nbl-1)/hbl
        dnorm=stable*sigg+(1.-stable)*min(sigg,epsilon)
c
        call wscale(i,j,hbl,dnorm,bfsfc,wm,ws)
c
        sigg=-zgrid(i,j,nbl-1)/hbl
        aa1=sigg-2.
        aa2=3.-2.*sigg
        aa3=sigg-1.
        gm=aa1+aa2*gat1(1)+aa3*dat1(1)
        gs=aa1+aa2*gat1(2)+aa3*dat1(2)
        gt=aa1+aa2*gat1(3)+aa3*dat1(3)
        dkm1(1)=hbl*wm*sigg*(1.+sigg*gm)
        dkm1(2)=hbl*ws*sigg*(1.+sigg*gs)
        dkm1(3)=hbl*ws*sigg*(1 +sigg*gt)
c
c --- now enhance diffusivity at interface nbl
c
c --- this procedure was altered for hycom to reduce diffusivity enhancement
c --- if the interface in question is located more than dp0enh below hbl.
c --- this prevents enhanced boundary layer mixing from penetrating too far
c --- below hbl when the m.l. base is located in a very thick layer
        k=nbl-1
        ka=k+1
        delta=(hbl+zgrid(i,j,k))/(zgrid(i,j,k)-zgrid(i,j,ka))
        q=1.-case*max(0.,min(1.,(p(i,j,ka)-hbl*onem-dp0enh)/dp0enh))
c
        dkmp2=case*vcty(i,j,ka)+(1.-case)*blmc(ka,1)
        dstar=(1.-delta)**2*dkm1(1)+delta**2*dkmp2      
        blmc(ka,1)=(1.-delta)*vcty(i,j,ka)+delta*dstar
c
        dkmp2=case*difs(i,j,ka)+(1.-case)*blmc(ka,2)
        dstar=(1.-delta)**2*dkm1(2)+delta**2*dkmp2    
        blmc(ka,2)=(1.-delta)*difs(i,j,ka)+delta*dstar
c
        dkmp2=case*dift(i,j,ka)+(1.-case)*blmc(ka,3)
        dstar=(1.-delta)**2*dkm1(3)+delta**2*dkmp2     
        blmc(ka,3)=(1.-delta)*dift(i,j,ka)+delta*dstar
c
        if (case.eq.1.) then
          blmc(ka,1)=max(vcty(i,j,ka),q*blmc(ka,1))
          blmc(ka,2)=max(difs(i,j,ka),q*blmc(ka,2))
          blmc(ka,3)=max(dift(i,j,ka),q*blmc(ka,3))
        endif
c
        if (nonloc) ghats(i,j,ka)=(1.-case)*ghats(i,j,ka)
c
c --- combine interior and boundary layer coefficients and nonlocal term
c --- zero scalar bottom diffusivities to assure no bottom flux
        do k=2,nbl
          vcty(i,j,k)=blmc(k,1)
          difs(i,j,k)=blmc(k,2)
          dift(i,j,k)=blmc(k,3)
        enddo
        do k=nbl+1,klist(i,j)
          ghats(i,j,k)=0.0
        enddo
        do k=klist(i,j)+1,kk+1
          vcty(i,j,k)=difmiw
          difs(i,j,k)=difsiw
          dift(i,j,k)=difsiw
          ghats(i,j,k)=0.0
        enddo
c
cdiag   if (i.eq.itest.and.j.eq.jtest) then
cdiag     write (lp,103) (nstep,iter,i+i0,j+j0,k,
cdiag.    hwide(k),1.e4*vcty(i,j,k),1.e4*dift(i,j,k),1.e4*difs(i,j,k),
cdiag.      ghats(i,j,k),k=1,kk+1)
cdiag     call flush(lp)
cdiag   endif
c
c --- save array dpbl=onem*hbl for output and diagnosis
        dpbl(i,j)=onem*hbl
c
        if (iter.lt.niter) then
c
c --- perform the vertical mixing at p points
c
          do k=1,klist(i,j)
            difft(k+1)=dift(i,j,k+1)
            diffs(k+1)=difs(i,j,k+1)
            diffm(k+1)=vcty(i,j,k+1)
            ghat(k+1)=ghats(i,j,k+1)
            t1do(k)=temp(i,j,k,n)
            s1do(k)=saln(i,j,k,n)
            if (trcout) tr1do(k)=tracer(i,j,k)
            u1do(k)=.5*(u(i,j,k,n)+u(i+1,j  ,k,n))
            v1do(k)=.5*(v(i,j,k,n)+v(i  ,j+1,k,n))
            hm(k)=hwide(k)
            zm(k)=zgrid(i,j,k)
          enddo
c
          nlayer=klist(i,j)
          k=nlayer+1
          ka=min(k,kk)
          difft(k)=0.0
          diffs(k)=0.0
          diffm(k)=0.0
          ghat(k)=0.0
          t1do(k)=temp(i,j,ka,n)
          s1do(k)=saln(i,j,ka,n)
          if (trcout) tr1do(k)=tracer(i,j,ka)
          u1do(k)=u1do(k-1)
          v1do(k)=v1do(k-1)
          zm(k)=zgrid(i,j,k)
c
c --- compute factors for coefficients of tridiagonal matrix elements.
c         tri(k=1:NZ,0) : dt/hwide(k)/ dzb(k-1)=z(k-1)-z(k)=dzabove)
c         tri(k=1:NZ,1) : dt/hwide(k)/(dzb(k  )=z(k)-z(k+1)=dzbelow)
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
          ghatflux=-(surflx(i,j)-sswflx(i,j))*thref/spcifh
          call tridcof(difft,tri,nlayer,tcu,tcc,tcl)
          call tridrhs(hm,t1do,difft,ghat,ghatflux,tri,nlayer,rhs)
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,t1do,t1dn,difft)
c
c --- s solution
          ghatflux=-salflx(i,j)*thref
          call tridcof(diffs,tri,nlayer,tcu,tcc,tcl)
          call tridrhs(hm,s1do,diffs,ghat,ghatflux,tri,nlayer,rhs)
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,s1do,s1dn,diffs)
c
c --- tracer solution
          if (trcout) then
            ghatflux=0.
            call tridrhs(hm,tr1do,diffs,ghat,ghatflux,tri,nlayer,rhs)
            call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,tr1do,tr1dn,diffs)
          endif
c
cdiag     if (i.eq.itest.and.j.eq.jtest) then
cdiag       write (lp,104) (nstep,iter,i+i0,j+j0,k,
cdiag.        hm(k),t1do(k),t1dn(k),s1do(k),s1dn(k),tr1do(k),tr1dn(k),
cdiag.        k=1,nlayer)
cdiag       call flush(lp)
cdiag     endif
c
c --- u solution
          call tridcof(diffm,tri,nlayer,tcu,tcc,tcl)
          do k=1,nlayer
            rhs(k)=u1do(k)
          enddo
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,u1do,u1dn,diffm)
c
c --- v solution
          do k=1,nlayer
            rhs(k)=v1do(k)
          enddo
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,v1do,v1dn,diffm)
c
cdiag     if (i.eq.itest.and.j.eq.jtest) then
cdiag       write (lp,105) (nstep,iter,i+i0,j+j0,k,
cdiag.        hm(k),u1do(k),u1dn(k),v1do(k),v1dn(k),k=1,nlayer)
cdiag       call flush(lp)
cdiag     endif
c
c --- reset old variables in preparation for next iteration
          do k=1,nlayer+1
            told(k)=t1dn(k)
            sold(k)=s1dn(k)
            thold(k)=sig(told(k),sold(k))-thbase
            if (iter.lt.niter) then
              uold(k)=u1dn(k)
              vold(k)=v1dn(k)
            endif
          enddo
        endif                         ! iter < niter
c
      enddo                           ! iteration loop
c
 101  format(i9,3i4,'absorbup,dn,dtemp,dsaln ',2f6.3,2f10.6)
 102  format(25x,'   thick      viscty    t diff    s diff  '
     .     /(i9,i2,3i4,2x,4f10.2))
 103  format(25x,'   thick      viscty    t diff    s diff   nonlocal'
     .     /(i9,i2,3i4,2x,4f10.2,f11.6))
 104  format(25x,
     .     '  thick   t old   t new   s old   s new trc old trc new'
     .     /(i9,i2,3i4,1x,f9.2,4f8.3,2f7.4))
 105  format(25x,'   thick   u old   u new   v old   v new'
     .     /(i9,i2,3i4,1x,f10.2,4f8.3))
c
      return
      end
c
      subroutine mxkppbij(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j
c
c -------------------------------------------------
c --- kpp vertical diffusion, single j-row (part B)
c -------------------------------------------------
c
c --- perform the final vertical mixing at p points
c
c --- local 1-d arrays for matrix solution
      real t1do(kdm+1),t1dn(kdm+1),s1do(kdm+1),s1dn(kdm+1),
     .     tr1do(kdm+1),tr1dn(kdm+1),
     .     difft(kdm+1),diffs(kdm+1),
     .     ghat(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     +     tcc(kdm),         ! central ...     (k  ) ..
     +     tcl(kdm),         ! lower .....     (k-1) ..
     +     rhs(kdm)          ! right-hand-side terms
c
      real    ghatflux
      integer k,ka,nlayer
c
      include 'stmt_fns.h'
c
      do k=1,klist(i,j)
        difft(k+1)=dift(i,j,k+1)
        diffs(k+1)=difs(i,j,k+1)
        ghat(k+1)=ghats(i,j,k+1)
        t1do(k)=temp(i,j,k,n)
        s1do(k)=saln(i,j,k,n)
        if (trcout) tr1do(k)=tracer(i,j,k)
        hm(k)=dp(i,j,k,n)/onem
        zm(k)=zgrid(i,j,k)
      enddo
c
      nlayer=klist(i,j)
      k=nlayer+1
      ka=min(k,kk)
      difft(k)=0.0
      diffs(k)=0.0
      ghat(k)=0.0
      t1do(k)=temp(i,j,ka,n)
      s1do(k)=saln(i,j,ka,n)
      if (trcout) tr1do(k)=tracer(i,j,ka)
      zm(k)=zgrid(i,j,k)
c
c --- compute factors for coefficients of tridiagonal matrix elements.
c     tri(k=1:NZ,0) : dt/hwide(k)/ dzb(k-1)=z(k-1)-z(k)=dzabove)
c     tri(k=1:NZ,1) : dt/hwide(k)/(dzb(k  )=z(k)-z(k+1)=dzbelow)
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
      ghatflux=-(surflx(i,j)-sswflx(i,j))*thref/spcifh
      call tridcof(difft,tri,nlayer,tcu,tcc,tcl)
      call tridrhs(hm,t1do,difft,ghat,ghatflux,tri,nlayer,rhs)
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,t1do,t1dn,difft)
c
c --- s solution
      ghatflux=-salflx(i,j)*thref
      call tridcof(diffs,tri,nlayer,tcu,tcc,tcl)
      call tridrhs(hm,s1do,diffs,ghat,ghatflux,tri,nlayer,rhs)
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,s1do,s1dn,diffs)
c
c --- tracer solution
      if (trcout) then
        ghatflux=0.
        call tridrhs(hm,tr1do,diffs,ghat,ghatflux,tri,nlayer,rhs)
        call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,tr1do,tr1dn,diffs)
      endif
c
c --- adjust t, s, th, arrays
      do k=1,kk
        if (k.le.klist(i,j)) then
          temp(i,j,k,n)=t1dn(k)
          saln(i,j,k,n)=s1dn(k)
          th3d(i,j,k,n)=sig(t1dn(k),s1dn(k))-thbase
          if (trcout) tracer(i,j,k)=tr1dn(k)
        endif
      enddo
c
      return
      end
c
      subroutine mxkppciju(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j
c
c ----------------------------------------------------------------------------
c --- kpp vertical diffusion, single j-row (part A), momentum at u grid points
c ----------------------------------------------------------------------------
c
c local variables for kpp mixing
c
c --- local 1-d arrays for matrix solution
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
      real presu
      integer k,ka,nlayer,kmask(idm)
c
      presu=0.
      kmask(1)=0
      do k=1,kk+1
        ka=min(k,kk)
        if (dpu(i,j,ka,n).le.tencm.or.k.eq.kk+1) kmask(1)=1
        if (presu.lt.depthu(i,j)-tencm.and.kmask(1).eq.0) then
          diffm(k+1)=.5*(vcty(i,j,k+1)+vcty(i-1,j,k+1))
          u1do(k)=u(i,j,k,n)
          hm(k)=dpu(i,j,k,n)/onem
          if (k.eq.1) then
            zm(k)=-.5*hm(k)
          else
            zm(k)=zm(k-1)-.5*(hm(k-1)+hm(k))
          endif
          presu=presu+dpu(i,j,k,n)
          nlayer=k
        else if (k.eq.nlayer+1) then
          diffm(k)=0.
          u1do(k)=u1do(k-1)
          zm(k)=zm(k-1)-.5*hm(k-1)
        endif
      enddo
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
cdiag if (i.eq.itest.and.j.eq.jtest) then
cdiag   write (lp,106) (nstep,i+i0,j+j0,k,
cdiag.    hm(k),u1do(k),u1dn(k),k=1,nlayer)
cdiag   call flush(lp)
cdiag endif
      return
 106  format(23x,'   thick   u old   u new'/(i9,3i4,1x,f10.3,2f8.3))
      end
      subroutine mxkppcijv(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j
c
c ----------------------------------------------------------------------------
c --- kpp vertical diffusion, single j-row (part A), momentum at v grid points
c ----------------------------------------------------------------------------
c
c local variables for kpp mixing
c
c --- local 1-d arrays for matrix solution
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
      real presv
      integer k,ka,nlayer,kmask(idm)
c
      presv=0.
      kmask(1)=0
      do k=1,kk+1
        ka=min(k,kk)
        if (dpv(i,j,ka,n).le.tencm.or.k.eq.kk+1) kmask(1)=1
        if (presv.lt.depthv(i,j)-tencm.and.kmask(1).eq.0) then
          diffm(k+1)=.5*(vcty(i,j,k+1)+vcty(i,j-1,k+1))
          v1do(k)=v(i,j,k,n)
          hm(k)=dpv(i,j,k,n)/onem
          if (k.eq.1) then
            zm(k)=-.5*hm(k)
          else
            zm(k)=zm(k-1)-.5*(hm(k-1)+hm(k))
          endif
          presv=presv+dpv(i,j,k,n)
          nlayer=k
        else if (k.eq.nlayer+1) then
          diffm(k)=0.
          v1do(k)=v1do(k-1)
          zm(k)=zm(k-1)-.5*hm(k-1)
        endif
      enddo
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
cdiag if (i.eq.itest.and.j.eq.jtest) then
cdiag   write (lp,107) (nstep,i+i0,j+j0,k,
cdiag.    hm(k),v1do(k),v1dn(k),k=1,nlayer)
cdiag   call flush(lp)
cdiag endif
      return
 107  format(23x,'   thick   v old   v new'/(i9,3i4,1x,f10.3,2f8.3))
      end
      subroutine wscale(i,j,zlevel,dnorm,bfsfc,wm,ws)
      use mod_xc  ! HYCOM communication interface
c
c -------------------------------------------------------------------------
c --- subroutine to compute turbulent velocity scales for kpp mixing scheme
c -------------------------------------------------------------------------
c
      implicit none
c
      include 'common_blocks.h'
c
c --- see inikpp for initialization of /kppltr/ and other constants.
c
      integer    nzehat,nustar
      parameter (nzehat=890,nustar=192)
c
      real, dimension (0:nzehat+1,0:nustar+1) ::
     & wmt            ! momentum velocity scale table
     &,wst            ! scalar   velocity scale table
      common/kppltr/ wmt,wst
      save  /kppltr/
c
      integer i,j
      real    zlevel,dnorm,bfsfc,wm,ws
c
      real    zdiff,udiff,zfrac,ufrac,
     .        wam,wbm,was,wbs,ucube,zehat
      integer iz,izp1,ju,jup1
c
c --- use lookup table for zehat < zmax  only;  otherwise use stable formulae
c
      zehat=-vonk*dnorm*zlevel*bfsfc
      if (zehat.le.zmax) then
        zdiff=zehat-zmin
        iz=int(zdiff/deltaz)
        iz=max(min(iz,nzehat),0)
        izp1=iz+1
c
        udiff=ustar(i,j)-umin
        ju=int(udiff/deltau)
        ju=max(min(ju,nustar),0)
        jup1=ju+1
c
        zfrac=zdiff/deltaz-iz
        ufrac=udiff/deltau-ju
c
        wam=(1.-zfrac)*wmt(iz,jup1)+zfrac*wmt(izp1,jup1)
        wbm=(1.-zfrac)*wmt(iz,ju  )+zfrac*wmt(izp1,ju  )
        wm =(1.-ufrac)*wbm         +ufrac*wam
c
        was=(1.-zfrac)*wst(iz,jup1)+zfrac*wst(izp1,jup1)
        wbs=(1.-zfrac)*wst(iz,ju  )+zfrac*wst(izp1,ju  )
        ws =(1.-ufrac)*wbs         +ufrac*was
c
      else
c
        ucube=ustar(i,j)**3
        wm=vonk*ustar(i,j)*ucube/(ucube+c11*zehat)
        ws=wm
c
      endif
c
      return
      end
c
c
c> Revision history:
c>
c> Jun  2000 - conversion to SI units.
c> Jul  2000 - included wscale in this file to facilitate in-lining
