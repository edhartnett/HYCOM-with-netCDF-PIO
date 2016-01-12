      subroutine mxkprf(m,n)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
c
c --- hycom version 2.1
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c ---------------------------------------------------------
c --- k-profile vertical mixing models
c ---   a) large, mc williams, doney kpp vertical diffusion
c ---   b) mellor-yamada 2.5 vertical diffusion
c ---------------------------------------------------------
c
      logical    lpipe_mxkprf
      parameter (lpipe_mxkprf=.false.)
c
      real      delp,sigmlj
      integer   i,j,k,l
      character text*12
c
      include 'stmt_fns.h'
c
      if (mxlmy) then
        call xctilr(u(      1-nbdy,1-nbdy,1,m),1,kk, 1,1, halo_uv)
        call xctilr(u(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
        call xctilr(v(      1-nbdy,1-nbdy,1,m),1,kk, 1,1, halo_vv)
        call xctilr(v(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_vv)
        call xctilr(p(      1-nbdy,1-nbdy,2  ),1,kk, 1,1, halo_ps)
        call xctilr(ubavg(  1-nbdy,1-nbdy,  m),1, 1, 1,1, halo_uv)
        call xctilr(vbavg(  1-nbdy,1-nbdy,  m),1, 1, 1,1, halo_vv)
      else
        call xctilr(u(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
        call xctilr(v(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_vv)
        call xctilr(p(      1-nbdy,1-nbdy,2  ),1,kk, 1,1, halo_ps)
      endif
c
      margin = 0  ! no horizontal derivatives
c
c --- diffisuvity/viscosity calculation
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call mxkprfaj(m,n, j)
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
      if     (lpipe .and. lpipe_mxkprf) then
c ---   compare two model runs.
        util6 = klist
        write (text,'(a12)') 'klist       '
        call pipe_compare_sym1(util6,ip,text)
        if (mxlmy) then
          do k= 0,kk+1
            write (text,'(a9,i3)') 'q2     k=',k
            call pipe_compare_sym1(q2( 1-nbdy,1-nbdy,k,n),ip,text)
            write (text,'(a9,i3)') 'q2l    k=',k
            call pipe_compare_sym1(q2l(1-nbdy,1-nbdy,k,n),ip,text)
            write (text,'(a9,i3)') 'difqmy k=',k
            call pipe_compare_sym1(difqmy(1-nbdy,1-nbdy,k),ip,text)
            write (text,'(a9,i3)') 'diftmy k=',k
            call pipe_compare_sym1(diftmy(1-nbdy,1-nbdy,k),ip,text)
            write (text,'(a9,i3)') 'vctymy k=',k
            call pipe_compare_sym1(vctymy(1-nbdy,1-nbdy,k),ip,text)
          enddo
        endif
        do k= 1,kk
          write (text,'(a9,i3)') 'vcty   k=',k
          call pipe_compare_sym1(vcty(1-nbdy,1-nbdy,k),ip,text)
          write (text,'(a9,i3)') 'dift   k=',k
          call pipe_compare_sym1(dift(1-nbdy,1-nbdy,k),ip,text)
          write (text,'(a9,i3)') 'difs   k=',k
          call pipe_compare_sym1(difs(1-nbdy,1-nbdy,k),ip,text)
        enddo
      endif
c
c ---   final mixing of variables at p points
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call mxkprfbj(m,n, j)
      enddo
!$OMP END PARALLEL DO
c
c --- final velocity mixing at u,v points
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call mxkprfcj(m,n, j)
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
     &                              onem*(-zgrid(i,j,k-1)+
     &               ((zgrid(i,j,k-1)-zgrid(i,j,k))*
     &               (th3d(i,j,1,n)+sigmlj-th3d(i,j,k-1,n)))/
     &               (th3d(i,j,k,n) +epsil -th3d(i,j,k-1,n))))
                  exit
                endif
              enddo
            enddo
          enddo
        enddo
c
!$OMP   END PARALLEL DO
c
        call xctilr(p(     1-nbdy,1-nbdy,2),1,kk, 1,1, halo_ps)
        call xctilr(dpmixl(1-nbdy,1-nbdy,n),1, 1, 1,1, halo_ps)
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
     &              -min(p(i,j,k  ),dpmixl(i,j,n))
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
     &             (min(p(i,j,k+1)+p(i-1,j,k+1),
     &                  dpmixl(i,j,n)+dpmixl(i-1,j,n))
     &             -min(p(i,j,k  )+p(i-1,j,k  ),
     &                  dpmixl(i,j,n)+dpmixl(i-1,j,n)))
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
     &             (min(p(i,j,k+1)+p(i,j-1,k+1),
     &                  dpmixl(i,j,n)+dpmixl(i,j-1,n))
     &             -min(p(i,j,k  )+p(i,j-1,k  ),
     &                  dpmixl(i,j,n)+dpmixl(i,j-1,n)))
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
      subroutine mxkprfaj(m,n, j)
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
      if (mxlkpp) then
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            call mxkppaij(m,n, i,j)
          enddo
        enddo
      else if (mxlmy) then
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            call mxmyaij(m,n, i,j)
          enddo
        enddo
      endif
c
      return
      end
      subroutine mxkprfbj(m,n, j)
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
          call mxkprfbij(m,n, i,j)
        enddo
      enddo
c
      return
      end
      subroutine mxkprfcj(m,n, j)
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
          call mxkprfciju(m,n, i,j)
        enddo
      enddo
c
      do l=1,isv(j)
        do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
          call mxkprfcijv(m,n, i,j)
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
c -------------------------------------------------------------
c --- kpp vertical diffusion, single j-row (part A)
c --- vertical coordinate is z negative below the ocean surface
c -------------------------------------------------------------
c
c local variables for kpp mixing
      real delta               ! fraction hbl lies beteen zgrid neighbors
      real zrefmn              ! nearsurface reference z, minimum
      real zref                ! nearsurface reference z
      real wref,qwref          ! nearsurface reference width,inverse
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
      real hwide(kdm)          ! layer thicknesses in m (minimum 1mm)
      real dpmm(kdm)           !     max(onemm,dp(i,j,:,n))
      real qdpmm(kdm)          ! 1.0/max(onemm,dp(i,j,:,n))
      real pij(kdm+1)          ! local copy of p(i,j,:)
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
      real tmn                 ! time averaged SST
      real smn                 ! time averaged SSS
      real dsgdt               ! dsigdt(tmn,smn)
      real buoyfl              ! total     surface buoyancy (into atmos.)
      real buoysw              ! shortwave surface buoyancy (into atmos.)
      real bfsfc               ! surface buoyancy forcing   (into atmos.)
c
      integer nbl              ! layer containing boundary layer base
      integer k1,k2            ! bulk richardson number indices
c
c --- local 1-d arrays for matrix solution
      real u1do(kdm+1),u1dn(kdm+1),v1do(kdm+1),v1dn(kdm+1),t1do(kdm+1),
     &     t1dn(kdm+1),s1do(kdm+1),s1dn(kdm+1),
     &     diffm(kdm+1),difft(kdm+1),diffs(kdm+1),
     &     ghat(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
c
c --- local 1-d arrays for iteration loops
      real uold(kdm+1),vold (kdm+1),told (kdm+1),
     &     sold(kdm+1),thold(kdm+1)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     &     tcc(kdm),         ! central ...     (k  ) ..
     &     tcl(kdm),         ! lower .....     (k-1) ..
     &     rhs(kdm)          ! right-hand-side terms
c
      real dtemp,dsaln,wq,wt,ratio,q,ghatflux,
     &     dvdzup,dvdzdn,viscp,difsp,diftp,f1,sigg,aa1,aa2,aa3,gm,gs,gt,
     &     dkmp2,dstar,hmonob,hblmin,hblmax,sflux1,bvsq,vtsq,
     &     vctyh,difsh,difth,zrefo,qspcifh
c
      integer k,ka,kb,nlayer,ksave,iter,jrlv
c
      include 'stmt_fns.h'
c
c --- locate lowest substantial mass-containing layer.
      pij(1)=p(i,j,1)
      do k=1,kk
        dpmm( k)  =max(onemm,dp(i,j,k,n))
        qdpmm(k)  =1.0/dpmm(k)
        pij(  k+1)=pij(k)+dp(i,j,k,n)
        p(i,j,k+1)=pij(k+1)
      enddo
      do k=kk,1,-1
        if (dpmm(k).gt.tencm) then
          exit
        endif
      enddo
      klist(i,j)=max(k,2)  !always consider at least 2 layers
c
c --- forcing of t,s by surface fluxes. shortwave flux penetration depends on
c --- jerlov water type, which is set in thermf.f.  flux positive into ocean
c
      jrlv=jerlov(i,j)
      qspcifh=1.0/spcifh
c
      do k=1,kk
        if (thermo) then
          if (k.eq.1) then
            if (pij(k+1).lt.tenm) then
              swfrac(k+1)=redfac(jrlv) *exp(-pij(k+1)*betard(jrlv))+
     &                (1.-redfac(jrlv))*exp(-pij(k+1)*betabl(jrlv))
            else
              swfrac(k+1)=(1.-redfac(jrlv))*
     &                    exp(-pij(k+1)*betabl(jrlv))
            endif
            sflux1=surflx(i,j)-sswflx(i,j)
            dtemp=(sflux1+(1.-swfrac(k+1))*sswflx(i,j))*
     &            delt1*g*qspcifh*qdpmm(k)
            dsaln=salflx(i,j)*
     &            delt1*g*        qdpmm(k) 
cdiag       if (i.eq.itest.and.j.eq.jtest) then
cdiag         write (lp,101) nstep,i+i0,j+j0,k, 
cdiag&          0.,1.-swfrac(k+1),dtemp,dsaln
cdiag         call flush(lp)
cdiag       endif
          elseif (k.le.klist(i,j)) then
            if (pij(k+1).lt.tenm) then
              swfrac(k+1)=redfac(jrlv) *exp(-pij(k+1)*betard(jrlv))+
     &                (1.-redfac(jrlv))*exp(-pij(k+1)*betabl(jrlv))
            else
              swfrac(k+1)=(1.-redfac(jrlv))*
     &                    exp(-pij(k+1)*betabl(jrlv))
            endif
            if(k.ne.klist(i,j)) then
              dtemp=(swfrac(k)-swfrac(k+1))*sswflx(i,j)*
     &              delt1*g*qspcifh*qdpmm(k)
            else
              dtemp= swfrac(k)             *sswflx(i,j)*
     &               delt1*g*qspcifh*qdpmm(k)
            endif
            dsaln=0.0
cdiag       if (i.eq.itest.and.j.eq.jtest) then
cdiag         write (lp,101) nstep,i+i0,j+j0,k,
cdiag&          1.-swfrac(k),1.-swfrac(k+1),dtemp
cdiag         call flush(lp)
cdiag       endif
          else !k.gt.klist(i,j)
            dtemp=0.0
            dsaln=0.0
          endif 
        else !.not.thermo
          dtemp=0.0
          dsaln=0.0
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
          hwide(k)=dpmm(k)*qonem
          zgrid(i,j,k)=-.5*hwide(k)
        else if (k.lt.klist(i,j)) then
          hwide(k)=dpmm(k)*qonem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
        else if (k.eq.klist(i,j)) then
          hwide(k)=dpmm(k)*qonem
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
c --- calculate nearsurface reference variables,
c --- averaged over -2*epsilon*zgrid, but no more than 8m.
        zrefmn = -4.0
        zrefo  =  1.0  ! impossible value
        do k=1,klist(i,j)
          zref=max(epsilon*zgrid(i,j,k),zrefmn)  ! nearest to zero
          if     (zref.ne.zrefo) then  ! new zref
            wref =-2.0*zref
            qwref=1.0/wref
            wq=min(hwide(1),wref)*qwref
            uref=uold(1)*wq
            vref=vold(1)*wq
            bref=-g*thref*(thold(1)+thbase)*wq
            wt=0.0
            do ka=2,k
              wt=wt+wq
              if (wt.ge.1.0) then
                exit
              endif
              wq=min(1.0-wt,hwide(ka)*qwref)
              uref=uref+uold(ka)*wq
              vref=vref+vold(ka)*wq
              bref=bref-g*thref*(thold(ka)+thbase)*wq
            enddo
          endif
          zrefo=zref
c
          ritop(k)=(zref-zgrid(i,j,k))*
     &             (bref+g*thref*(thold(k)+thbase))
          dvsq(k)=(uref-uold(k))**2+(vref-vold(k))**2
*
*         if (i.eq.itest.and.j.eq.jtest) then
*           if     (k.eq.1) then
*             write(lp,'(3a)')
*    &          ' k        z  zref',
*    &          '      u   uref      v   vref',
*    &          '      b   bref    ritop   dvsq'
*           endif
*           write(lp,'(i2,f9.2,f6.2,4f7.3,2f7.3,f9.4,f7.4)')
*    &         k,zgrid(i,j,k),zref,
*    &         uold(k),uref,vold(k),vref,
*    &         -g*thref*(thold(k)+thbase),bref,
*    &         ritop(k),dvsq(k)
*           call flush(lp)
*         endif
c
          if (-zgrid(i,j,k).lt.10.0) then
            swfrac(k)=redfac(jrlv)*
     &                exp(zgrid(i,j,k)*onem*betard(jrlv))+
     &                (1.-redfac(jrlv))*
     &                exp(zgrid(i,j,k)*onem*betabl(jrlv))
          else
            swfrac(k)=(1.-redfac(jrlv))*
     &                exp(zgrid(i,j,k)*onem*betabl(jrlv))
          endif
        enddo  !k=1,klist
c
c --- calculate interface variables required to estimate interior diffusivities
        do k=1,klist(i,j)
          ka=k+1
          shsq  (ka)=(uold(k)-uold(ka))**2+(vold(k)-vold(ka))**2
          alfadt(ka)=.5*(dsigdt(told(k ),sold(k ))+
     &                   dsigdt(told(ka),sold(ka)))*(told(k)-told(ka))
          betads(ka)=.5*(dsigds(told(k ),sold(k ))+
     &                   dsigds(told(ka),sold(ka)))*(sold(k)-sold(ka))
          dbloc (ka)=-g*thref*(thold(k)-thold(ka))
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
        do k=2,klist(i,j)+1
          if (shinst) then
            rigr=max(0.,dbloc(k)*(zgrid(i,j,k-1)-zgrid(i,j,k))/
     &          (shsq(k)+epsil))
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
        enddo 
c
c --- double-diffusion (salt fingering and diffusive convection)
        if (dbdiff) then
          do k=2,klist(i,j)+1
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
     &         .and. -alfadt(k).gt.betads(k)) then
              rrho=-alfadt(k)/betads(k)
              diffdd=1.5e-6*9.*.101*exp(4.6*exp(-.54*(1./rrho-1.)))
              if (rrho.gt.0.5) then
                prandtl=(1.85-.85/rrho)*rrho
              else
                prandtl=.15*rrho
              endif
              dift(i,j,k)=dift(i,j,k)+diffdd
              difs(i,j,k)=difs(i,j,k)+prandtl*diffdd
            endif
          enddo
        endif
c
cdiag   if (i.eq.itest.and.j.eq.jtest) then
cdiag      write (lp,102) (nstep,iter,i+i0,j+j0,k,
cdiag&     hwide(k),1.e4*vcty(i,j,k),1.e4*dift(i,j,k),1.e4*difs(i,j,k),
cdiag&       k=1,kk+1)
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
c --- minimum hbl is hwide() + 1 cm, maximum is bottom mid-layer - 1 cm.
c
        hblmin=hwide(1)+0.01
        hblmax=-zgrid(i,j,klist(i,j))-0.01
c
c --- buoyfl = buoyancy flux (m**2/sec**3) into atmos.
c --- note: surface density increases (column is destabilized) if buoyfl > 0
c --- buoysw = shortwave radiation buoyancy flux (m**2/sec**3) into atmos.
c --- salflx, sswflx and surflx are positive into the ocean
      tmn=.5*(temp(i,j,1,m)+temp(i,j,1,n))
      smn=.5*(saln(i,j,1,m)+saln(i,j,1,n))
      dsgdt=dsigdt(tmn,smn)
      buoyfl=g*thref*(dsigds(tmn,smn)*salflx(i,j)*thref+
     &                dsgdt          *surflx(i,j)*thref/spcifh)
      buoysw=g*thref*(dsgdt          *sswflx(i,j)*thref/spcifh)
c 
c --- diagnose the new boundary layer depth as the depth where a bulk
c --- richardson number exceeds ric
c
c --- initialize hbl and nbl to bottomed out values
        k1=1
        k2=2
        rib(k1)=0.0
        nbl=klist(i,j)
        hbl=hblmax
c
c --- diagnose hbl and nbl
        do k=2,nbl
          case=-zgrid(i,j,k)
          bfsfc=buoyfl-swfrac(k)*buoysw
          stable=.5+sign(.5,-bfsfc)
          dnorm=stable+(1.-stable)*epsilon
c
c --- compute turbulent velocity scales at dnorm, for
c --- hbl = case = -zgrid(i,j,k)
          call wscale(i,j,case,dnorm,bfsfc,wm,ws)
c
c --- compute the turbulent shear contribution to rib
          bvsq=.5*(dbloc(k  )/(zgrid(i,j,k-1)-zgrid(i,j,k  ))+
     &             dbloc(k+1)/(zgrid(i,j,k  )-zgrid(i,j,k+1)))
          vtsq=-zgrid(i,j,k)*ws*sqrt(abs(bvsq))*vtc
c
c --- compute bulk richardson number at new level
c --- linearly interpolate to find hbl as the depth where rib = ricr
          rib(k2)=ritop(k)/(dvsq(k)+vtsq+epsil)
          if (nbl.eq.klist(i,j).and.rib(k2).ge.ricr) then
            hbl=-zgrid(i,j,k-1)+(zgrid(i,j,k-1)-zgrid(i,j,k))*
     &             (ricr-rib(k1))/(rib(k2)-rib(k1)+epsil)
            nbl=k
            if (hbl.lt.hblmin) then
              hbl=hblmin
              nbl=2
            endif
            if (hbl.gt.hblmax) then
              hbl=hblmax
              nbl=klist(i,j)
            endif
          endif
c
          ksave=k1
          k1=k2
          k2=ksave
        enddo  !k=1,nbl
c
c --- calculate swfrml, the fraction of solar radiation absorbed by depth hbl
        k=nbl
        q=(zgrid(i,j,k-1)+hbl)/(zgrid(i,j,k-1)-zgrid(i,j,k))
        swfrml=swfrac(k-1)+q*(swfrac(k)-swfrac(k-1))
c
c --- limit check on hbl for negative (stablizing) surface buoyancy forcing
        bfsfc=buoyfl-swfrml*buoysw
        stable=.5+sign(.5,-bfsfc)
        bfsfc=bfsfc-stable*epsil                 !insures bfsfc never=0
        if (bfsfc.lt.0.) then
          hmonob=-cmonob*ustar(i,j)**3/(vonk*bfsfc)
          hbl=max(hblmin,
     &            min(hbl,
     &                hekman(i,j),
     &                hmonob,
     &                hblmax))
        endif
c
c--- find new nbl and re-calculate swfrml
        nbl=klist(i,j)
        do k=2,klist(i,j)
          if (-zgrid(i,j,k).gt.hbl) then
            nbl=k
            q=(zgrid(i,j,k-1)+hbl)/(zgrid(i,j,k-1)-zgrid(i,j,k))
            swfrml=swfrac(k-1)+q*(swfrac(k)-swfrac(k-1))
            exit
          endif
        enddo
c
c --- find forcing stability and buoyancy forcing for final hbl values
c --- determine case (for case=0., hbl lies between -zgrid(i,j,nbl)
c --- and the interface above. for case=1., hbl lies between 
c --- -zgrid(i,j,nbl-1) and the interface below)
c
c --- velocity scales at hbl
        bfsfc=buoyfl-swfrml*buoysw
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
        q=(hbl*onem-p(i,j,ka))*qdpmm(ka)
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
cdiag&    hwide(k),1.e4*vcty(i,j,k),1.e4*dift(i,j,k),1.e4*difs(i,j,k),
cdiag&      ghats(i,j,k),k=1,kk+1)
cdiag     call flush(lp)
cdiag   endif
c
c --- save array dpbl=onem*hbl for ice, output and diagnosis
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
c --- salflx, sswflx and surflx are positive into the ocean
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
cdiag     if (i.eq.itest.and.j.eq.jtest) then
cdiag       write (lp,104) (nstep,iter,i+i0,j+j0,k,
cdiag&        hm(k),t1do(k),t1dn(k),s1do(k),s1dn(k),
cdiag&        0.0,0.0,
cdiag&        k=1,nlayer)
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
cdiag&        hm(k),u1do(k),u1dn(k),v1do(k),v1dn(k),k=1,nlayer)
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
     &     /(i9,i2,3i4,2x,4f10.2))
 103  format(25x,'   thick      viscty    t diff    s diff   nonlocal'
     &     /(i9,i2,3i4,2x,4f10.2,f11.6))
 104  format(25x,
     &     '  thick   t old   t new   s old   s new trc old trc new'
     &     /(i9,i2,3i4,1x,f9.2,4f8.3,2f7.4))
 105  format(25x,'   thick   u old   u new   v old   v new'
     &     /(i9,i2,3i4,1x,f10.2,4f8.3))
c
      return
      end
c
      subroutine mxmyaij(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j
c
c ---------------------------------------------------------------
c --- mellor-yamada 2.5 vertical diffusion, single j-row (part A)
c --- vertical coordinate is z negative below the ocean surface
c ---------------------------------------------------------------
c
c --- arrays q2 and q2l are prognostic variables representing tke and
c --- tke multiplied by the turbulent eddy length scale. these arrays
c --- are calculated on interfaces in a special vertical grid where
c --- interfaces are centered at the surface, bottom, and at mid-depths
c --- of each hycom layer (kdm+2 interfaces). this enables the q2 and
c --- q2l arrays to be advected and diffused in subroutine tsadvc.
c
c --- local variables for my2.5 mixing
c
      real sm(kdm+2),sh(kdm+2),prod(kdm+2),stf(kdm+2)
      real dtef(kdm+2),gh(kdm+2),ee(kdm+2),gg(kdm+2),turlen(kdm+2)
      real thstr(kdm),sfht(kdm),z(kdm+2),zz(kdm+2),dz(kdm+2),dzz(kdm+2)
      real th1d(kdm+1),u1d(kdm+1),v1d(kdm+1)
c
      real dbloc(kdm+2)        ! buoyancy jump across interface
      real swfrac(kdm+1)       ! fractional surface shortwave radiation flux
      real dpmm(kdm)           !     max(onemm,dp(i,j,:,[nm]))
      real qdpmm(kdm)          ! 1.0/max(onemm,dp(i,j,:,[nm]))
      real pij(kdm+1)          ! local copy of p(i,j,:)
c
      real dh,oneeta,akn,coef1,coef2,coef3,dtemp,dsaln,sflux1,pmid,div
      real wusurf,wvsurf,wubot,wvbot,delu,delv,ubav,vbav,qspcifh
c
      integer k,k1,khy1,khy2,kmy1,kmy2,jrlv
c
      include 'stmt_fns.h'
c
c --- set mid-time pressure array
c --- locate lowest substantial mass-containing layer.
      pij(1)=p(i,j,1)
      do k=1,kk
        dpmm( k)  =max(onemm,dp(i,j,k,m))
        pij(  k+1)=pij(k)+dp(i,j,k,m)
        p(i,j,k+1)=pij(k+1)
      enddo
      do k=kk,1,-1
        if (dpmm(k).gt.tencm) then
          exit
        endif
      enddo
      klist(i,j)=max(k,2)  !always consider at least 2 layers
c
c --- calculate surface elevation anomaly and dh:
      if (tbaric) then
        thstr(kk)=th3d(i,j,kk,m)
     &           +kappaf(temp(i,j,kk,m),saln(i,j,kk,m),pij(kk+1))
      else
        thstr(kk)=th3d(i,j,kk,m)
      endif
      sfht(kk)=(psikk(i,j)+
     &         (pij(kk+1)*(thkk(i,j)-thstr(kk))
     &         -pbavg(i,j,m)*(thstr(kk)+thbase))*thref**2)/(thref*onem)
      do k=kk-1,1,-1
        if (tbaric) then
          thstr(k)=th3d(i,j,k,m)
     &            +kappaf(temp(i,j,k,m),saln(i,j,k,m),pij(k))
        else
          thstr(k)=th3d(i,j,k,m)
        endif
        oneeta=1.+pbavg(i,j,m)/pij(kk+1)
        sfht(k)=sfht(k+1)+(pij(k+1)*oneeta
     &        *(thstr(k+1)-thstr(k))*thref)*qonem
      enddo
      dh=depths(i,j)+sfht(1)
c
c --- generate the scaled m-y vertical grid
c --- calculate z (interface z), dz, zz (central layer z), dzz in m
      khy1=klist(i,j)
      khy2=khy1+1
      kmy1=khy2
      kmy2=kmy1+1
c
      z(1)=0.0
      dz(1)=0.5*dp(i,j,1,m)/p(i,j,khy2)
      z(2)=-dz(1)
      do k=3,kmy1
        dz(k-1)=0.5*(dp(i,j,k-2,m)+dp(i,j,k-1,m))/pij(khy2)
        z(k)=z(k-1)-dz(k-1)
      enddo
      dz(kmy1)=0.5*dp(i,j,khy1,m)/pij(khy2)
      z(kmy2)=z(kmy1)-dz(kmy1)
      dz(kmy2)=0.0
c
      do k=1,kmy1
        zz(k)=0.5*(z(k)+z(k+1))
      enddo
      zz(kmy2)=2.0*z(kmy2)-zz(kmy1)
c
      do k=1,kmy1
        dzz(k)=zz(k)-zz(k+1)
      enddo
      dzz(kmy2)=0.0
c
c --- calculate 1-d arrays on m-y vertical grid
      th1d(1)=th3d(i,j,1,m)
      ubav=0.5*(ubavg(i,j,m)+ubavg(i+1,j  ,m))
      vbav=0.5*(vbavg(i,j,m)+vbavg(i  ,j+1,m))
      u1d(1)=0.5*(u(i,j,1,m)+u(i+1,j  ,1,m))+ubav
      v1d(1)=0.5*(v(i,j,1,m)+v(i  ,j+1,1,m))+vbav
      do k=2,khy1
        th1d(k)=0.5*(th3d(i,j,k-1,m)+th3d(i,j,k,m))
        u1d(k)=0.25*(u(i  ,j  ,k-1,m)+u(i  ,j  ,k,m)
     &              +u(i+1,j  ,k-1,m)+u(i+1,j  ,k,m))+ubav
        v1d(k)=0.25*(v(i  ,j  ,k-1,m)+v(i  ,j  ,k,m)
     &              +v(i  ,j+1,k-1,m)+v(i  ,j+1,k,m))+vbav
      enddo
      if (khy1.lt.kk) then
        th1d(kmy1)=max(th1d(khy1),
     &             0.5*(th3d(i,j,khy1,m)+th3d(i,j,kmy1,m)))
      else
        th1d(kmy1)=max(th1d(khy1),th3d(i,j,khy1,m))
      endif
      u1d(kmy1)=0.5*(u(i,j,khy1,m)+u(i+1,j  ,khy1,m))+ubav
      v1d(kmy1)=0.5*(v(i,j,khy1,m)+v(i  ,j+1,khy1,m))+vbav
c
c --- make sure background tke maintains minimum value
      do k=0,kk+1
        q2( i,j,k,m)=max(smll,q2( i,j,k,m))
        q2l(i,j,k,m)=max(smll,q2l(i,j,k,m))
        q2( i,j,k,n)=max(smll,q2( i,j,k,n))
        q2l(i,j,k,n)=max(smll,q2l(i,j,k,n))
      enddo
c
c --- calculate sm,sh coefficients
      do k=2,kmy1
        sm(k)=-delt1*0.5*(difqmy(i,j,k-1)+difqmy(i,j,k  )+2.0*difsiw)
     &       /(dzz(k-1)*dz(k  )*dh*dh)
        sh(k)=-delt1*0.5*(difqmy(i,j,k-2)+difqmy(i,j,k-1)+2.0*difsiw)
     &       /(dzz(k-1)*dz(k-1)*dh*dh)
      enddo
c
c ------------------------------------------------------------------
c --- solve delt1*(difq*q2(n)')' - q2(n)*(2.*delt1*dtef+1.) = -q2(m)
c --- for q2(n) (tke)
c ------------------------------------------------------------------
c
c --- surface and bottom stress boundary conditions
      if (windf) then
        wusurf=thref*(taux(i,j,l0)*w0+taux(i,j,l1)*w1
     &               +taux(i,j,l2)*w2+taux(i,j,l3)*w3)
        wvsurf=thref*(tauy(i,j,l0)*w0+tauy(i,j,l1)*w1
     &               +tauy(i,j,l2)*w2+tauy(i,j,l3)*w3)
      else
        wusurf=0.0
        wvsurf=0.0
      endif
      wubot=-0.5*thkbot*onem*u1d(kmy1)*drag(i,j)*thref/g
      wvbot=-0.5*thkbot*onem*v1d(kmy1)*drag(i,j)*thref/g
      ee(1)=0.0
      gg(1)=const1*sqrt(wusurf*wusurf+wvsurf*wvsurf)
      q2(i,j,kmy1,n)=const1*sqrt(wubot*wubot+wvbot*wvbot)
c
c --- calculate vertical buoyancy gradient at interfaces
      dbloc(1)=0.0
      do k=2,kmy1
        dbloc(k)=thref*g*(th1d(k-1)-th1d(k))/(dzz(k-1)*dh)
      enddo
      dbloc(kmy2)=0.0
c
c --- calculate turbulent length scale and richardson number gh at interfaces
c
      turlen(1)=0.0
      gh(1)=0.
      do k=2,kmy1
        turlen(k)=q2l(i,j,k-1,m)/q2(i,j,k-1,m)
        gh(k)=min(0.028,turlen(k)**2/q2(i,j,k-1,m)*dbloc(k))
      enddo
      turlen(kmy2)=0.   
      gh(kmy2)=0.
c
c --- calculate tke production at interfaces
      prod(1)=0.0
      do k=2,kmy1
        delu=u1d(k)-u1d(k-1)
        delv=v1d(k)-v1d(k-1)
        prod(k)=diftmy(i,j,k-1)*dbloc(k)+0.5*vctymy(i,j,k-1)*sef
     &         *(delu**2+delv**2)/(dzz(k-1)*dh)**2
      enddo
      prod(kmy2)=0.0
c
c --- solve the equation
      do k=1,kmy2
        stf(k)=1.0
        if (gh(k).lt.0. ) stf(k)=1.0-0.9*(gh(k)/ghc)**1.5
        if (gh(k).lt.ghc) stf(k)=0.1
        dtef(k)=q2(i,j,k-1,m)**1.5/(b1my*q2l(i,j,k-1,m)+smll)*stf(k)
      enddo
      do k=2,kmy1
        gg(k)=1.0/(sm(k)+sh(k)*(1.0-ee(k-1))-(2.0*delt1*dtef(k)+1.0))
        ee(k)=sm(k)*gg(k)
        gg(k)=(-2.0*delt1*prod(k)+sh(k)*gg(k-1)-q2(i,j,k-1,n))*gg(k)
      enddo
      do k=kmy1,1,-1
        q2(i,j,k-1,n)=ee(k)*q2(i,j,k,n)+gg(k)
      enddo
c
c ------------------------------------------------------------------
c --- solve delt1*(difq*q2l(n)')' - q2l(n)*(delt1*dtef+1.) = -q2l(m)
c --- for q2l(n) (tke times turbulent length scale)
c ------------------------------------------------------------------
c
      q2(i,j,1,n)=max(smll,q2(i,j,1,n))
      ee(2)=0.0
      gg(2)=-vonk*z(2)*dh*q2(i,j,1,n)
      q2l(i,j,kmy1,n)=0.0
      do k=3,kmy1
        dtef(k)=dtef(k)*(1.+e2my*((1.0/abs(z(k)-z(1))+
     &         1.0/abs(z(k)-z(kmy2)))*turlen(k)/(dh*vonk))**2)
        gg(k)=1.0/(sm(k)+sh(k)*(1.0-ee(k-1))-(delt1*dtef(k)+1.0))
        ee(k)=sm(k)*gg(k)
        gg(k)=(delt1*(-prod(k)*turlen(k)*e1my)+sh(k)*gg(k-1)
     &       -q2l(i,j,k-1,n))*gg(k)
      enddo
      do k=kmy1,2,-1
        q2l(i,j,k-1,n)=ee(k)*q2l(i,j,k,n)+gg(k)
      enddo
      do k=0,kmy1
        q2( i,j,k,n)=max(smll,q2( i,j,k,n))
        q2l(i,j,k,n)=max(smll,q2l(i,j,k,n))
      enddo
c
c ----------------------------------------------------
c --- calculate the viscosity and diffusivity profiles
c ----------------------------------------------------
c
c --- note that sm and sh limit to infinity when gh approaches 0.0288
      do k=1,kmy2
        coef1=a2my*(1.0-6.0*a1my/b1my*stf(k))
        coef2=3.0*a2my*b2my/stf(k)+18.0*a1my*a2my
        coef3=a1my*(1.0-3.0*c1my-6.0*a1my/b1my*stf(k))
        sh(k)=coef1/(1.0-coef2*gh(k))
        sm(k)=coef3+sh(k)*coef4*gh(k)
        sm(k)=sm(k)/(1.0-coef5*gh(k))
        akn=turlen(k)*sqrt(abs(q2(i,j,k-1,n)))
        difqmy(i,j,k-1)=max(difsiw,(akn*0.41*sh(k)+difqmy(i,j,k-1))*0.5)
        vctymy(i,j,k-1)=max(difmiw,(akn*     sm(k)+vctymy(i,j,k-1))*0.5)
        diftmy(i,j,k-1)=max(difsiw,(akn*     sh(k)+diftmy(i,j,k-1))*0.5)
      enddo
c
c --- set diffusivity/viscosty on hycom interfaces
c --- save surface boundary layer array (dpbl) for ice, output and diagnosis
      vcty(i,j,1)=0.0
      dift(i,j,1)=0.0
      difs(i,j,1)=0.0
      dpbl(i,j)  =0.0
      do k=2,khy1
        vcty(i,j,k)=0.5*(vctymy(i,j,k-1)+vctymy(i,j,k))
        dift(i,j,k)=0.5*(diftmy(i,j,k-1)+diftmy(i,j,k))
        difs(i,j,k)=dift(i,j,k)
        if     (dpbl(i,j).eq.0.0 .and.
     &          dift(i,j,k).lt.20.0*difsiw) then
          dpbl(i,j) = p(i,j,k+1)
        endif
      enddo
      if     (dpbl(i,j).eq.0.0) then
        dpbl(i,j) = p(i,j,khy1+1)
      endif
      vcty(i,j,khy2)=0.0
      dift(i,j,khy2)=0.0
      difs(i,j,khy2)=0.0
c
c --- set new time pressure array
c --- locate lowest substantial mass-containing layer.
      pij(1)=p(i,j,1)
      do k=1,kk
        dpmm( k)  =max(onemm,dp(i,j,k,n))
        qdpmm(k)  =1.0/dpmm(k)
        pij(  k+1)=pij(k)+dp(i,j,k,n)
        p(i,j,k+1)=pij(k+1)
      enddo
      do k=kk,1,-1
        if (dpmm(k).gt.tencm) then
          exit
        endif
      enddo
      klist(i,j)=max(k,2)  !always consider at least 2 layers
c
      if (klist(i,j).gt.khy1) then
        do k=khy1,klist(i,j)-1
          vcty(i,j,k+1)=vcty(i,j,k)
          dift(i,j,k+1)=dift(i,j,k)
          difs(i,j,k+1)=difs(i,j,k)
        enddo
        khy2=klist(i,j)+1
        vcty(i,j,khy2)=0.0
        dift(i,j,khy2)=0.0
        difs(i,j,khy2)=0.0
      endif
      khy1=klist(i,j)
      khy2=klist(i,j)+1
c
cdiag if (i.eq.itest.and.j.eq.jtest) then
cdiag    write (lp,101) (nstep,i+i0,j+j0,k,pij(k)*qonem,
cdiag.   1.e4*vcty(i,j,k),1.e4*dift(i,j,k),1.e4*difs(i,j,k),
cdiag.   1.e4*difqmy(i,j,k-1),k=1,khy2)
cdiag    call flush(lp)
cdiag endif
c
c --- calculate zgrid
      zgrid(i,j,1)=-0.5*dp(i,j,1,n)*qonem
      do k=2,khy1
        zgrid(i,j,k)=zgrid(i,j,k-1)
     &              -0.5*(dp(i,j,k-1,n)+dp(i,j,k,n))*qonem
      enddo
      zgrid(i,j,khy2)=-pij(khy2)*qonem
c
c --- forcing of t,s by surface fluxes. shortwave flux penetration depends on
c --- jerlov water type, which is set in thermf.f.  flux positive into ocean
c
      jrlv=jerlov(i,j)
      qspcifh=1.0/spcifh
c
      do k=1,khy1
        if (thermo) then
          if (k.eq.1) then
            if (pij(k+1).lt.tenm) then
              swfrac(k+1)=redfac(jrlv) *exp(-pij(k+1)*betard(jrlv))+
     &                (1.-redfac(jrlv))*exp(-pij(k+1)*betabl(jrlv))
            else
              swfrac(k+1)=(1.-redfac(jrlv))*
     &                    exp(-pij(k+1)*betabl(jrlv))
            endif
            sflux1=surflx(i,j)-sswflx(i,j)
            dtemp=(sflux1+(1.-swfrac(k+1))*sswflx(i,j))*
     &            delt1*g*qspcifh*qdpmm(k)
            dsaln=salflx(i,j)*
     &            delt1*g*        qdpmm(k)
cdiag       if (i.eq.itest .and. j.eq.jtest) then
cdiag         write (lp,102) nstep,i+i0,j+j0,k, 
cdiag&          0.,1.-swfrac(k+1),dtemp,dsaln
cdiag         call flush(lp)
cdiag       endif
          else
            if (pij(k+1).lt.tenm) then
              swfrac(k+1)=redfac(jrlv) *exp(-pij(k+1)*betard(jrlv))+
     &                (1.-redfac(jrlv))*exp(-pij(k+1)*betabl(jrlv))
            else
              swfrac(k+1)=(1.-redfac(jrlv))*
     &                    exp(-pij(k+1)*betabl(jrlv))
            endif
            if(k.ne.khy1) then
              dtemp=(swfrac(k)-swfrac(k+1))*sswflx(i,j)*
     &              delt1*g*qspcifh*qdpmm(k)
            else
              dtemp= swfrac(k)             *sswflx(i,j)*
     &               delt1*g*qspcifh*qdpmm(k)
            endif
            dsaln=0.
cdiag       if (i.eq.itest .and. j.eq.jtest) then
cdiag          write (lp,102) nstep,i+i0,j+j0,k,
cdiag&         1.-swfrac(k),1.-swfrac(k+1),dtemp
cdiag         call flush(lp)
cdiag       endif
          endif
        else !.not.thermo
          dtemp=0.0
          dsaln=0.0
        endif
        temp(i,j,k,n)=temp(i,j,k,n)+dtemp
        saln(i,j,k,n)=saln(i,j,k,n)+dsaln
        th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
      enddo
c
 101  format(25x,'   thick      viscty    t diff    s diff    q diff  '
     &     /(i9,3i4,2x,5f10.2))
 102  format(i9,3i4,'absorbup,dn,dtemp,dsaln ',2f6.3,2f10.6)
c
      return   
      end
c
      subroutine mxkprfbij(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j
c
c -------------------------------------------------------------
c --- k-profile vertical diffusion, single j-row (part B)
c --- vertical coordinate is z negative below the ocean surface
c -------------------------------------------------------------
c
c --- perform the final vertical mixing at p points
c
c --- local 1-d arrays for matrix solution
      real t1do(kdm+1),t1dn(kdm+1),s1do(kdm+1),s1dn(kdm+1),
     &     tr1do(kdm+1,mxtrcr),tr1dn(kdm+1,mxtrcr),
     &     difft(kdm+1),diffs(kdm+1),
     &     ghat(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     &     tcc(kdm),         ! central ...     (k  ) ..
     &     tcl(kdm),         ! lower .....     (k-1) ..
     &     rhs(kdm)          ! right-hand-side terms
c
      real    ghatflux
      integer k,ka,ktr,nlayer
c
      include 'stmt_fns.h'
c
      nlayer=klist(i,j)
c
      do k=1,nlayer
        difft(k+1)=dift(i,j,k+1)
        diffs(k+1)=difs(i,j,k+1)
        ghat(k+1)=ghats(i,j,k+1)
        t1do(k)=temp(i,j,k,n)
        s1do(k)=saln(i,j,k,n)
        do ktr= 1,ntracr
          tr1do(k,ktr)=tracer(i,j,k,n,ktr)
        enddo
        hm(k)=max(onemm,dp(i,j,k,n))*qonem
        zm(k)=zgrid(i,j,k)
      enddo
c
      k=nlayer+1
      ka=min(k,kk)
      difft(k)=0.0
      diffs(k)=0.0
      ghat(k)=0.0
      t1do(k)=temp(i,j,ka,n)
      s1do(k)=saln(i,j,ka,n)
      do ktr= 1,ntracr
        tr1do(k,ktr)=tracer(i,j,ka,n,ktr)
      enddo
      zm(k)=zm(k-1)-0.001
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
c --- salflx, sswflx and surflx are positive into the ocean
c
c --- t solution
      ghatflux=-(surflx(i,j)-sswflx(i,j))*thref/spcifh
      call tridcof(difft,tri,nlayer,tcu,tcc,tcl)
      call tridrhs(hm,t1do,difft,ghat,ghatflux,tri,nlayer,rhs)
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,t1do,t1dn,difft)
c
c --- t-like tracer solution
      do ktr= 1,ntracr
        if     (trcflg(ktr).eq.2) then
          ghatflux=-(surflx(i,j)-sswflx(i,j))*thref/spcifh
          call tridrhs(hm,
     &                 tr1do(1,ktr),difft,ghat,ghatflux,tri,nlayer,rhs)
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,
     &                 tr1do(1,ktr),tr1dn(1,ktr),difft)
        endif
      enddo
c
c --- s solution
      ghatflux=-salflx(i,j)*thref
      call tridcof(diffs,tri,nlayer,tcu,tcc,tcl)
      call tridrhs(hm,s1do,diffs,ghat,ghatflux,tri,nlayer,rhs)
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,s1do,s1dn,diffs)
c
c --- standard tracer solution
      do ktr= 1,ntracr
        if     (trcflg(ktr).ne.2) then
          ghatflux=0.
          call tridrhs(hm,
     &                 tr1do(1,ktr),diffs,ghat,ghatflux,tri,nlayer,rhs)
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,
     &                 tr1do(1,ktr),tr1dn(1,ktr),diffs)
        endif
      enddo
c
c --- adjust t, s, th, arrays
      do k=1,klist(i,j)
        temp(i,j,k,n)=t1dn(k)
        saln(i,j,k,n)=s1dn(k)
        th3d(i,j,k,n)=sig(t1dn(k),s1dn(k))-thbase
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr)=tr1dn(k,ktr)
        enddo
      enddo
c
      return
      end
c
      subroutine mxkprfciju(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j
c
c -------------------------------------------------------------------------
c --- k-profile vertical diffusion, single j-row, momentum at u grid points
c --- vertical coordinate is z negative below the ocean surface
c -------------------------------------------------------------------------
c
c local variables for kpp mixing
c
c --- local 1-d arrays for matrix solution
      real u1do(kdm+1),u1dn(kdm+1),
     &     diffm(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     &     tcc(kdm),         ! central ...     (k  ) ..
     &     tcl(kdm),         ! lower .....     (k-1) ..
     &     rhs(kdm)          ! right-hand-side terms
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
          hm(k)=max(onemm,dpu(i,j,k,n))*qonem
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
cdiag&    hm(k),u1do(k),u1dn(k),k=1,nlayer)
cdiag   call flush(lp)
cdiag endif
      return
 106  format(23x,'   thick   u old   u new'/(i9,3i4,1x,f10.3,2f8.3))
      end
      subroutine mxkprfcijv(m,n, i,j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, i,j
c
c --------------------------------------------------------------------------
c --- k-profile vertical diffusion, single j-row , momentum at v grid points
c --- vertical coordinate is z negative below the ocean surface
c --------------------------------------------------------------------------
c
c local variables for kpp mixing
c
c --- local 1-d arrays for matrix solution
      real v1do(kdm+1),v1dn(kdm+1),
     &     diffm(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
c
c --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),         ! upper coeff for (k-1) on k line of trid.matrix
     &     tcc(kdm),         ! central ...     (k  ) ..
     &     tcl(kdm),         ! lower .....     (k-1) ..
     &     rhs(kdm)          ! right-hand-side terms
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
          hm(k)=max(onemm,dpv(i,j,k,n))*qonem
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
cdiag&    hm(k),v1do(k),v1dn(k),k=1,nlayer)
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
c --- vertical coordinate is z negative below the ocean surface
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
     &        wam,wbm,was,wbs,ucube,zehat
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
c> May  2002 - buoyfl (into the atmos.), calculated here
