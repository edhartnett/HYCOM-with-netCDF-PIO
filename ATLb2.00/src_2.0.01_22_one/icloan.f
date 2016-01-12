      subroutine icloan(m,n)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- 'energy loan' ice model. no advection, no dynamics. ice amount
c --- represents energy 'loaned' to water column to prevent wintertime
c --- cooling below freezing level. loan is paid back in summer.
c
      real, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & icex
c
      integer i,j,l
      real    tmxl,paybak,borrow,hice,thicmx
c
cdiag real*8  d1
cdiag real    qmax
cdiag integer imax,jmax
c
c --- thkice = average ice thickness ( = hice x covice(i,j))
c --- covice = ice coverage (rel.units, 0.0 to 1.0)
c --- temice = ice surface temperature
c --- surflx = net total heatflux between atm and ocean or ice (W/m^2)
c --- salflx = salt flux (implied by fresh water flux)
c
c --- hice   = actual ice thickness (m)
c
c --- tice   = melting point (deg)
c --- tgrad  = vert. temperature gradient inside ice (deg/m)
c --- thin   = min. ice thickness (m)
c --- thkmax = max. ice thickness (m)
c --- rhoice = density of ice (kg/m**3)
c --- fusion = latent heat of fusion (J/kg)
c --- rate   = max. ice melting rate (m/sec)
c --- saldif = salinity difference water minus ice
c
      real       tice,tgrad,thin,thkmax,rhoice,fusion,rate,saldif
      parameter (tice=-1.8, tgrad=2.0, thin=0.5, thkmax=10.0)
      parameter (rhoice=917.0, fusion=334.e3, rate=5.e-6, saldif=25.0)
c
c --- energy loan: add extra energy to the ocean to keep SST from dropping
c --- below tice in winter. return this borrowed energy to the 'energy bank'
c --- in summer as quickly as surflx < 0 allows.
c
c --- salt loan: analogous to energy loan.
c
      margin = 0  ! no horizontal derivatives
c
      thicmx=0.0
!$OMP PARALLEL DO PRIVATE(j,l,i,tmxl,borrow,paybak)
!$OMP&            REDUCTION(MAX:thicmx)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 10 j=1-margin,jj+margin
      do 11 l=1,isp(j)
      do 11 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
cdiag util1(i,j)=thkice(i,j)*scp2(i,j)
c
c --- calculate hypothetical mixed-layer temp due to diab. forcing 
      tmxl=temp(i,j,1,n)+surflx(i,j)*delt1*g/(spcifh*dpbl(i,j)) 
c
      borrow=(tice-tmxl)*spcifh*dpbl(i,j)/(delt1*g)  ! W/m^2
      if (tmxl.lt.tice) then
c
c --- add energy to bring tmxl back to tice (only if tmxl < tice)
c
cdiag   if (thkice(i,j).le.0.) write (lp,'(i9,a,2i5)') nstep,
cdiag.    '  new ice forms at i,j=',i,j
        surflx(i,j)=surflx(i,j)+borrow
        salflx(i,j)=salflx(i,j)+saldif*borrow/fusion
        thkice(i,j)=thkice(i,j)+borrow*delt1/(fusion*rhoice)
c
cdiag   if (i.eq.itest .and. j.eq.jtest)
cdiag.   write (lp,'(i9,2i5,a,4f9.3)') nstep,i,j,'  sst,tmxl,ice,loan:',
cdiag.    temp(i,j,1,n),tmxl,thkice(i,j),borrow
c
c --- return the borrowed amount whenever tmxl > tice
c
      else if (thkice(i,j).gt.0.) then  ! tmxl > tice
c
        paybak=min(thkice(i,j)*fusion*rhoice/delt1,-borrow,
     &                                      rate*fusion*rhoice)
        surflx(i,j)=surflx(i,j)-paybak
        salflx(i,j)=salflx(i,j)-saldif*paybak/fusion
        thkice(i,j)=thkice(i,j)-paybak*delt1/(fusion*rhoice)
c
cdiag   if (thkice(i,j).le.0.) write (lp,'(i9,a,2i5)') nstep,
cdiag.     '  ice depleted at  i,j=',i,j
c
cdiag   if (i.eq.itest .and. j.eq.jtest)
cdiag.   write (lp,'(i9,2i5,a,4f9.3)') nstep,i,j,'  sst,tmxl,ice,loan:',
cdiag.    temp(i,j,1,n),tmxl,thkice(i,j),-paybak
c
      end if
c
      icex(i,j)=max(thkice(i,j)-thkmax,0.)  ! ice exceeding thkmax
      thicmx=max(thicmx,thkice(i,j))
 11   continue
 10   continue
!$OMP END PARALLEL DO
c
      call xcmaxr(thicmx)
c
cdiag call xcsum(d1, util1,ip)
cdiag if (d1.gt.0.0 .and. mnproc.eq.1) then
cdiag   write(lp,'(a,2f11.4)') 'total ice:',d1/area
cdiag   call flush(lp)
cdiag endif
c
c --- spread out portion of ice thicker than thkmax
      if (thicmx.gt.thkmax) then
        call psmooth(icex, margin)
      endif
c
!$OMP PARALLEL DO PRIVATE(j,l,i,hice)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            thkice(i,j)=icex(i,j)+min(thkice(i,j),thkmax)
c
            if (thkice(i,j).lt.1.e-5*thin) then
              covice(i,j)=0.
              temice(i,j)=temp(i,j,1,n)
            else
c --- compute fractional ice coverage for energy flux calculation
              covice(i,j)=min(1.,thkice(i,j)/thin)
c
c --- compute ice surface temperature (minimum -80 deg)
              hice=thkice(i,j)/covice(i,j)
              temice(i,j)=max(-80.,tice-tgrad*hice)
            end if
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
C
cdiag imax=0
cdiag jmax=0
cdiag qmax=0
cdiag do 20 j=1-margin,jj+margin
cdiag do 20 l=1,isp(j)
cdiag do 20 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
cdiag if (thkice(i,j).gt.qmax) then
cdiag   qmax=thkice(i,j)
cdiag   imax=i
cdiag   jmax=j
cdiag end if
 20   continue
cdiag if (qmax.gt.0.) write (lp,'(a,2i5,f11.3)') 'max. ice thknss at',
cdiag.   imax,jmax,qmax
c
      return
      end subroutine icloan
c
c
c> Revision history
c>
c> June 2000 - conversion to SI units
c> July 2000 - switched sign convention for vertical fluxes (now >0 if down)
