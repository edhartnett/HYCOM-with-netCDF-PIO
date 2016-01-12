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
c --- presence of ice can effect the atmos. heat flux, see thermf.
c
      real, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & icex
c
      integer i,j,l
      real    tfrz,tsur,tmxl,smxl,paybak,borrow,hice,thicmx,f2t
c
cdiag real*8  d1
cdiag real    qmax
cdiag integer imax,jmax
c
c --- hice   = actual ice thickness (m), local variable
c
c --- thkice = average ice thickness ( = hice x covice(i,j))
c --- covice = ice coverage (rel. units, 0.0 to 1.0)
c --- temice = ice surface temperature (degC)
c --- surflx = net total heatflux between atm and ocean or ice (W/m^2)
c --- salflx = salt flux (implied by fresh water flux)
c
      real       tfrz_n,ticemn,ticemx,salice,rhoice,fusion,meltmx
      parameter (tfrz_n= -1.79, ! nominal ice melting point (degC)
     &           ticemn=-50.0,  ! minimum ice surface temperature (degC)
     &           ticemx=  0.0,  ! maximum ice surface temperature (degC)
     &           salice=  5.0,  ! salinity of ice (psu)
     &           rhoice=917.0,  ! density  of ice (kg/m**3)
     &           fusion=334.e3, ! latent heat of fusion (J/kg)
     &           meltmx=  5.e-6)! max. ice melting rate (m/sec), 0.432 m/day
c
      real       fluxmx
      parameter (fluxmx=meltmx*fusion*rhoice)  !max. ice melting flux (W/m^2)
c
      include 'stmt_fns.h'
c
c --- energy loan: add extra energy to the ocean to keep SST from dropping
c --- below tfrz in winter. return this borrowed energy to the 'energy bank'
c --- in summer as quickly as surflx < 0 allows.
c
c --- salt loan: analogous to energy loan.
c
      margin = 0  ! no horizontal derivatives
c
      thicmx=0.0
!$OMP PARALLEL DO PRIVATE(j,l,i,tsur,tmxl,smxl,tfrz,borrow,paybak,f2t)
!$OMP&            REDUCTION(MAX:thicmx)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 10 j=1-margin,jj+margin
      do 11 l=1,isp(j)
      do 11 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
cdiag util1(i,j)=thkice(i,j)*scp2(i,j)
c
c --- calculate hypothetical mixed-layer temp due to diab. forcing 
      f2t =(delt1*g)/(spcifh*dpbl(i,j))
      tmxl=0.5*(temp(i,j,1,m)+temp(i,j,1,n)) + surflx(i,j)*f2t
      smxl=0.5*(saln(i,j,1,m)+saln(i,j,1,n))
      tfrz=tfrz_0 + smxl*tfrz_s  !salinity dependent freezing point
c
      borrow=(tfrz-tmxl)/f2t  ! W/m^2
c
c --- 8 time step e-folding time scale back to tfrz
c --- limit heat flux range (for both forming and melting ice)
      borrow=max( -fluxmx, min( fluxmx, 0.125*borrow ) )
c
      if (tmxl.lt.tfrz) then
c
c ---   add energy to move tmxl towards tfrz (only if tmxl < tfrz)
c
cdiag   if (thkice(i,j).le.0.0) then
cdiag     write (lp,'(i9,a,2i5)')
cdiag&      nstep,'  new ice forms at i,j=',i+i0,j+j0
cdiag   endif
        thkice(i,j)=thkice(i,j)+borrow*(delt1/(fusion*rhoice))
        surflx(i,j)=surflx(i,j)+borrow
        salflx(i,j)=salflx(i,j)+borrow*(smxl-salice)*(1.0/fusion)
c
cdiag   if (i.eq.itest .and. j.eq.jtest)then
cdiag     write (lp,'(i9,2i5,a,4f9.3)')
cdiag&      nstep,i+i0,j,+j0'  sst,tmxl,ice,loan:',
cdiag&      temp(i,j,1,n),tmxl,thkice(i,j),borrow
cdiag   endif
c
      else if (thkice(i,j).gt.0.0) then  ! tmxl > tfrz
c
c ---   ice, so return the borrowed amount whenever tmxl > tfrz
c
        paybak=min( -borrow, thkice(i,j)*(fusion*rhoice/delt1) )
        thkice(i,j)=thkice(i,j)-paybak*(delt1/(fusion*rhoice))
        surflx(i,j)=surflx(i,j)-paybak
        salflx(i,j)=salflx(i,j)-paybak*(smxl-salice)*(1.0/fusion)
c
cdiag   if (thkice(i,j).le.0.0) then
cdiag     write (lp,'(i9,a,2i5)')
cdiag&      nstep,'  ice depleted at  i,j=',i+i0,j+j0
cdiag   endif
c
cdiag   if (i.eq.itest .and. j.eq.jtest)then
cdiag     write (lp,'(i9,2i5,a,4f9.3)')
cdiag&      nstep,i+i0,j+j0,'  sst,tmxl,ice,loan:',
cdiag&      temp(i,j,1,n),tmxl,thkice(i,j),-paybak
cdiag   endif
c
      elseif (icmflg.eq.2) then  ! tmxl > tfrz & thkice(i,j) == 0.0
c
c ---   no ice, so add extra cooling under the ice mask (tsur<=tfrz_n)
c ---   don't allow a new tsur maximum, to preserve sea ice
c
        if     (yrflag.lt.2) then
          tsur = min( max( surtmp(i,j,l0), surtmp(i,j,l1),
     &                     surtmp(i,j,l2), surtmp(i,j,l3) ),
     &                surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1+
     &                surtmp(i,j,l2)*w2+surtmp(i,j,l3)*w3   )
        else
          tsur = min( max( surtmp(i,j,l0), surtmp(i,j,l1) ),
     &                surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1   )
        endif
        if     (tsur.le.tfrz_n) then
          surflx(i,j)=surflx(i,j)+borrow
        endif
      end if
c
      icex(i,j)=max(thkice(i,j)-hicemx,0.0)  ! ice exceeding hicemx
      thicmx=max(thicmx,thkice(i,j))
 11   continue
 10   continue
!$OMP END PARALLEL DO
c
      call xcmaxr(thicmx)
c
cdiag call xcsum(d1, util1,ip)
cdiag if (d1.gt.0.0) then
cdiag   if (mnproc.eq.1) then
cdiag   write(lp,'(a,2f11.4)') 'total ice:',d1/area
cdiag   call flush(lp)
cdiag   endif
cdiag endif
c
c --- spread out portion of ice thicker than hicemx
      if (thicmx.gt.hicemx) then
        call psmooth(icex, margin)
      endif
c
!$OMP PARALLEL DO PRIVATE(j,l,i,hice,smxl,tfrz)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            thkice(i,j)=icex(i,j)+min(thkice(i,j),hicemx)
c
c ---       compute fractional ice coverage for energy flux calculation
            if (thkice(i,j).lt.1.e-5*hicemn) then
              covice(i,j)=0.0
            else
              covice(i,j)=min(1.0,thkice(i,j)*(1.0/hicemn))
              hice=thkice(i,j)/covice(i,j)  !minimum of hicemn
            end if
c
c ---       compute ice surface temperature
            if     (covice(i,j).eq.0.0) then
              temice(i,j)=ticemx
            elseif (ticegr.eq.0.0) then  !use surtmp
              temice(i,j)=max( ticemn,
     &                         min( ticemx,
     &                              surtmp(i,j,l0)*w0+
     &                              surtmp(i,j,l1)*w1+
     &                              surtmp(i,j,l2)*w2+
     &                              surtmp(i,j,l3)*w3  ) )
            else
              temice(i,j)=max( ticemn, ticemx-ticegr*hice )
            endif
          enddo !i
        enddo !l
      enddo !j
!$OMP END PARALLEL DO
c
cdiag thicmx=0.0
cdiag do j=1,jj
cdiag   do l=1,isp(j)
cdiag     do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
cdiag       thicmx=max(thicmx,thkice(i,j))
cdiag     enddo
cdiag   enddo
cdiag enddo
cdiag call xcmaxr(thicmx)
cdiag if (thicmx.gt.0.0) then
cdiag   do j=1,jj
cdiag     do l=1,isp(j)
cdiag       do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
cdiag         if (thkice(i,j).eq.thicmx) then
cdiag           write(lp,'(a,2i5,f11.3)')
cdiag&            'max. ice thknss at',i+i0,j+j0,thicmx
cdiag         endif
cdiag       enddo
cdiag     enddo
cdiag   enddo
cdiag   call xcsync(flush_lp)
cdiag endif
c
      return
      end subroutine icloan
c
c
c> Revision history
c>
c> June 2000 - conversion to SI units
c> July 2000 - switched sign convention for vertical fluxes (now >0 if down)
c> May  2003 - added option to impose an ice mask
c> June 2003 - added 8 time step e-folding time scale
c> June 2003 - limited rate of ice formation
c> June 2003 - replaced constant saldif with smxl-salice
c> Mar. 2005 - freezing point linearly dependent on salinity
c> Mar. 2005 - ice surface temperature optionally from surtmp
