      subroutine geopar
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
c
c --- set up model parameters related to geography
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      real      realat,rlmin,rlmax,dp0kf,sigscl,dpsmin,dpsmax,dpm,dpms
      real      hmina,hminb,hmaxa,hmaxb
      integer   i,ia,ib,j,ja,jb,k,l
      integer   idim,jdim,length
      character preambl(5)*79,cline*80
c
      include 'stmt_fns.h'
c
c --- read basin depth array
c
      if     (mnproc.eq.1) then
      write (lp,'(3a)') ' reading bathymetry file from ',
     &                  flnmdep(1:len_trim(flnmdep)),'.[ab]'
      endif
      call xcsync(flush_lp)
      open (unit=9,file=flnmdep(1:len_trim(flnmdep))//'.b',status='old')
      read (     9,'(a79)')  preambl
      read (     9,'(a)')    cline
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      close(unit=9)
      if     (mnproc.eq.1) then
      write (lp,'(/(1x,a))') preambl,cline
      endif
c
      call zaiopf(flnmdep(1:len_trim(flnmdep))//'.a','old', 9)
      call zaiord(depths,ip,.false., hmina,hmaxa, 9)
      call zaiocl(9)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        endif
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1,jj
        do i= 1,ii
          if     (depths(i,j).gt.0.5*huge) then
            depths(i,j) = 0.0
          endif
        enddo
      enddo
c
      call xctilr(depths,1,1, nbdy,nbdy, halo_ps)
c
c --- introduce logorithmic k-dependence of dp0
      dp00 =onem*dp00
      dp00x=onem*dp00x
      if     (isopyc) then
        dp0k(1)=thkmin
      else
        dp0k(1)=dp00
      endif
      dpm  = dp0k(1)/onem
      dpms = dpm
      if     (mnproc.eq.1) then
      write(lp,*)
      write(lp,130) 1,dp0k(1)/onem,dpm,dpms
      endif
 130  format('dp0k(',i2,') =',f6.1,' m',
     .          '    thkns =',f6.1,' m',
     .          '    depth =',f7.1,' m')
      call xcsync(flush_lp)
c
      dp0kf=1.0
      do k=2,kk
        dp0kf=dp0kf*dp00f
        if     (k.le.nhybrd) then
          dp0k(k)=min(dp00*dp0kf,dp00x)
        else
          dp0k(k)=0.0
        endif
        dpm  = dp0k(k)/onem
        dpms = dpms + dpm
        if     (mnproc.eq.1) then
        write(lp,130) k,dp0k(k)/onem,dpm,dpms
        endif
        if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
          write(6,*) 'geopar: dp0kf  = ',dp0kf,    mnproc
          write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
        endif
        call xcsync(flush_lp)
      enddo
c
c --- sigma dependence of dp0
      if     (isopyc) then
        dp00s =onem*9999.0
        sigscl=onem
      else
        dp00s =onem*dp00s
        sigscl=onem/nsigma
      endif
      dpsmin=huge
      dpsmax=0.0
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          dp0sig(i,j) = max(dp00s,depths(i,j)*sigscl)
          dpsmin=min(dpsmin,dp0sig(i,j))
          dpsmax=max(dpsmax,dp0sig(i,j))
        enddo
      enddo
      call xcminr(dpsmin)
      call xcmaxr(dpsmax)
      if     (mnproc.eq.1) then
      write(lp,140) dpsmin/onem,dpsmax/onem
 140  format('dp0sig min,max =',f8.1,' to',f8.1,' m' /)
      endif
c
c --- determine do-loop limits for u,v,p,q points
      call bigrid(depths, util1,util2,util3)
ccc      do 3 i=1-nbdy,ii+nbdy
ccc 3    write (lp,'('' i='',i3,'' jfp,jlp='',7(1x,2i5))') i,
ccc     . (jfp(i,l),jlp(i,l),l=1,jsp(i))
ccc      do 5 j=1-nbdy,jj+nbdy
ccc 5    write (lp,'('' j='',i3,'' ifp,ilp='',7(1x,2i5))') j,
ccc     . (ifp(j,l),ilp(j,l),l=1,isp(j))
c
ccc     call prtmsk(ip,depths,util1,idm,ii,jj,0.0,1.0,
ccc    .     'bottom depth (m)')
c
c --- define coriolis parameter and grid size
      rlmin =  huge
      rlmax = -huge
      if     (mapflg.le.1) then
c --- mercator projection
      do 56 j=1,jj
      do 56 i=1,ii
      realat=alat((j0+j-0.5)-ypivn,gridn)
      rlmin = min(rlmin,realat)
      rlmax = max(rlmax,realat)
      corio(i,j)=sin(realat)*4.*pi/86164.0  ! sidereal day
c
c --- scux,scuy,scvx,scvy: grid scale at u,v points in x,y dir. respectively
      scux(i,j)=111.2e3*(alat((j0+j+0.5)-ypivn,gridn)
     .                  -alat((j0+j-0.5)-ypivn,gridn))*radian
      scvx(i,j)=111.2e3*(alat( j0+j     -ypivn,gridn)
     .                  -alat((j0+j-1.0)-ypivn,gridn))*radian
      scuy(i,j)=scux(i,j)		!  (assuming square grid cells)
      scvy(i,j)=scvx(i,j)		!  (assuming square grid cells)
c
c --- size of grid cells (length x width) at u,v,p,q points resp.
      scu2(i,j)=scux(i,j)*scuy(i,j)
      scv2(i,j)=scvx(i,j)*scvy(i,j)
      scp2(i,j)=scvx(i,j)**2            !  (assuming square grid cells)
      scq2(i,j)=scux(i,j)**2		!  (assuming square grid cells)
c
      scuxi(i,j)=1.0/scux(i,j)
      scvyi(i,j)=1.0/scvy(i,j)
      scp2i(i,j)=1.0/scp2(i,j)
      scq2i(i,j)=1.0/scq2(i,j)
c
*     if (i.eq.1) then
*       if     (mnproc.eq.1) then
*       write (lp,101) 'j',j0+j,realat*radian,
*    .                 corio(i,j),scux(i,j),scvx(i,j)
*       endif
*     endif
 101  format(a1,' =',i4,' lat =',f8.3,' corio,scux,scvx =',3e12.5)
 56   continue
      elseif (mapflg.le.3) then
c --- uniform lat/lon or square beta-plane grid
      do j=1,jj
        realat=alatu((j0+j-0.5)-ypivn,gridn)
        rlmin = min(rlmin,realat)
        rlmax = max(rlmax,realat)
        do i=1,ii
          corio(i,j)=sin(realat)*4.*pi/86164.0  ! sidereal day

          scux(i,j)=111.2e3*gridn
          scvx(i,j)=111.2e3*gridn
          if     (mapflg.eq.2) then
            scuy(i,j)=111.2e3*grido*cos(alatu( j0+j     -ypivn,gridn))
            scvy(i,j)=111.2e3*grido*cos(alatu((j0+j-0.5)-ypivn,gridn))
          else
            scuy(i,j)=scux(i,j)   !  (square beta-plane grid)
            scvy(i,j)=scvx(i,j)   !  (square beta-plane grid)
          endif
          scu2(i,j)=scux(i,j)*scuy(i,j)
          scv2(i,j)=scvx(i,j)*scvy(i,j)
          scp2(i,j)=scvx(i,j)*scuy(i,j)
          scq2(i,j)=scux(i,j)*scvy(i,j)
c
          scuxi(i,j)=1.0/scux(i,j)
          scvyi(i,j)=1.0/scvy(i,j)
          scp2i(i,j)=1.0/scp2(i,j)
          scq2i(i,j)=1.0/scq2(i,j)
        enddo
*       if     (mnproc.eq.1) then
*       write (lp,101) 'j',j0+j,realat*radian,
*    .                   corio(1,j),scux(1,j),scvx(1,j)
*     endif
      enddo
      endif
c
      call xcminr(rlmin)
      call xcmaxr(rlmax)
      if     (mnproc.eq.1) then
      write (lp,102) rlmin*radian,rlmax*radian
      endif
 102  format(/ ' latitude range (vorticity grid) = ',2f9.3 /)
c
      call xctilr(corio, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scux,  1,1, nbdy,nbdy, halo_us)
      call xctilr(scuy,  1,1, nbdy,nbdy, halo_us)
      call xctilr(scvx,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(scvy,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(scu2,  1,1, nbdy,nbdy, halo_us)
      call xctilr(scv2,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(scp2,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(scq2,  1,1, nbdy,nbdy, halo_qs)
      call xctilr(scuxi, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scvyi, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scp2i, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(scq2i, 1,1, nbdy,nbdy, halo_qs)
c
cdiag call zebra(scp2,idm,ii,jj)
cdiag write (lp,'('' shown above: (mesh size)^2 at mass points'')')
cdiag call zebra(corio,idm,ii,jj)
cdiag write (lp,'('' shown above: coriolis parameter'')')
cdiag call prmsk1(iq,corio,util1,idm,ii,jj,0.,1.e+5,'coriolis parm.')
c
      do j=1,jj
        do i=1,ii
          util1(i,j)=depths(i,j)*scp2(i,j)
        enddo
      enddo
      call xcsum(avgbot, util1,ip)
      call xcsum(area,   scp2, ip)
      avgbot=avgbot/area
      if     (mnproc.eq.1) then
      write (lp,100) avgbot,area
 100  format(' mean basin depth (m) and area (10^6 km^2):',f9.1,
     .       -12p,f10.2)
      endif
      call xcsync(flush_lp)
c
c --- initialize some arrays
c --- set depthu,dpu,utotn,pgfx,depthv,dpv,vtotn,pgfy to zero everywhere,
c --- so that they can be used at "lateral neighbors" of u and v points.
c --- similarly for pbot,dp at neighbors of q points.
c
!$OMP PARALLEL DO PRIVATE(j,i,k)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          p(     i,j,1)=0.0
          pu(    i,j,1)=0.0
          pv(    i,j,1)=0.0
          utotn( i,j)=0.0
          vtotn( i,j)=0.0
          pgfx(  i,j)=0.0
          pgfy(  i,j)=0.0
          depthu(i,j)=0.0
          depthv(i,j)=0.0
          pbot(  i,j)=0.0
c
          ubavg( i,j,1)=huge
          ubavg( i,j,2)=huge
          ubavg( i,j,3)=huge
          vbavg( i,j,1)=huge
          vbavg( i,j,2)=huge
          vbavg( i,j,3)=huge
          utotm( i,j)=huge
          vtotm( i,j)=huge
          uflux( i,j)=huge
          vflux( i,j)=huge
          uflux1(i,j)=huge
          vflux1(i,j)=huge
          uflux2(i,j)=huge
          vflux2(i,j)=huge
          uflux3(i,j)=huge
          vflux3(i,j)=huge
          uja(   i,j)=huge
          ujb(   i,j)=huge
          via(   i,j)=huge
          vib(   i,j)=huge
          do k=1,kk
            dp( i,j,k,1)=0.0
            dp( i,j,k,2)=0.0
            dpu(i,j,k,1)=0.0
            dpu(i,j,k,2)=0.0
            dpv(i,j,k,1)=0.0
            dpv(i,j,k,2)=0.0
c
            u(  i,j,k,1)=huge
            u(  i,j,k,2)=huge
            v(  i,j,k,1)=huge
            v(  i,j,k,2)=huge
c
            uflx(  i,j,k)=huge
            vflx(  i,j,k)=huge
c
            dpav(  i,j,k)=0.0
            uflxav(i,j,k)=0.0
            vflxav(i,j,k)=0.0
            tracer(i,j,k)=0.0
            diaflx(i,j,k)=0.0
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l)+1)
            ubavg(i,j,1)=0.0
            ubavg(i,j,2)=0.0
            ubavg(i,j,3)=0.0
            utotm (i,j)=0.0
            uflux (i,j)=0.0
            uflux2(i,j)=0.0
            uflux3(i,j)=0.0
            uja(i,j)=0.0
            ujb(i,j)=0.0
c
            do k=1,kk
              uflx(i,j,k)=0.0
              u(i,j,k,1)=0.0
              u(i,j,k,2)=0.0
            enddo
          enddo
        enddo
      enddo
c
      call xctilr(ubavg,    1,   3, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(utotm,    1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uflux,    1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uflux2,   1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uflux3,   1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uja,      1,   1, nbdy,nbdy, halo_us)
      call xctilr(ujb,      1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflx,     1,  kk, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(u,        1,  kk, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(u,     kk+1,2*kk, nbdy,nbdy, halo_us)  ! note scalar
c
!$OMP PARALLEL DO PRIVATE(i,l,j,k)
!$OMP&         SCHEDULE(STATIC)
      do i=1,ii
        do l=1,jsp(i)
          do j=max(1,jfp(i,l)),min(jj,jlp(i,l)+1)
            vbavg(i,j,1)=0.0
            vbavg(i,j,2)=0.0
            vbavg(i,j,3)=0.0
            vtotm (i,j)=0.0
            vflux (i,j)=0.0
            vflux2(i,j)=0.0
            vflux3(i,j)=0.0
            via(i,j)=0.0
            vib(i,j)=0.0
c
            do k=1,kk
              vflx(i,j,k)=0.0
              v(i,j,k,1)=0.0
              v(i,j,k,2)=0.0
            enddo
          enddo
        enddo
      enddo
c
      call xctilr(vbavg,    1,   3, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vtotm,    1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflux,    1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflux2,   1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflux3,   1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(via,      1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vib,      1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflx,     1,  kk, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(v,        1,  kk, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(v,     kk+1,2*kk, nbdy,nbdy, halo_vs)  ! note scalar
c
      return
      end
c
c
c> Revision history:
c>
c> May  1997 - extended list of variables set to 'huge' on land
c> Oct. 1999 - added code that defines the vertical distribution of dp0
c>             used in hybgen
c> Jan. 2000 - added mapflg logic for different projections
c> Feb. 2000 - added dp00f for logorithmic z-level spacing
c> Mar. 2000 - added dp00s for sigma-spacing in shallow water
c> May  2000 - conversion to SI units (still wrong corio)
c> Feb. 2001 - removed rotated grid option
