      subroutine thermf(m,n)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- ---------------
c --- thermal forcing
c --- ---------------
c
      integer i,j,k,nm,l
      real    pwl
      real*8  t1mean,s1mean,tmean,smean,pmean,rmean,
     &        rareac,runsec,secpyr
      real*8  d1,d2,d3,d4
c
      real*8  tmean0,smean0,rmean0
      save    tmean0,smean0,rmean0
c
      include 'stmt_fns.h'
c
      margin = 0  ! no horizontal derivatives
c
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do k=1,kk
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
            enddo
          enddo
        enddo
      enddo
c
c --- ----------------------------
c --- thermal forcing at sidewalls
c --- ----------------------------
c
      if (relax .and. nstep.gt.2) then
c
      do 53 j=1-margin,jj+margin
      do 53 l=1,isp(j)
      do 53 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
        if (rmu(i,j).ne.0.0) then
          k=1
          saln(i,j,k,n)=saln(i,j,k,n)+delt1*rmu(i,j)*
     &       (swall(i,j,k,lr0)*wr0+swall(i,j,k,lr1)*wr1
     &       +swall(i,j,k,lr2)*wr2+swall(i,j,k,lr3)*wr3
     &       - saln(i,j,k,n))
          temp(i,j,k,n)=temp(i,j,k,n)+delt1*rmu(i,j)*
     &       (twall(i,j,k,lr0)*wr0+twall(i,j,k,lr1)*wr1
     &       +twall(i,j,k,lr2)*wr2+twall(i,j,k,lr3)*wr3
     &       - temp(i,j,k,n))
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
c
          if     (hybrid) then
            do k=kk,2,-1
              pwl=pwall(i,j,k,lr0)*wr0+pwall(i,j,k,lr1)*wr1
     &           +pwall(i,j,k,lr2)*wr2+pwall(i,j,k,lr3)*wr3
              if     (pwl.gt.p(i,j,kk+1)-tencm) then
                pwl=p(i,j,kk+1)
              endif
              p(i,j,k)=min(p(i,j,k+1),
     &                     p(i,j,k)+delt1*rmu(i,j)*(pwl-p(i,j,k)))
              dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
c
              saln(i,j,k,n)=saln(i,j,k,n)+delt1*rmu(i,j)*
     &           (swall(i,j,k,lr0)*wr0+swall(i,j,k,lr1)*wr1
     &           +swall(i,j,k,lr2)*wr2+swall(i,j,k,lr3)*wr3
     &           - saln(i,j,k,n))
              if     (k.le.nhybrd .and. pwl.lt.p(i,j,kk+1)) then
                temp(i,j,k,n)=temp(i,j,k,n)+delt1*rmu(i,j)*
     &             (twall(i,j,k,lr0)*wr0+twall(i,j,k,lr1)*wr1
     &             +twall(i,j,k,lr2)*wr2+twall(i,j,k,lr3)*wr3
     &             - temp(i,j,k,n))
                th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
              else
                th3d(i,j,k,n)=       theta(k)
                temp(i,j,k,n)=tofsig(theta(k)+thbase,saln(i,j,k,n))
              endif
            enddo
          else  ! isopyc
            do k=kk,2,-1
              saln(i,j,k,n)=saln(i,j,k,n)+delt1*rmu(i,j)*
     &           (swall(i,j,k,lr0)*wr0+swall(i,j,k,lr1)*wr1
     &           +swall(i,j,k,lr2)*wr2+swall(i,j,k,lr3)*wr3
     &           - saln(i,j,k,n))
              temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
              if (k.ge.3) then
                pwl=pwall(i,j,k,lr0)*wr0+pwall(i,j,k,lr1)*wr1
     &             +pwall(i,j,k,lr2)*wr2+pwall(i,j,k,lr3)*wr3
                pwl=max(p(i,j,2),pwl)
                if     (pwl.gt.p(i,j,kk+1)-tencm) then
                  pwl=p(i,j,kk+1)
                endif
                p(i,j,k)=min(p(i,j,k+1),
     &                       p(i,j,k)+delt1*rmu(i,j)*(pwl-p(i,j,k)))
              endif
              dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
            enddo
          endif
        endif
 53   continue
c
      endif              		!  relax = .true.
c
c --- --------------------------------
c --- thermal forcing of ocean surface
c --- --------------------------------
c
      if (thermo) then
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call thermfj(m,n, j)
      enddo
      call xcsum(d1, util1,ip)
      call xcsum(d2, util2,ip)
      call xcsum(d3, util3,ip)
      call xcsum(d4, util4,ip)
      watcum=watcum+d1
      empcum=empcum+d2
      t1mean=d3
      s1mean=d4
c
      if (nstep.eq.nstep1+1 .or. diagno) then
!$OMP   PARALLEL DO PRIVATE(j,k,nm,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          k=1
            nm=n
            if (nstep.eq.nstep1+1) nm=m
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=               dp(i,j,k,nm)*scp2(i,j)
                util2(i,j)=temp(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
                util3(i,j)=saln(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
                util4(i,j)=th3d(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
              enddo
            enddo
          do k=2,kk
            nm=n
            if (nstep.eq.nstep1+1) nm=m
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=util1(i,j)+    dp(i,j,k,nm)*scp2(i,j)
                util2(i,j)=util2(i,j)+
     &                     temp(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
                util3(i,j)=util3(i,j)+
     &                     saln(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
                util4(i,j)=util4(i,j)+
     &                     th3d(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
              enddo
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xcsum(d1, util1,ip)
        call xcsum(d2, util2,ip)
        call xcsum(d3, util3,ip)
        call xcsum(d4, util4,ip)
        pmean=d1
        tmean=d2/pmean
        smean=d3/pmean
        rmean=d4/pmean
        if     (mnproc.eq.1) then
        write (lp,'(i9,a,3f9.3)') 
     &    nstep,' mean basin temp, saln, dens ',
     &    tmean,smean,rmean+thbase
        endif !1st tile
        if     (nstep.eq.nstep1+1) then
c
c ---     save initial basin means.
          tmean0=tmean
          smean0=smean
          rmean0=rmean
        else
c
c ---     diagnostic printout of fluxes.
          rareac=1.0/(area*(nstep-nstep1))
          runsec=   baclin*(nstep-nstep1)
          if      (yrflag.eq.0) then
            secpyr=360.00d0*86400.0d0
          elseif (yrflag.lt.3) then
            secpyr=366.00d0*86400.0d0
          elseif (yrflag.ge.3) then
            secpyr=365.25d0*86400.0d0
          endif
          if     (mnproc.eq.1) then
          write (lp,'(i9,a,2f9.3)') 
     &     nstep,' mean surface temp and saln  ',
     &     t1mean/area,s1mean/area
          write (lp,'(i9,a,2f9.3,a)') 
     &     nstep,' energy residual (atmos,tot) ',
     &     watcum*rareac,
     &     (tmean-tmean0)*(spcifh*avgbot/thref)/runsec,
     &    ' (w/m**2)'
c ---     note that empcum is now salflx cum.
          write (lp,'(i9,a,2f9.3,a)')
     &     nstep,'  e - p residual (atmos,tot) ',
     &     empcum*(thref/saln0)*rareac*100.0*secpyr,
     &     (smean-smean0)/(saln0*runsec)*avgbot*100.0*secpyr,
     &    ' (cm/year)'
          write (lp,'(i9,a,2f9.3)') 
     &     nstep,' temp drift per century      ',
     &     (watcum*rareac/(spcifh*avgbot/thref))*(secpyr*100.0d0),
     &     (tmean-tmean0)*(secpyr*100.0d0)/runsec
          write (lp,'(i9,a,2f9.3)') 
     &     nstep,' saln drift per century      ',
     &     (empcum*rareac/(       avgbot/thref))*(secpyr*100.0d0),
     &     (smean-smean0)*(secpyr*100.0d0)/runsec
          write (lp,'(i9,a,9x,f9.3)') 
     &     nstep,' dens drift per century      ',
     &     (rmean-rmean0)*(secpyr*100.0d0)/runsec
          endif !1st tile
          call xcsync(flush_lp)
        endif
      endif
c
      endif				!  thermo = .true.
c
      return
      end subroutine thermf
c
      subroutine thermfj(m,n, j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, j
c
      real    radfl,swfl,wind,airt,vpmx,prcp,evap,
     &        emnp,snsibl,dsgdt,tmn,smn,rmus,rmut
      real    cd0,clh,cl0,cl1,csh,
     &        rair,slat,ssen,tdif,wsph
      integer i,l
c
c --- 'ustrmn' = minimum ustar
c --- 'cekman' = constant for calculating thickness of ekman layer
c --- 'csubp'  = specific heat of air at constant pressure (j/kg/deg)
c --- 'evaplh' = latent heat of evaporation (j/kg)
C --- 'csice'  = ice-air sensible exchange coefficient
c
      real       ustrmn,cekman,csubp,evaplh,csice
      parameter (ustrmn=1.0e-5, cekman=0.7,
     &           csubp =1005.7, evaplh=2.47e6,
     &           csice=0.0006)
c
c --- 'airdns' = air density at sea level (kg/m**3)
c --- 'cd'     = drag coefficient
c --- 'ctl'    = thermal transfer coefficient (latent)
c --- 'cts1'   = thermal transfer coefficient (sensible, stable)
c --- 'cts2'   = thermal transfer coefficient (sensible, unstable)
c
      real       airdns,cd,ctl,cts1,cts2
      parameter (airdns=1.2)
      parameter (cd  =0.0013, ctl =0.0012,
     &           cts1=0.0012, cts2=0.0012)
c
c --- 'pairc'  = air pressure (mb) * 100
c --- 'rgas'   = gas constant (j/kg/k)
c --- 'tzero'  = celsius to kelvin temperature offset
c --- 'clmin'  = minimum allowed cl
c --- 'clmax'  = maximum allowed cl
c --- 'wsmin'  = minimum allowed wind speed (for cl and cd)
c --- 'wsmax'  = maximum allowed wind speed (for cl and cd)
c
      real       pairc,rgas,tzero,clmin,clmax,wsmin,wsmax
      parameter (pairc=1013.0*100.0,
     &           rgas =287.1,   tzero=273.16, 
     &           clmin=0.0003,  clmax=0.002,
     &           wsmin=3.5,     wsmax=27.5)
c
      real qsatur
      include 'stmt_fns.h'
c
c --- saturation specific humidity (lowe, j.appl.met., 16, 100-103, 1976)
      qsatur(t)=.622e-3*(6.107799961e+00+t*(4.436518521e-01
     &               +t*(1.428945805e-02+t*(2.650648471e-04
     &               +t*(3.031240396e-06+t*(2.034080948e-08
     &               +t* 6.136820929e-11))))))
c
c --- salinity relaxation coefficient
      rmus=1./(30.*86400.)
c
c --- temperature relaxation coefficient
      rmut=1./(30.*86400.)
c
c --- --------------------------------
c --- thermal forcing of ocean surface
c --- --------------------------------
c
      do 851 l=1,isp(j)
c
      do 851 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c --- wind = wind speed (m/s)
      wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
     &    +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3
c --- radfl= net       radiative thermal flux (w/m**2) into ocean
      radfl=radflx(i,j,l0)*w0+radflx(i,j,l1)*w1
     &     +radflx(i,j,l2)*w2+radflx(i,j,l3)*w3
c --- swfl = shortwave radiative thermal flux (w/m**2) into ocean
      swfl =swflx (i,j,l0)*w0+swflx (i,j,l1)*w1
     &     +swflx (i,j,l2)*w2+swflx (i,j,l3)*w3
c --- airt = air temperature (C)
      airt=airtmp(i,j,l0)*w0+airtmp(i,j,l1)*w1
     &    +airtmp(i,j,l2)*w2+airtmp(i,j,l3)*w3
c --- vpmx = water vapor mixing ratio (kg/kg)
      vpmx=vapmix(i,j,l0)*w0+vapmix(i,j,l1)*w1
     &    +vapmix(i,j,l2)*w2+vapmix(i,j,l3)*w3
c --- prcp = precipitation (m/sec; > 0 since into ocean)
      prcp=precip(i,j,l0)*w0+precip(i,j,l1)*w1
     &    +precip(i,j,l2)*w2+precip(i,j,l3)*w3
c
      if     (flxflg.eq.1) then
c
c ---   MICOM bulk air-sea flux parameterization
        if (temp(i,j,1,n).lt.airt) then
          csh=cts1
        else
          csh=cts2
        endif
c ---   evap   = evaporation (w/m**2) into atmos.
c ---   snsibl = sensible heat flux   into atmos.
        if (.not.icegln) then
          evap  =ctl*airdns*evaplh*wind*
     &           max(0.,qsatur(temp(i,j,1,n))-vpmx)
          snsibl=csh*airdns*csubp*wind*(temp(i,j,1,n)-airt)
        else  ! ice
          evap  =(1.-covice(i,j))*
     &           ctl*airdns*evaplh*wind*
     &           max(0.,qsatur(temp(i,j,1,n))-vpmx)
          snsibl=airdns*csubp*wind*
     &           ( (1.-covice(i,j))*csh  *(temp(i,j,1,n)-airt) +
     &                 covice(i,j) *csice*(temice(i,j)  -airt)  )
        endif
c ---   surflx = thermal energy flux (w/m**2) into ocean
        surflx(i,j)=radfl - snsibl - evap
c
c ---   ustar and ekman depth
        ustar(i,j)=max(ustrmn,sqrt(thref*cd*airdns)*wind)
        hekman(i,j)=ustar(i,j)*(cekman*4.0)/
     &             (abs(corio(i,j  ))+abs(corio(i+1,j  ))+
     &              abs(corio(i,j+1))+abs(corio(i+1,j+1)))
      else
c
c ---   Kara,A.B. P.A.Rochford H.E.Hurlburt 2000: Efficient and accurate
c ---   bulk parameterizations of air-sea fluxes for use in general
c ---   circulation models. J. Atmos. Ocean Tech. 17, 1421-1438.
c
        tdif = temp(i,j,1,n) - airt
        wsph = min( wsmax, max( wsmin, wind ) )
        rair = pairc / (rgas * ( tzero + airt ))
        slat = evaplh*rair
        ssen = csubp *rair
        cl0  =  0.885e-3 + 0.0748e-3 * wsph - 0.00143e-3 * wsph**2
        cl1  = -0.113e-4 + 4.89e-4   / wsph
        clh  = min( clmax, max( clmin, cl0 + cl1 * tdif ) )
        csh  = 0.9554*clh
c
c ---   evap   = evaporation         (w/m**2) into atmos.
c ---   snsibl = sensible heat flux  (w/m**2) into atmos.
c ---   surflx = thermal energy flux (w/m**2) into ocean
        if (.not.icegln) then
          evap   = slat*clh*wind*(qsatur(temp(i,j,1,n))-vpmx)
          snsibl = ssen*csh*wind* tdif
        else  ! ice
          evap   = (1.-covice(i,j))*
     &             slat*clh*wind*(qsatur(temp(i,j,1,n))-vpmx)
          snsibl = ssen*wind*
     &             ( (1.-covice(i,j))*csh  * tdif +
     &                   covice(i,j) *csice*(temice(i,j)-airt) )
        endif
        surflx(i,j) = radfl - snsibl - evap
c
c ---   ustar and ekman depth
        cd0 = 0.862e-3 + 0.088e-3 * wsph - 0.00089e-3 * wsph**2
        ustar( i,j)=max(ustrmn,sqrt(thref*cd0*rair)*wind)
        hekman(i,j)=ustar(i,j)*(cekman*4.0)/
     &             (abs(corio(i,j  ))+abs(corio(i+1,j  ))+
     &              abs(corio(i,j+1))+abs(corio(i+1,j+1)))
c
cdiag   if     (i.eq.itest.and.j.eq.jtest) then
cdiag     write(lp,'(i9,2i5,a,4f8.5)')
cdiag.    nstep,i0+i,j0+j,' cl0,cl,cs,cd    = ',cl0,clh,csh,cd0
cdiag     write(lp,'(i9,2i5,a,2f8.2,f8.5)')
cdiag.    nstep,i0+i,j0+j,' wsph,tdif,ustar = ',wsph,tdif,ustar(i,j)
cdiag     call flush(lp)
cdiag   endif
      endif
c
c --- relax to surface temperature
      if     (trelax) then
c ---   use a reference relaxation thickness (min. mixed layer depth)
c ---   actual e-folding time is (dp(i,j,1,n)/(thkmin*onem))/rmut
        surflx(i,j)=surflx(i,j)+
     &   (rmut*spcifh*thkmin*onem/g)*
     &   ( twall(i,j,1,lr0)*wr0+twall(i,j,1,lr1)*wr1
     &    +twall(i,j,1,lr2)*wr2+twall(i,j,1,lr3)*wr3
     &    - temp(i,j,1,n) )
      endif
c --- sswflx = shortwave radiative energy flux (w/m**2) into ocean
      if (.not.icegln) then
        sswflx(i,j)=swfl
      else
        sswflx(i,j)=swfl*(1.-covice(i,j))
      endif
c --- emnp = evaporation minus precipitation (m/sec) into atmos.
      if     (pcipf) then
        emnp=evap*thref/evaplh - prcp
      else
        emnp=0.0
      endif
c
c --- relax to surface salinity
      if     (srelax) then
c ---   use a reference relaxation thickness (min. mixed layer depth)
c ---   actual e-folding time is (dp(i,j,1,n)/(thkmin*onem))/rmus
        emnp=emnp+
     &   (thref/saln0)*
     &   (rmus*thkmin*onem/g)*
     &   ( swall(i,j,1,lr0)*wr0+swall(i,j,1,lr1)*wr1
     &    +swall(i,j,1,lr2)*wr2+swall(i,j,1,lr3)*wr3
     &    - saln(i,j,1,n) )
      endif
c --- salflx = salt flux (10**-3 kg/m**2/sec) into ocean
      salflx(i,j)=saln(i,j,1,n)*emnp/thref
c
c --- n o t e : t/s changes due to surflx/salflx are computed in mxlayr routine
c
c --- buoyfl = buoyancy flux, w_prime_buoyancy_prime_bar (m**2/sec**3)
c --- note: surface density increases (column is destabilized) if buoyfl > 0
      tmn=.5*(temp(i,j,1,m)+temp(i,j,1,n))
      smn=.5*(saln(i,j,1,m)+saln(i,j,1,n))
      dsgdt=dsigdt(tmn,smn)
      buoyfl(i,j)=g*thref*
     &                   (dsigds(tmn,smn)*salflx(i,j)*thref+
     &                    dsgdt          *surflx(i,j)*thref/spcifh)
c --- buoysw = shortwave radiation buoyancy flux
      buoysw(i,j)=g*thref*dsgdt          *swfl       *thref/spcifh
c
      util1(i,j)=  surflx(i,j)*scp2(i,j)
      util2(i,j)=  salflx(i,j)*scp2(i,j)
      util3(i,j)=temp(i,j,1,n)*scp2(i,j)
      util4(i,j)=saln(i,j,1,n)*scp2(i,j)
c
cdiag if     (i.eq.itest.and.j.eq.jtest) then
cdiag   write (lp,100) 
cdiag.  nstep,i0+i,j0+j,
cdiag.  '   radfl       swfl        wind        airt         sst',
cdiag.  radfl,swfl,wind,airt,temp(i,j,1,n),
cdiag.  nstep,i0+i,j0+j,
cdiag.  '   vpmx       snsibl       evap        prcp        emnp',
cdiag.  vpmx,snsibl,evap,prcp,emnp,
cdiag.  nstep,i0+i,j0+j,
cdiag.  '  surflx      sswflx      salflx      hekman       ustar',
cdiag.  surflx(i,j),sswflx(i,j),salflx(i,j),hekman(i,j),ustar(i,j),
cdiag.  nstep,i0+i,j0+j,
cdiag.  '  buoyfl      buoysw',
cdiag.  buoyfl(i,j),buoysw(i,j)
cdiag   call flush(lp)
cdiag endif
 100  format(i9,2i5,a/18x,1p5e12.4/
     &       i9,2i5,a/18x,1p5e12.4/
     &       i9,2i5,a/18x,1p5e12.4/
     &       i9,2i5,a/18x,1p2e12.4)
 851  continue
      return
      end subroutine thermfj
c
c
c> Revision history:
c>
c> Oct. 1999 - surface flux calculations modified for kpp mixed layer model,
c>             including penetrating solar radiation based on jerlov water type
c> Apr. 2000 - conversion to SI units
c> Oct  2000 - added thermfj to simplify OpenMP logic
c> Dec  2000 - modified fluxes when ice is present
c> Dec  2000 - added Kara bulk air-sea flux parameterization (flxflg=2)
