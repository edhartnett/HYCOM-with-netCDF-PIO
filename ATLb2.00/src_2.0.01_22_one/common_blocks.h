c-----------------------------------------------------------------------------
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2) ::
     & u,v,           ! velocity components
     & dp,dpu,dpv,    ! layer thickness
     & temp,          ! temperature
     & saln,          ! salinity
     & th3d           ! potential density

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
     & p,pu,pv        ! interface pressure

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     & dpold,         ! layer thickness
     & dpoldm,        ! layer thickness
     & thstar,        ! virtual potential density
     & diaflx,        ! time integral of diapyc.flux
     & tracer         ! inert tracer (optional)

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & corio,         ! coriolis parameter
     & psikk,         ! montg.pot. in bottom layer
     & potvor,        ! potential vorticity
     & thkk           ! virtual potential density in bottom layer

      common/hycom1r/ u,v,dp,dpold,dpoldm,dpu,dpv,p,pu,pv,
     &                corio,psikk,thkk,potvor,
     .                temp,saln,th3d,thstar,diaflx,tracer
      save  /hycom1r/
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     & montg,         ! montgomery potential
     & uflx,vflx,     ! mass fluxes
     & uflxav,vflxav, ! average fluxes
     & dpav           ! average fluxes

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,3) ::
     & ubavg,vbavg,   ! barotropic velocity
     & pbavg          ! barotropic pressure

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & defor1,defor2, ! deformation components
     & ubrhs, vbrhs,  ! rhs of barotropic u,v eqns.
     & utotm, vtotm,  ! total (barotrop.+baroclin.)..
     & utotn, vtotn,  ! ..velocities at 2 time levels
     & uflux, vflux,  ! horizontal mass fluxes
     & uflux1,vflux1, ! more mass fluxes
     & uflux2,vflux2, ! more mass fluxes
     & uflux3,vflux3  ! more mass fluxes

      common/hycom2r/ montg,defor1,defor2,ubavg,vbavg,pbavg,
     &                ubrhs,vbrhs,utotm,vtotm,utotn,vtotn,
     &                uflux,vflux,uflux1,vflux1,uflux2,vflux2,
     &                uflux3,vflux3,uflx,vflx,uflxav,vflxav,dpav
      save  /hycom2r/
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & util1,util2,   ! arrays for temporary storage
     & util3,util4,   ! arrays for temporary storage
     & scux, scuy,    ! mesh size at u pts in x,y dir.
     & scvx, scvy,    ! mesh size at v pts in x,y dir.
     & scu2, scv2,    ! grid box size at u,v pts
     & scp2, scq2,    ! grid box size at p,q pts
     & scuxi,scvyi,   ! inverses of scux,scvy
     & scp2i,scq2i,   ! inverses of scpx,scqy
     & pgfx, pgfy,    ! horiz. presssure gradient
     & gradx,grady,   ! horiz. presssure gradient
     & depthu,depthv, ! bottom pres. at u,v points
     & pvtrop,        ! pot.vort. of barotropic flow
     & depths,        ! water depth
     & drag           ! bottom drag

      common/hycom3r/ util1,util2,util3,util4,
     &                scux,scuy,scvx,scvy,scuxi,scvyi,
     &                scu2,scv2,scp2,scq2,scp2i,scq2i,
     &                pgfx,pgfy,gradx,grady,
     &                depthu,depthv,pvtrop,depths,drag
      save  /hycom3r/
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & uja,   ujb,    ! velocities at lateral ..
     & via,   vib,    !       .. neighbor points
     & pbot,          ! bottom pressure at t=0
     & sgain,         ! salin.changes from diapyc.mix.
     & surflx,        ! surface net thermal energy flux
     & sswflx,        ! surface swv thermal energy flux
     & salflx,        ! surface salinity flux
     & buoyfl,        ! net surface buoyancy flux
     & buoysw,        ! shortwave   buoyancy flux
     & ustar,         ! friction velocity
     & turgen,        ! turb.kin.energ. generation
     & thkice,        ! grid-cell avg. ice thknss (m)
     & covice,        ! ice coverage (rel.units)
     & temice         ! ice surf.temp.

      common/hycom4r/ uja,ujb,via,vib,pbot,
     &                sgain,surflx,sswflx,salflx,buoyfl,buoysw,
     &                ustar,turgen,thkice,covice,temice
      save  /hycom4r/
c
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & klist,         ! k-index
     & jerlov         ! jerlov water type 1-5

      common/hycom4i/ klist,jerlov
      save  /hycom4i/
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     & dpmixl,        ! mixed layer depth
     & t1sav,         ! upper sublayer temperature
     & s1sav,         ! upper sublayer salinity
     & tmlb,          ! temp in lyr. containing mlb.
     & smlb           ! saln in lyr. containing mlb

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & hekman,        ! ekman layer thickness
     & dpbl,          ! turbulent boundary layer depth
     & dpmold,        ! mixed layer depth
     & tmix,          ! mixed layer temperature
     & smix,          ! mixed layer salinity
     & thmix,         ! mixed layer potential density
     & umix,  vmix,   ! mixed layer velocity
     & dp0sig         ! minimum sigma   separation

      real, dimension (kdm) ::
     & dp0k           ! minimum z-layer separation

      common/hycom5r/ hekman,dpbl,dpmold,dpmixl,tmix,smix,thmix,
     &                t1sav,s1sav,tmlb,smlb,umix,vmix,dp0sig,dp0k
      save  /hycom5r/
c
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     & nmlb           ! layer containing mlb.

      common/hycom5i/ nmlb
      save  /hycom5i/
c
c ---  s w i t c h e s    (if set to .true., then...)
c --- diagno      output model fields and diagnostic messages
c --- thermo      use thermodynamic forcing (flxflg>0)
c --- windf       use wind stress   forcing (wndflg>0)
c --- pcipf       use evap-precip surface salinity flux
c --- relax       activate lateral boundary nudging
c --- srelax      activate surface salinity nudging
c --- trelax      activate surface temperature nudging
c --- relaxf      input relaxation fields (relax.or.srelax.or.trelax)
c --- hybrid      use hybrid vertical coordinates
c --- isopyc      use isopycnic vertical coordinates (MICOM mode)
c --- tbaric      include thermobaricity (kappaf)
c --- icegln      use energy loan ice model (iceflg==1)
c --- mxlkta      KT:  activate    original mixed layer model (mlflag==2)
c --- mxlktb      KT:  activate alternative mixed layer model (mlflag==3)
c --- mxlkrt      KT:  activate MICOM or HYCOM Kraus-Turner (mlflag==2,3)
c --- pensol      KT:  activate penetrating solar radiation
c --- mxlkpp      KPP: activate mixed layer model (mlflag==1)
c --- shinst      KPP: activate shear instability mixing
c --- dbdiff      KPP: activate double diffusion  mixing
c --- nonloc      KPP: activate nonlocal b. layer mixing
c --- difsmo      KPP: activate horiz smooth diff coeffs
c --- trcrin      initialize tracer from restart file
c --- trcout      advect tracer and save results in history/restart file
c
      logical       diagno,thermo,windf,pcipf,
     .              relax,srelax,trelax,relaxf,
     .              hybrid,isopyc,tbaric,icegln,
     .              mxlkta,mxlktb,mxlkrt,pensol,
     .              mxlkpp,shinst,dbdiff,nonloc,difsmo,
     .              trcrin,trcout
      common/swtchs/diagno,thermo,windf,pcipf,
     .              relax,srelax,trelax,relaxf,
     .              hybrid,isopyc,tbaric,icegln,
     .              mxlkta,mxlktb,mxlkrt,pensol,
     .              mxlkpp,shinst,dbdiff,nonloc,difsmo,
     .              trcrin,trcout
      save  /swtchs/
c
c ---  t e x t
c ---  ctitle     four lines describing the simulation
c
      character*80    ctitle
      common/hycom1c/ ctitle(4)
      save  /hycom1c/
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,4) ::
     & pwall,         ! pressure b.c. at sidewalls
     & swall,         ! salinity b.c. at sidewalls
     & twall          ! temp.    b.c. at sidewalls
 
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) ::
     & taux,          ! wind stress in x direction
     & tauy,          ! wind stress in y direction
     & wndspd,        ! wind speed (tke source)
     & airtmp,        ! air temperature
     & vapmix,        ! atmosph. vapor mixing ratio
     & precip,        ! precipitation
     & radflx,        ! net solar radiation
     & swflx          ! net shortwave radiation

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & rmu            ! weights for s.w.b.c. relax

      real, dimension (5) ::
     & betard,        ! red  extinction coefficient
     & betabl,        ! blue extinction coefficient
     & redfac         ! fract. of penetr. red light

      common/frcing/ taux,tauy,wndspd,
     &               airtmp,vapmix,precip,radflx,swflx,
     &               rmu,pwall,swall,twall,betard,betabl,redfac
      save  /frcing/
c
c --- kpp variables
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
     . zgrid          !  grid levels in centimeters
     .,vcty           !  vert. viscosity coefficient
     .,difs           !  vert. scalar diffusivity
     .,dift           !  vert. temperature diffusivity
     .,ghats          !  nonlocal transport

      real ::
     . vonk           !  von karman constant
     .,zmin,zmax      !  zehat limits for table
     .,umin,umax      !  ustar limits for table
     .,epsilon        ! vertical coordinate scale factor
     .,cmonob         ! constant for calculating monin-obukov length
     .,rinfty         ! KPP: value for calculating rshear instability
     .,difm0          ! KPP: max viscosity   due to shear instability
     .,difs0          ! KPP: max diffusivity due to shear instability
     .,difmiw         ! KPP: background/internal wave viscosity   (m^2/s)
     .,difsiw         ! KPP: background/internal wave diffusivity (m^2/s)
     .,dsfmax         ! KPP: salt fingering diffusivity factor    (m^2/s)
     .,rrho0          ! KPP: salt fingering rp=(alpha*delT)/(beta*delS)
     .,ricr           ! KPP: critical bulk richardson number
     .,cs             ! KPP: value for nonlocal flux term
     .,cstar          ! KPP: value for nonlocal flux term
     .,cv             ! KPP: value for turb shear contributon to bulk rich. no.
     .,c11            ! KPP: value for turb velocity scale
     .,deltaz         ! delta zehat in table
     .,deltau         ! delta ustar in table
     .,vtc            ! constant for estimating background shear in rib calc.
     .,cg             ! constant for estimating nonlocal flux term of diff. eq.
     .,dp0enh         ! dist. for tapering diff. enhancement at interface nbl-1

      common/kppr/ zgrid,vcty,dift,difs,ghats,
     &             vonk,zmin,zmax,umin,umax,
     &             rinfty,difm0,difs0,difmiw,difsiw,
     &             rrho0,dsfmax,ricr,epsilon,cmonob,cs,cstar,cv,c11,
     &             deltaz,deltau,vtc,cg,dp0enh
      save  /kppr/
c
      integer ::
     & niter          ! KPP: iterations for semi-implicit soln. (2 recomended)
      common/kppi/ niter
      save  /kppi/
c
      real            time,delt1,dlt,
     &                w0,w1,w2,w3,wr0,wr1,wr2,wr3
      common/varblsr/ time,delt1,dlt,
     &                w0,w1,w2,w3,wr0,wr1,wr2,wr3
      save  /varblsr/
c
      real*8           area,avgbot,watcum,empcum
      common/varblsd/  area,avgbot,watcum,empcum
      save  /varblsd/
c
      integer         nstep,nstep1,nstep2,lstep,
     &                l0,l1,l2,l3,lr0,lr1,lr2,lr3
      common/varblsi/ nstep,nstep1,nstep2,lstep,
     &                l0,l1,l2,l3,lr0,lr1,lr2,lr3
      save  /varblsi/
c
c --- 'saln0'  = initial salinity value
c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c --- 'veldff' = diffusion velocity (m/s) for momentum dissipation
c --- 'temdff' = diffusion velocity (m/s) for temp/salin. mixing
c --- 'thkdff' = diffusion velocity (m/s) for thickness diffusion
c --- 'viscos' is nondimensional, used in deformation-dependent viscosity
c --- 'biharm' = fraction of diffusion that is biharmonic (0.0 to 1.0)
c --- 'vertmx' = diffusion velocity (m/s) for mom.mixing across mix.layr.base
c --- 'diapyc' = KT: diapycnal diffusivity x buoyancy freq. (m**2/s**2)
c --- 'dtrate' = KT: maximum permitted m.l. detrainment rate (m/day)
c --- 'h1'     = depth interval used in lateral weighting of hor.pres.grad.
c --- 'slip'   = +1 for free-slip, -1  for non-slip boundary conditions
c --- 'cb'     = coefficient of quadratic bottom friction
c --- 'cbar'   = rms flow speed (m/s) for linear bottom friction law
c --- 'dsurfq' = number of days between model diagnostics at the surface
c --- 'diagfq' = number of days between model diagnostics
c --- 'rstrfq' = number of days between model restart output
c --- 'wuv1/2' = weights for time smoothing of u,v field
c --- 'wts1/2' = weights for time smoothing of t,s field
c --- 'wbaro'  = weight for time smoothing of barotropic u,v,p field
c --- 'thkmin' = minimum mixed-layer thickness
c --- 'thkbot' = thickness of bottom boundary layer
c --- 'sigjmp' = minimum density jump across interfaces   (theta units)
c --- 'tmljmp' = equivalent temperature jump across the mixed layer (degC)
c --- 'salmin' = minimum salinity allowed in an isopycnic layer
c --- 'dp00s'  = sigma   spacing minimum thickness (m)
c --- 'dp00'   = z-level spacing minimum thickness (m)
c --- 'dp00x'  = z-level spacing maximum thickness (m)
c --- 'dp00f'  = z-level spacing stretching factor (1.0=const.spacing)
c --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'hybflg' = hybrid generator flag (0=T&S,1=th&S,2=th&T)
c --- 'advflg' = scalar advection flag (0=T&S,1=th&S,2=th&T)
c --- 'mixfrq' = KT: number of time steps between diapycnal mixing calcs
c --- 'ntracr' = number of time steps between tracer transport
c --- 'clmflg' = climatology frequency flag (6=bimonthly,12=monthly)
c --- 'dypflg' = KT: diapycnal mixing flag (0=none,1=KPP,2=explicit)
c --- 'iniflg' = initial state flag (0=level,1=zonal,2=clim.,3=archv)
c --- 'lbflag' = lateral barotropic bndy flag (0=none,1=port,2=input)
c --- 'mapflg' = map flag (0=mercator,2=uniform,3=beta-plane,4=input)
c --- 'yrflag' = days in year flag (0=360,1=366)
c --- 'iversn' = hycom version number x10
c --- 'iexpt'  = experiment number x10
c --- 'jerlv0' = initial jerlov water type (1 to 5)
c --- 'iceflg' = ice model flag (0=none,1=energy loan model)
c --- 'wndflg' = wind stress input flag (0=none,1=on u/v grid,2=on p grid)
c --- 'flxflg' = thermal forcing flag (0=none,1=orig,2=new-flux-calc)
c
      real           sigma,theta,thbase,saln0,baclin,batrop,
     &               veldff,temdff,thkdff,viscos,biharm,vertmx,diapyc,
     &               dtrate,h1,slip,cb,cbar,
     &               dsurfq,diagfq,rstrfq,wuv1,wuv2,wts1,wts2,wbaro,
     &               thkmin,thkbot,sigjmp,tmljmp,
     &               salmin,dp00,dp00f,dp00x,dp00s
      common/parms1r/sigma(kdm),theta(kdm),thbase,saln0,baclin,batrop,
     &               veldff,temdff,thkdff,viscos,biharm,vertmx,diapyc,
     &               dtrate,h1,slip,cb,cbar,
     &               dsurfq,diagfq,rstrfq,wuv1,wuv2,wts1,wts2,wbaro,
     &               thkmin,thkbot,sigjmp,tmljmp,
     &               salmin(kdm),dp00,dp00f,dp00x,dp00s
      save  /parms1r/
c
      integer        mixfrq,nhybrd,nsigma,hybflg,advflg,ntracr,
     &               clmflg,dypflg,iniflg,lbflag,mapflg,yrflag,
     &               iversn,iexpt,jerlv0,iceflg,wndflg,flxflg
      common/parms1i/mixfrq,nhybrd,nsigma,hybflg,advflg,ntracr,
     &               clmflg,dypflg,iniflg,lbflag,mapflg,yrflag,
     &               iversn,iexpt,jerlv0,iceflg,wndflg,flxflg
      save  /parms1i/
c
c --- 'tenm,onem,...' = pressure thickness values corresponding to 10m,1m,...
c --- 'g'      = gravity acceleration
c --- 'thref'  = reference value of specific volume (m**3/kg)
c --- 'spcifh' = specific heat of sea water (j/kg/deg)
c --- 'epsil'  = small nonzero number used to prevent division by zero
c --- 'huge'   = large number used to indicate land points
c
      real          tenm,onem,tencm,onecm,onemm,
     &              g,thref,spcifh,epsil,huge,radian,pi
      common/consts/tenm,onem,tencm,onecm,onemm,
     &              g,thref,spcifh,epsil,huge,radian,pi
      save  /consts/
c
c --- grid point where detailed diagnostics are desired:
      integer       itest,jtest,ittest,jttest
      common/testpt/itest,jtest,ittest,jttest
      save  /testpt/
c
c --- standard  m a p   p r o j e c t i o n
c
c --- conventional mercator projection
c --- (i and x pointing east, j and y pointing north)
c
c --- ypivn = the j-index of the equator
c --- grido = mesh size of lat/lon grid in degrees
c --- gridn = mesh size of actual model grid in deg. longitude
c --- (xpivo,ypivo,xpivn not used).
c
      real         xpivo,ypivo,grido,xpivn,ypivn,gridn
      common/pivot/xpivo,ypivo,grido,xpivn,ypivn,gridn
      save  /pivot/
c
c --- filenames.
      character*48  flnmdep,flnmrsi,flnmrso,flnmflx,
     .              flnmarc,flnmovr,flnmfor,flnmforw
      common/iovars/flnmdep,flnmrsi,flnmrso,flnmflx,
     .              flnmarc,flnmovr,flnmfor,flnmforw
      save  /iovars/
c
c
c> Revision history:
c>
c> Feb. 2001 - added halo and converted to f90 declarations
c-----------------------------------------------------------------------------
