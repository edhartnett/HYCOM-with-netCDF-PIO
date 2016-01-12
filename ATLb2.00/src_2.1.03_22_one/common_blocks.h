c-----------------------------------------------------------------------------
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2) ::
     & u,v,           ! velocity components
     & dp,dpu,dpv,    ! layer thickness
     & temp,          ! temperature
     & saln,          ! salinity
     & th3d           ! potential density

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2,mxtrcr) ::
     & tracer         ! inert tracers

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
     & p,pu,pv        ! interface pressure

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     & dpold,         ! layer thickness
     & dpoldm,        ! layer thickness
     & thstar,        ! virtual potential density
     & diaflx         ! time integral of diapyc.flux

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & corio,         ! coriolis parameter
     & psikk,         ! montg.pot. in bottom layer
     & potvor,        ! potential vorticity
     & thkk           ! virtual potential density in bottom layer

      common/hycom1r/ u,v,dp,dpold,dpoldm,dpu,dpv,p,pu,pv,
     &                corio,psikk,thkk,potvor,
     &                temp,saln,th3d,thstar,diaflx,tracer
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
     &                uflux, vflux, uflux1,vflux1,
     &                uflux2,vflux2,uflux3,vflux3,
     &                uflx,vflx,uflxav,vflxav,dpav
      save  /hycom2r/
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & util1,util2,   ! arrays for temporary storage
     & util3,util4,   ! arrays for temporary storage
     & util5,util6,   ! arrays for temporary storage
     & plon, plat,    ! lon,lat at p pts
     & scux, scuy,    ! mesh size at u pts in x,y dir.
     & scvx, scvy,    ! mesh size at v pts in x,y dir.
     & scpx, scpy,    ! mesh size at p pts in x,y dir.
     & scqx, scqy,    ! mesh size at q pts in x,y dir.
     & scu2, scv2,    ! grid box area at u,v pts
     & scp2, scq2,    ! grid box area at p,q pts
     & scp2i,scq2i,   ! inverses of scp2,scq2
     & scuxi,scvyi,   ! inverses of scux,scvy
     & aspux,aspuy,   ! u-grid aspect ratios for diffusion
     & aspvx,aspvy,   ! v-grid aspect ratios for diffusion
     & pgfx, pgfy,    ! horiz. presssure gradient
     & gradx,grady,   ! horiz. presssure gradient
     & depthu,depthv, ! bottom pres. at u,v points
     & pvtrop,        ! pot.vort. of barotropic flow
     & depths,        ! water depth
     & drag           ! bottom drag

      common/hycom3r/ util1,util2,util3,util4,util5,util6,
     &                plon,plat,
     &                scux,scuy,scvx,scvy,scuxi,scvyi,
     &                scpx,scpy,scqx,scqy,
     &                scu2,scv2,scp2,scq2,scp2i,scq2i,
     &                aspux,aspuy,aspvx,aspvy,
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
     & ustar,         ! friction velocity
     & turgen,        ! turb.kin.energ. generation
     & thkice,        ! grid-cell avg. ice thknss (m)
     & covice,        ! ice coverage (rel.units)
     & temice         ! ice surf.temp.

      common/hycom4r/ uja,ujb,via,vib,pbot,
     &                sgain,surflx,sswflx,salflx,
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
     & umix,  vmix    ! mixed layer velocity

      real, dimension (kdm) ::
     & dp0k,          ! minimum deep    z-layer separation
     & ds0k,          ! minimum shallow z-layer separation
     & dssk           ! sigma depth scale factor

      common/hycom5r/ hekman,dpbl,dpmold,dpmixl,tmix,smix,thmix,
     &                t1sav,s1sav,tmlb,smlb,umix,vmix,dp0k,ds0k,dssk
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
c --- priver      use river precip bogas
c --- rivera      annual-only river precip bogas
c --- relax       activate lateral boundary T/S/p  climatological nudging
c --- srelax      activate surface salinity        climatological nudging
c --- trelax      activate surface temperature     climatological nudging
c --- trcrlx      activate lateral boundary tracer climatological nudging
c --- relaxf      input T/S/p   relaxation fields
c --- relaxs      input surface relaxation fields only
c --- relaxt      input tracer  relaxation fields
c --- hybrid      use hybrid vertical coordinates
c --- isopyc      use isopycnic vertical coordinates (MICOM mode)
c --- tbaric      include thermobaricity (kappaf)
c --- icegln      use energy loan ice model (iceflg==1)
c --- mxlkta      KT:    activate    original mixed layer model (mlflag==2)
c --- mxlktb      KT:    activate alternative mixed layer model (mlflag==3)
c --- mxlkrt      KT:    activate MICOM or HYCOM Kraus-Turner (mlflag==2,3)
c --- pensol      KT:    activate penetrating solar radiation
c --- mxlkpp      KPP:   activate mixed layer model (mlflag==1)
c --- shinst      KPP:   activate shear instability mixing
c --- dbdiff      KPP:   activate double diffusion  mixing
c --- nonloc      KPP:   activate nonlocal b. layer mixing
c --- difsmo      KPP:   activate horiz smooth diff coeffs
c --- mxlmy       MY2.5: activate mixed layer model (mlflag==5)
c --- mxlpwp      PWP:   activate mixed layer model (mlflag==4)
c --- trcrin      initialize tracer from restart file
c --- trcout      advect tracer and save results in history/restart file
c
      logical       diagno,thermo,windf,pcipf,priver,rivera,
     &              relax,srelax,trelax,trcrlx,relaxf,relaxs,relaxt,
     &              hybrid,isopyc,tbaric,icegln,
     &              mxlkta,mxlktb,mxlkrt,pensol,
     &              mxlkpp,shinst,dbdiff,nonloc,difsmo,
     &              mxlmy,mxlpwp,trcrin,trcout
      common/swtchs/diagno,thermo,windf,pcipf,priver,rivera,
     &              relax,srelax,trelax,trcrlx,relaxf,relaxs,relaxt,
     &              hybrid,isopyc,tbaric,icegln,
     &              mxlkta,mxlktb,mxlkrt,pensol,
     &              mxlkpp,shinst,dbdiff,nonloc,difsmo,
     &              mxlmy,mxlpwp,trcrin,trcout
      save  /swtchs/
c
c ---  t e x t
c ---  ctitle     four lines describing the simulation
c
      character*80    ctitle
      common/hycom1c/ ctitle(4)
      save  /hycom1c/
c
c --- atmospheric forcing fields
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) ::
     & taux,          ! wind stress in x direction
     & tauy,          ! wind stress in y direction
     & wndspd,        ! wind speed (tke source)
     & airtmp,        ! air temperature
     & vapmix,        ! atmosph. vapor mixing ratio
     & precip,        ! precipitation
     & radflx,        ! net solar radiation
     & swflx,         ! net shortwave radiation
     & rivers         ! river inflow bogused to surface precipitation

      real, dimension (5) ::
     & betard,        ! red  extinction coefficient
     & betabl,        ! blue extinction coefficient
     & redfac         ! fract. of penetr. red light

      common/frcing/ taux,tauy,wndspd,
     &               airtmp,vapmix,precip,radflx,swflx,
     &               rivers,
     &               betard,betabl,redfac
      save  /frcing/
c
c --- surface and sidewall and nestwall boundary fields
c ---  (kkwall and kknest are either kdm or, if inactive, 1).
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkwall,4) ::
     & pwall,         ! pressure b.c. at sidewalls
     & swall,         ! salinity b.c. at sidewalls
     & twall          ! temp.    b.c. at sidewalls

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkwall,4,
     &                                                   mxtrcr) ::
     & trwall         ! tracer   b.c. at sidewalls

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kknest,2) ::
     & pnest,         ! pressure b.c. at nestwalls
     & snest,         ! salinity b.c. at nestwalls
     & tnest,         ! temp.    b.c. at nestwalls
     & unest,         ! u-vel.   b.c. at nestwalls
     & vnest          ! v-vel.   b.c. at nestwalls
 
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     & ubnest,        ! barotropic u-velocity at nestwalls
     & vbnest,        ! barotropic v-velocity at nestwalls
     & ubpnst,        ! barotropic u-velocity at nestwalls on p-grid
     & vbpnst,        ! barotropic v-velocity at nestwalls on p-grid
     & pbnest         ! barotropic pressure   at nestwalls

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & rmu,           ! weights for s.w.b.c. relax
     & rmutr,         ! weightc for trc.b.c. relax
     & rmun           ! weights for n.w.b.c. relax

      common/wall1r/ rmu,  pwall,swall,twall,
     &               rmutr,trwall,
     &               rmun, pnest,snest,tnest,unest,vnest,
     &               ubnest,vbnest,ubpnst,vbpnst,pbnest
      save  /wall1r/

      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & maskbc         ! mask for nested barotropic boundary condition

      common/wall1i/ maskbc
      save  /wall1i/
c
c --- pwp variables
      real ::
     & rigc           ! PWP: critical gradient richardson number
     &,ribc           ! PWP: critical bulk richardson number

      common/pwpr/ rigc,ribc
      save  /pwpr/
c
c --- m-y 2.5 variables
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1,2) ::
     & q2             !  tke
     &,q2l            !  tke * turbulent length scale

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1) ::
     & difqmy         !  tke diffusivity
     &,vctymy         !  viscosity on mellor-yamada vertical grid
     &,diftmy         !  temperature diffusivity on mellor-yamada vertical grid

      real ::
     & ghc            !  constant for calculating tke production
     &,sef            !  constant for calculating tke production
     &,smll           !  constant for calculating tke
     &,const1         !  coefficient for estimating surface and bottom bc's
     &,coef4          !  coefficient for calculating  viscosity/diffusivity
     &,coef5          !  coefficient for calculating  viscosity/diffusivity
     &,a1my           !  coefficient for calculating  viscosity/diffusivity
     &,b1my           !  coefficient for calculating  viscosity/diffusivity
     &,a2my           !  coefficient for calculating  viscosity/diffusivity
     &,b2my           !  coefficient for calculating  viscosity/diffusivity
     &,c1my           !  coefficient for calculating  viscosity/diffusivity
     &,e1my           !  coefficient for calculating  viscosity/diffusivity
     &,e2my           !  coefficient for calculating  viscosity/diffusivity
     &,e3my           !  coefficient for calculating  viscosity/diffusivity
c
      common/myr/ q2,q2l,difqmy,vctymy,diftmy,
     &            ghc,sef,smll,const1,
     &            coef4,coef5,a1my,b1my,a2my,b2my,c1my,e1my,e2my,e3my
      save  /myr/
c
c --- kpp variables
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
     & zgrid          !  grid levels in centimeters
     &,vcty           !  vert. viscosity coefficient
     &,difs           !  vert. scalar diffusivity
     &,dift           !  vert. temperature diffusivity
     &,ghats          !  nonlocal transport

      real ::
     & vonk           !  von karman constant
     &,zmin,zmax      !  zehat limits for table
     &,umin,umax      !  ustar limits for table
     &,epsilon        ! vertical coordinate scale factor
     &,cmonob         ! constant for calculating monin-obukov length
     &,rinfty         ! KPP: value for calculating rshear instability
     &,difm0          ! KPP: max viscosity   due to shear instability
     &,difs0          ! KPP: max diffusivity due to shear instability
     &,difmiw         ! KPP: background/internal wave viscosity   (m^2/s)
     &,difsiw         ! KPP: background/internal wave diffusivity (m^2/s)
     &,dsfmax         ! KPP: salt fingering diffusivity factor    (m^2/s)
     &,rrho0          ! KPP: salt fingering rp=(alpha*delT)/(beta*delS)
     &,ricr           ! KPP: critical bulk richardson number
     &,cs             ! KPP: value for nonlocal flux term
     &,cstar          ! KPP: value for nonlocal flux term
     &,cv             ! KPP: value for turb shear contributon to bulk rich. no.
     &,c11            ! KPP: value for turb velocity scale
     &,deltaz         ! delta zehat in table
     &,deltau         ! delta ustar in table
     &,vtc            ! constant for estimating background shear in rib calc.
     &,cg             ! constant for estimating nonlocal flux term of diff. eq.
     &,dp0enh         ! dist. for tapering diff. enhancement at interface nbl-1

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
      real*8           area,avgbot,watcum,empcum
      common/varblsd/  area,avgbot,watcum,empcum
      save  /varblsd/
c
      real            time,delt1,dlt,
     &                w0, w1, w2, w3,  ! wind  interp. scale factors
     &                wr0,wr1,wr2,wr3, ! river interp. scale factors
     &                wc0,wc1,wc2,wc3, ! clim. interp. scale factors
     &                wn0,wn1,         ! nest  interp. scale factors
     &                wb0,wb1          ! baro. interp. scale factors
      common/varblsr/ time,delt1,dlt,
     &                w0,w1,w2,w3,wr0,wr1,wr2,wr3,wc0,wc1,wc2,wc3,
     &                wn0,wn1,wb0,wb1
      save  /varblsr/
c
      integer         nstep,nstep1,nstep2,lstep,
     &                l0, l1, l2, l3,  ! wind  indexes
     &                lr0,lr1,lr2,lr3, ! river indexes
     &                lc0,lc1,lc2,lc3, ! clim. indexes
     &                ln0,ln1,         ! nest  indexes
     &                lb0,lb1          ! baro. indexes
      common/varblsi/ nstep,nstep1,nstep2,lstep,
     &                l0,l1,l2,l3,lr0,lr1,lr2,lr3,lc0,lc1,lc2,lc3,
     &                ln0,ln1,lb0,lb1
      save  /varblsi/
c
c --- 'saln0'  = initial salinity value
c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c --- 'visco2' = deformation-dependent Laplacian  viscosity factor
c --- 'visco4' = deformation-dependent biharmonic viscosity factor
c --- 'veldf2' = diffusion velocity (m/s) for Laplacian  momentum dissipation
c --- 'veldf4' = diffusion velocity (m/s) for biharmonic momentum dissipation
c --- 'temdf2' = diffusion velocity (m/s) for Laplacian  temp/salin. mixing
c --- 'thkdf2' = diffusion velocity (m/s) for Laplacian  thickness diffusion
c --- 'thkdf4' = diffusion velocity (m/s) for biharmonic thickness diffusion
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
c --- 'bnstfq' = number of days between baro. nesting archive input
c --- 'nestfq' = number of days between 3-d   nesting archive input
c --- 'wuv1/2' = weights for time smoothing of u,v field
c --- 'wts1/2' = weights for time smoothing of t,s field
c --- 'wbaro'  = weight for time smoothing of barotropic u,v,p field
c --- 'thkmlr' = reference mixed-layer thickness for SST/SSS relaxation
c --- 'thkmin' = KT/PWP: minimum mixed-layer thickness
c --- 'thkbot' = thickness of bottom boundary layer
c --- 'sigjmp' = minimum density jump across interfaces   (theta units)
c --- 'tmljmp' = equivalent temperature jump across the mixed layer (degC)
c --- 'salmin' = minimum salinity allowed in an isopycnic layer
c --- 'dp00'   = deep    z-level spacing minimum thickness (m)
c --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
c --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
c --- 'ds00'   = shallow z-level spacing minimum thickness (m)
c --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
c --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
c --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'hybflg' = hybrid generator flag (0=T&S,1=th&S,2=th&T)
c --- 'advflg' = scalar advection flag (0=T&S,1=th&S,2=th&T)
c --- 'mixfrq' = KT: number of time steps between diapycnal mixing calcs
c --- 'ntracr' = number of tracers (<=mxtrcr)
c --- 'trcflg' = tracer type flag (one per tracer)
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
     &               visco2,visco4,veldf2,veldf4,
     &               temdf2,thkdf2,thkdf4,vertmx,diapyc,
     &               dtrate,h1,slip,cb,cbar,
     &               dsurfq,diagfq,rstrfq,bnstfq,nestfq,
     &               wuv1,wuv2,wts1,wts2,wbaro,
     &               thkmlr,thkmin,thkbot,sigjmp,tmljmp,
     &               salmin,dp00,dp00f,dp00x,ds00,ds00f,ds00x
      common/parms1r/sigma(kdm),theta(kdm),thbase,saln0,baclin,batrop,
     &               visco2,visco4,veldf2,veldf4,
     &               temdf2,thkdf2,thkdf4,vertmx,diapyc,
     &               dtrate,h1,slip,cb,cbar,
     &               dsurfq,diagfq,rstrfq,bnstfq,nestfq,
     &               wuv1,wuv2,wts1,wts2,wbaro,
     &               thkmlr,thkmin,thkbot,sigjmp,tmljmp,
     &               salmin(kdm),dp00,dp00f,dp00x,ds00,ds00f,ds00x
      save  /parms1r/
c
      integer        mixfrq,nhybrd,nsigma,hybflg,advflg,
     &               ntracr,trcflg(mxtrcr),
     &               clmflg,dypflg,iniflg,lbflag,mapflg,yrflag,
     &               iversn,iexpt,jerlv0,iceflg,wndflg,flxflg
      common/parms1i/mixfrq,nhybrd,nsigma,hybflg,advflg,
     &               ntracr,trcflg,
     &               clmflg,dypflg,iniflg,lbflag,mapflg,yrflag,
     &               iversn,iexpt,jerlv0,iceflg,wndflg,flxflg
      save  /parms1i/
c
c --- 'tenm,onem,...' = pressure thickness values corresponding to 10m,1m,...
c --- 'qonem'  = 1/onem
c --- 'g'      = gravity acceleration
c --- 'thref'  = reference value of specific volume (m**3/kg)
c --- 'qthref' = 1/thref
c --- 'spcifh' = specific heat of sea water (j/kg/deg)
c --- 'epsil'  = small nonzero number used to prevent division by zero
c --- 'huge'   = large number used to indicate land points
c
      real          tenm,onem,tencm,onecm,onemm,qonem,
     &              g,thref,qthref,spcifh,epsil,huge,radian,pi
      common/consts/tenm,onem,tencm,onecm,onemm,qonem,
     &              g,thref,qthref,spcifh,epsil,huge,radian,pi
      save  /consts/
c
c --- grid point where detailed diagnostics are desired:
      integer       itest,jtest,ittest,jttest
      common/testpt/itest,jtest,ittest,jttest
      save  /testpt/
c
c --- filenames.
      character*48  flnmdep,flnmgrd,flnmrsi,flnmrso,flnmflx,
     &              flnmarc,flnmovr,flnmfor,flnmforw
      common/iovars/flnmdep,flnmgrd,flnmrsi,flnmrso,flnmflx,
     &              flnmarc,flnmovr,flnmfor,flnmforw
      save  /iovars/
c
c> Revision history:
c>
c> Feb. 2001 - added halo and converted to f90 declarations
c> Jan. 2002 - added curvilinear grid arrays and deleted /pivot/
c-----------------------------------------------------------------------------
