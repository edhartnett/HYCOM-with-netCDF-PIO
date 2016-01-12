      subroutine blkdat
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      real      day1,hybrlx
      integer   k,mlflag,thflag,trcflg1
      character sigfmt*26
c
      include 'stmt_fns.h'
c
c --- initialize common variables.
c
      open(unit=99,file='blkdat.input')
c
c --- 'lp' = logical unit number for printer output
      lp = 6
c
c --- 'g'      = gravitational acceleration (m/s**2)
c --- 'thref'  = reference value of specific volume (m**3/kg)
      g     =   9.806
      thref =   1.0e-3
c
c --- layer thicknesses in units of pressure:
      tenm   = 98060.0    ! g/thref *10.0
      onem   =  9806.0    ! g/thref
      tencm  =   980.6    ! g/thref * 0.1
      onecm  =    98.06   ! g/thref * 0.01
      onemm  =     9.806  ! g/thref * 0.001
c
      qonem  = 1.0/onem   ! thref/g
      qthref = 1.0/thref
c
c --- pi-related values
      pi = 4.d0*atan(1.d0)
      radian=pi/180.0
c
c --- four lines (80-characters) describing the simulation
      read( 99,'(a80/a80/a80/a80)') ctitle
      if (mnproc.eq.1) then
      write(lp,*)
      write(lp,'(a80/a80/a80/a80)') ctitle
      call flush(lp)
      endif !1st tile
c
c --- 'iversn' = hycom version number x10
c --- 'iexpt'  = experiment number x10
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(iversn,'iversn')
      call blkini(iexpt, 'iexpt ')
c
      if (iversn.lt.21 .or. iversn.gt.21) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3,a,i3 /)') 
     &    'error - iversn must be between',21,' and',21
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
c
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(itest ,'idm   ')
      call blkini(jtest ,'jdm   ')
c
      if     (itest.ne.itdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i5 /)') 
     &    'error - expected idm =',itdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if     (jtest.ne.jtdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i5 /)') 
     &    'error - expected jdm =',jtdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
c
c --- 'itest,jtest' = grid point where detailed diagnostics are desired
c ---                 itest=jtest=0 turns off all detailed diagnostics
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(ittest,'itest ')
      call blkini(jttest,'jtest ')
c
      if (ittest.gt.itdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i5 /)') 
     &    'error - maximum itest is',itdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (jttest.gt.jtdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i5 /)') 
     &    'error - maximum jtest is',jtdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
c
c --- map global ittest,jttest to local itest,jtest
      if     (ittest.gt.i0 .and. ittest.le.i0+ii .and.
     &        jttest.gt.j0 .and. jttest.le.j0+jj      ) then
        itest = ittest - i0
        jtest = jttest - j0
      else
        itest = -99
        jtest = -99
      endif
c
*     if (mnproc.eq.1) then
*     write(lp,*)
*     endif !1st tile
*     call xcsync(flush_lp)
*     do k= 1,ijpr
*       if     (mnproc.eq.k) then
*         write(lp,'(a,3i5)') 'mnproc,[ij]test =',mnproc,itest,jtest
*       endif
*       call xcsync(flush_lp)
*     enddo !k
c
c --- 'kdm   ' = number of layers
c --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'dp00'   = deep    z-level spacing minimum thickness (m)
c --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
c --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
c --- 'ds00'   = shallow z-level spacing minimum thickness (m)
c --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
c --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
c --- 'dp00i'  = deep iso-pycnal spacing minimum thickness (m)
c --- 'isotop' = shallowest depth for isopycnal layers     (m), <0 from file
c
c --- the above specifies a vertical coord. that is isopycnal or:
c ---  near surface z in    deep water, based on dp00,dp00x,dp00f
c ---  near surface z in shallow water, based on ds00,ds00x,ds00f and nsigma
c ---               sigma between them, based on ds00,ds00x,ds00f and nsigma
c
c --- d[ps].k/d[ps].1 = d[ps]00f**(k-1) unless limited by d[ps]00x.
c --- if d[ps]00f>1, d[ps]00x (> d[ps]00) is the maximum thickness.
c --- if d[ps]00f<1, d[ps]00x (< d[ps]00) is the minimum thickness.
c
c --- near the surface (i.e. shallower than isotop), layers are always fixed
c --- depth (z or sigma).  layer 1 is always fixed, so isotop=0.0 is not
c --- allowed.  if isotop<0.0 then isotop(1:idm,1:jdm) is a spacially
c --- varying array, input from the file iso.top.[ab].
c
c --- away from the surface, the minimum layer thickness is dp00i.
c --- to recover original hybrid behaviour,   set dp00i=dp00x
c --- for z-only, sigma-only or sigma-z only, set dp00i=dp00x
c
c --- for z-only set nsigma=0 (and ds00,ds00x,ds00f=dp00,dp00x,dp00f)
c --- for sigma-z (shallow-deep) use a very small ds00
c ---  (pure sigma-z also has ds00f=dp00f and ds00x=dp00x*ds00/dp00)
c --- for z-sigma (shallow-deep) use a very large dp00 (not recommended)
c --- for sigma-only set nsigma=kdm, dp00 large, and ds00 small
c
c --- for an entirely fixed vertical coordinate (no isopycnal layers), set
c --- isotop large or make all target densities (sigma(k), below) very small.
c
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(k,     'kdm   ')
      call blkini(nhybrd,'nhybrd')
      call blkini(nsigma,'nsigma')
      call blkinr(dp00,  'dp00  ','(a6," =",f10.4," m")')
      call blkinr(dp00x, 'dp00x ','(a6," =",f10.4," m")')
      call blkinr(dp00f, 'dp00f ','(a6," =",f10.4," ")')
      call blkinr(ds00,  'ds00  ','(a6," =",f10.4," m")')
      call blkinr(ds00x, 'ds00x ','(a6," =",f10.4," m")')
      call blkinr(ds00f, 'ds00f ','(a6," =",f10.4," ")')
      call blkinr(dp00i, 'dp00i ','(a6," =",f10.4," m")')
      call blkinr(isotop,'isotop','(a6," =",f10.4," m")')
c
c --- isopycnal (MICOM-like) iff nhybrd is 0
      isopyc = nhybrd .eq. 0
      hybrid = .not. isopyc
      if (hybrid .and. nsigma.le.1) then
        nsigma=1
      endif
c
      if     (k.ne.kdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3 /)') 
     &    'error - expected kdm =',kdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if     (nhybrd.gt.kdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3 /)') 
     &    'error - maximum nhybrd is kdm =',kdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if     (nsigma.gt.nhybrd) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3 /)') 
     &    'error - maximum nsigma is nhybrd =',nhybrd
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (isopyc .and. max(dp00,dp00x).ne.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - must have dp00x==dp00==0.0 for isopycnal case'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (dp00f.eq.1.0 .and. dp00.ne.dp00x) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - must have dp00x==dp00 for dp00f==1.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (dp00f.gt.1.0 .and. dp00.ge.dp00x) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - dp00x must be > dp00 for dp00f>1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (dp00f.lt.1.0 .and. dp00.le.dp00x) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - dp00x must be < dp00 for dp00f<1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (ds00.gt.dp00) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - must have ds00 <= dp00'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (ds00.le.0.0 .and. .not.isopyc) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - must have ds00>0.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (ds00f.eq.1.0 .and. ds00.ne.ds00x) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - must have ds00x==ds00 for ds00f==1.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (ds00f.gt.1.0 .and. ds00.ge.ds00x) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - ds00x must be > ds00 for ds00f>1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (ds00f.lt.1.0 .and. ds00.le.ds00x) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - ds00x must be < ds00 for ds00f<1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (isotop.eq.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - isotop cannot be 0.0 (layer 1 never isopycnal)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
c
c --- 'saln0'  = initial salinity value (psu), only used for iniflg<2
c --- 'locsig' = locally-referenced potential density for stability (0=F,1=T)
c --- 'kapref' = thermobaric reference state (-1=input,0=none,1,2,3=constant)
c ---              1=Arctic/Antarctic; 2=Atlantic; 3=Mediterranean
c --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
c ---            this is a check on the compile-time stmt_funcs.h setup.
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinr(saln0, 'saln0 ','(a6," =",f10.4," psu")')
      call blkinl(locsig,'locsig')
      call blkini(kapref,'kapref')
      call blkini(thflag,'thflag')
c --- kapnum is number of thermobaric reference states (1 or 2)
      if     (kapref.ne.-1) then
        kapnum=1
      else
        kapnum=2
      endif
c
      if     (kapref.ne.0 .and. thflag.ne.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')
     &    'error - kapref!=0 only allowed for sigma2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !thermobaric error
      if     (kapref.lt.-1 .or. kapref.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i1 /)')
     &    'error - kapref must be between -1 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !kapref error
      if     (thflag.eq.0) then
        if     (abs(sig(0.0,0.0)+0.136471).gt.0.001) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     &      'error - thflag not consistent with sig()'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
      elseif (thflag.eq.2) then
        if     (abs(sig(0.0,0.0)-9.77093 ).gt.0.001) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     &      'error - thflag not consistent with sig()'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
      else
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - thflag must be 0 or 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !thflag
      if     (thflag.eq.0) then
        sigfmt = '(a6," =",f10.4," sigma-0")'
      elseif (thflag.eq.2) then
        sigfmt = '(a6," =",f10.4," sigma-2")'
      endif
c
c --- 'thbase' = reference density (sigma units)
      call blkinr(thbase,'thbase',sigfmt)
c
c --- 'vsigma' = spacially varying isopycnal layer target densities (0=F,1=T)
c ---            if true, target densities input from file iso.sigma.[ab]
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinl(vsigma,'vsigma')
c
c --- 'sigma ' = isopycnal layer target densities (sigma units)
      do k=1,kdm
        call blkinr(sigma(k),'sigma ',sigfmt)
c
        if     (k.gt.1) then
          if      (sigma(k).le.sigma(k-1)) then
            if (mnproc.eq.1) then
            write(lp,'(/ a,i3 /)') 
     &        'error - sigma is not stabally stratified'
            call flush(lp)
            endif !1st tile
            call xcstop('(blkdat)')
                   stop '(blkdat)'
          endif !sigma(k).le.sigma(k-1)
        endif !k>1
      enddo !k
c
c --- 'iniflg' = initial state flag (0=level,1=zonal,2=climatology)
c --- 'jerlv0' = initial jerlov water type (1 to 5; 0 to use kpar)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(iniflg,'iniflg')
      call blkini(jerlv0,'jerlv0')
c
      if (iniflg.lt.0 .or. iniflg.gt.3) then  !inicon==3 ok, for old .inputs
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - iniflg must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (jerlv0.lt.0 .or. jerlv0.gt.5) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - jerlv0 must be 0 or between 1 and 5'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
c --- red and blue light extinction coefficients (1/pressure units)
c --- for jerlov water types 1 to 5 - fraction of penetrating red light
      betard(1) = 1.0/( 0.35*onem)
      betard(2) = 1.0/( 0.6 *onem)
      betard(3) = 1.0/( 1.0 *onem)
      betard(4) = 1.0/( 1.5 *onem)
      betard(5) = 1.0/( 1.4 *onem)
      betabl(1) = 1.0/(23.0 *onem)
      betabl(2) = 1.0/(20.0 *onem)
      betabl(3) = 1.0/(17.0 *onem)
      betabl(4) = 1.0/(14.0 *onem)
      betabl(5) = 1.0/( 7.9 *onem)
      redfac(1) = 0.58
      redfac(2) = 0.62
      redfac(3) = 0.67
      redfac(4) = 0.77
      redfac(5) = 0.78
c
c --- 'yrflag' = days in year flag (0=360,1=366,2=366Jan1,3=actual)
c --- 'dsurfq' = number of days between model diagnostics at the surface
c --- 'diagfq' = number of days between model diagnostics
c --- 'meanfq' = number of days between model diagnostics (time averaged)
c ---             NOT YET IMPLEMENTED - must be 0.0
c --- 'rstrfq' = number of days between model restart output
c --- 'bnstfq' = number of days between baro nesting archive input
c ---             (0.0 if lbflag is not 2)
c --- 'nestfq' = number of days between 3-d  nesting archive input
c ---             (0.0 turns off relaxation to 3-d nesting input)
c --- 'baclin' = baroclinic time step (seconds), int. divisor of 86400
c --- 'batrop' = barotropic time step (seconds), int. divisor of baclin/2
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(yrflag,'yrflag')
      call blkinr(dsurfq,'dsurfq','(a6," =",f10.4," days")')
      call blkinr(diagfq,'diagfq','(a6," =",f10.4," days")')
      call blkinr(meanfq,'meanfq','(a6," =",f10.4," days")')
      call blkinr(rstrfq,'rstrfq','(a6," =",f10.4," days")')
      call blkinr(bnstfq,'bnstfq','(a6," =",f10.4," days")')
      call blkinr(nestfq,'nestfq','(a6," =",f10.4," days")')
      call blkinr(baclin,'baclin','(a6," =",f10.4," sec")')
      call blkinr(batrop,'batrop','(a6," =",f10.4," sec")')
c
      if (yrflag.lt.0 .or. yrflag.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - yrflag must be between 0 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (yrflag.le.1) then  ! monthly forcing
        if     (abs(nint(86400.0/baclin)-86400.0/baclin).gt.0.01) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')
     &      'error - baclin not an integer divisor of 86400 secs'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      else  ! high frequency forcing
        if     (abs(nint(21600.0/baclin)-21600.0/baclin).gt.0.01) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')
     &      'error - baclin not an integer divisor of 21600 secs'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      endif
      if     (abs(nint(0.5*baclin/batrop)-
     &                 0.5*baclin/batrop  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')
     &    'error - batrop not an integer divisor of baclin/2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (abs(nint((dsurfq*86400.d0)/baclin)-
     &                 (dsurfq*86400.d0)/baclin  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)')
     &    'error - ',
     &    'dsurfq is not a whole number of baroclinic time steps'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      else  ! make it exact
*       write(lp,*) 'old dsurfq = ',dsurfq
        dsurfq = nint((dsurfq*86400.d0)/baclin)*(baclin/86400.d0)
*       write(lp,*) 'new dsurfq = ',dsurfq
      endif
      if     (abs(nint((diagfq*86400.d0)/baclin)-
     &                 (diagfq*86400.d0)/baclin  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)')
     &    'error - ',
     &    'diagfq is not a whole number of baroclinic time steps'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      else  ! make it exact
*       write(lp,*) 'old diagfq = ',diagfq
        diagfq = nint((diagfq*86400.d0)/baclin)*(baclin/86400.d0)
*       write(lp,*) 'new diagfq = ',diagfq
      endif
      if     (abs(nint((meanfq*86400.d0)/baclin)-
     &                 (meanfq*86400.d0)/baclin  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)')
     &    'error - ',
     &    'meanfq is not a whole number of baroclinic time steps'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      else  ! make it exact
*       write(lp,*) 'old meanfq = ',meanfq
        meanfq = nint((meanfq*86400.d0)/baclin)*(baclin/86400.d0)
*       write(lp,*) 'new meanfq = ',diagfq
      endif
      if     (meanfq.ne.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)')
     &    'error - ',
     &    'meanfq must be 0.0 (i.e. not yet implemented)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (kknest.ne.kdm .and. nestfq.gt.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - kknest (dimensions.h) must be kdm when nesting'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
c --- 'hybrlx' = HYBGEN: inverse relaxation coefficient (time steps)
c                 (1.0 for no relaxation)
c --- 'hybflg' = hybrid  generator flag (0=T&S, 1=th&S,   2=th&T)
c --- 'advflg' = thermal advection flag (0=T&S, 1=th&S,   2=th&T)
c --- 'advtyp' = scalar  advection type (0=PCM, 1=MPDATA, 2=FCT2, 4=FCT4)
c --- 'slip'   = +1 for free-slip, -1  for non-slip boundary conditions
c --- 'visco2' = deformation-dependent Laplacian  viscosity factor
c --- 'visco4' = deformation-dependent biharmonic viscosity factor
c --- 'veldf2' = diffusion velocity (m/s) for Laplacian  momentum dissipation
c ---             (negative to input spacially varying diffusion velocity)
c --- 'veldf4' = diffusion velocity (m/s) for biharmonic momentum dissipation
c ---             (negative to input spacially varying diffusion velocity)
c --- 'thkdf2' = diffusion velocity (m/s) for Laplacian  thickness diffusion
c --- 'thkdf4' = diffusion velocity (m/s) for biharmonic thickness diffusion
c ---             (negative to input spacially varying diffusion velocity)
c --- 'temdf2' = diffusion velocity (m/s) for Laplacian  temp/saln diffusion
c --- 'temdfc' = temp diffusion conservation (0.0,1.0 all dens,temp resp.)
c --- 'vertmx' = diffusion velocity (m/s) for mom.mixing at mix.layr.base
c ---             (vertmx only used in MICOM-like isopycnal mode)
c --- 'cbar'   = rms flow speed     (m/s) for linear bottom friction
c --- 'cb'     = coefficient of quadratic bottom friction
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinr(hybrlx,'hybrlx','(a6," =",f10.4," time steps")')
      call blkini(hybflg,'hybflg')
      call blkini(advflg,'advflg')
      call blkini(advtyp,'advtyp')
      call blkinr(slip,  'slip  ',
     &            '(a6," =",f10.4," (-1=no-slip, +1=free-slip)")')
      call blkinr(visco2,'visco2','(a6," =",f10.4," ")')
      call blkinr(visco4,'visco4','(a6," =",f10.4," ")')
      call blkinr(veldf2,'veldf2','(a6," =",f10.4," m/s")')
      call blkinr(veldf4,'veldf4','(a6," =",f10.4," m/s")')
      call blkinr(thkdf2,'thkdf2','(a6," =",f10.4," m/s")')
      call blkinr(thkdf4,'thkdf4','(a6," =",f10.4,
     &                          " m/s (-ve if variable)")')
      call blkinr(temdf2,'temdf2','(a6," =",f10.4," m/s")')
      call blkinr(temdfc,'temdfc',
     &    '(a6," =",f10.4," (0.0,1.0 conserve dens,temp resp.)")')
      call blkinr(vertmx,'vertmx','(a6," =",f10.4," m/s")')
      call blkinr(cbar,  'cbar  ','(a6," =",f10.4," m/s")')
      call blkinr(cb,    'cb    ','(a6," =",f10.4," ")')
c
      if (hybrlx.lt.1.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - hybrlx must be at least 1.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      qhybrlx = 1.0/hybrlx
c
      if (hybflg.lt.0 .or. hybflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - hybflg must be 0 (T&S) or 1 (th&S) or 2 (th&T)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (advflg.lt.0 .or. advflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - advflg must be 0 (T&S) or 1 (th&S) or 2 (th&T)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (advflg.eq.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - advflg==2 (th&T) not yet implemented'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (advflg.ne.hybflg) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - advflg and hybflg must be identical'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (advtyp.ne.0 .and. advtyp.ne.1 .and.
     &    advtyp.ne.2 .and. advtyp.ne.4      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)') 
     &   'error - advtyp must be 0,1,2,4',
     &   ' (PCM,MPDATA,FCT2,FCT4)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (slip.ne.-1.0 .and. slip.ne.1.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - slip must be -1.0 (no-slip) or +1.0 (free-slip)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (thkdf2.ne.0.0 .and. thkdf4.ne.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - only one of thkdf2 and thkdf4 can be non-zero'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (isopyc .and. temdfc.ne.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - isopycnal mode must have temdfc=0.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (temdfc.lt.0.0 .or. temdfc.gt.1.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - temdfc must be between 0.0 and 1.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
c --- 'h1'     = depth interval used in lateral weighting of hor.pres.grad.
c --- 'thkbot' = thickness of bottom boundary layer (m)
c --- 'sigjmp' = minimum density jump across interfaces  (kg/m**3)
c --- 'tmljmp' = equivalent temperature jump across mixed-layer (degC)
c --- 'thkmls' = reference mixed-layer thickness for SSS relaxation (m)
c --- 'thkmlt' = reference mixed-layer thickness for SST relaxation (m)
c --- 'thkriv' = nominal thickness of river inflow (m)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      h1     = tenm
      call blkinr(thkbot,'thkbot','(a6," =",f10.4," m")')
      call blkinr(sigjmp,'sigjmp','(a6," =",f10.4," kg/m**3")')
      call blkinr(tmljmp,'tmljmp','(a6," =",f10.4," degC")')
      call blkinr(thkmls,'thkmls','(a6," =",f10.4," m")')
      call blkinr(thkmlt,'thkmlt','(a6," =",f10.4," m")')
      call blkinr(thkriv,'thkriv','(a6," =",f10.4," m")')
c
c --- 'iceflg' = ice model flag    (0=none,1=energy loan model)
c --- 'tfrz_0' = ENLN: ice melting point (degC) at S=0psu
c --- 'tfrz_s' = ENLN: gradient of ice melting point (degC/psu)
c --- 'ticegr' = ENLN: vertical temperature gradient inside ice (deg/m)
c ---                   (0.0 to get ice surface temp. from atmos. surtmp)
c --- 'hicemn' = ENLN: minimum ice thickness (m)
c --- 'hicemx' = ENLN: maximum ice thickness (m)
c
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(iceflg,'iceflg')
      icegln = iceflg.eq.1
c
      if (iceflg.lt.0 .or. iceflg.gt.1) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - iceflg must be between 0 and 1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
c --- to recover original EL model use:
c --- tfrz_0=-1.8, tfrz_s=0.0,   ticegr=2.0, hicemn=0.5, hicemx=10.0.
c --- for freezing point linear in S, typically use:
c --- tfrz_0=0.0, tfrz_s=-0.055, ticegr=0.0, hicemn=0.5, hicemx=10.0.
c
      call blkinr(tfrz_0,'tfrz_0','(a6," =",f10.4," degC")')
      call blkinr(tfrz_s,'tfrz_s','(a6," =",f10.4," degC/psu")')
      call blkinr(ticegr,'ticegr','(a6," =",f10.4," degC/m")')
      call blkinr(hicemn,'hicemn','(a6," =",f10.4," m")')
      call blkinr(hicemx,'hicemx','(a6," =",f10.4," m")')
c
c --- 'ntracr' = number of tracers (0=none,negative to initialize)
c --- 'trcflg' = tracer flags      (one digit per tracer, most sig. replicated)
c ---              0: passive, 100% at surface
c ---              1: passive, psudo-silicate
c ---              2: passive, temperature
c ---              3: passive
c ---            4-8: unused
c ---              9: default biology (NPZ-3,NPZD-4,Chai-9)
c
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(ntracr,'ntracr')
      trcrin = ntracr.gt.0  ! positive from restart, otherwise initialize
      trcout = ntracr.ne.0
      ntracr = abs(ntracr)
c
      if (ntracr.gt.mxtrcr) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3, a /)') 
     &    'error - maximum ntracr is',mxtrcr,
     &    '  (recompile with larger mxtrcr)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      call blkini(trcflg1,'trcflg')
      do k= 1,ntracr
        trcflg(k) = mod(trcflg1,10)  ! least significant decimal digit
        if     (trcflg1.ge.10) then
          trcflg1 = trcflg1/10  ! shift by one decimal digit
        else
          ! replicate last decimal digit across remaining tracers
        endif
        if (mnproc.eq.1) then
        write(lp,'(a,i3,i2)') '    k,trcflg =',k,trcflg(k)
        endif !1st tile
        if     (trcflg(k).gt.3 .and. trcflg(k).lt.9) then !not 0,1,2,3,9
          if (mnproc.eq.1) then
          write(lp,'(/ a,i3 /)') 
     &      'error - unknown tracer type for tracer',k
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      enddo !k
c
c --- 'tsofrq' = number of time steps between anti-drift offset calcs
c --- 'tofset' = temperature anti-drift offset (degC/century)
c --- 'sofset' = salnity     anti-drift offset  (psu/century)
c
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(tsofrq,'tsofrq')
      call blkinr(tofset,'tofset','(a6," =",f10.4," degC/century")')
      call blkinr(sofset,'sofset','(a6," =",f10.4,"  psu/century")')
c
c --- convert from per century to per second.
      tofset = tofset / (100.d0*365.d0*86400.d0)
      sofset = sofset / (100.d0*365.d0*86400.d0)
c
c --- 'mlflag' = mixed layer flag  (0=none,1=KPP,2-3=KTa-b,4=PWP,5=MY2.5,6=GISS)
c ---    'mxlkpp' = KPP:     activate mixed layer model (mlflag==1)
c ---    'mxlkrt' = KT:      MICOM or HYCOM Kraus-Turner (mlflag==2,3)
c ---    'mxlkta' = KT:      activate    original mixed layer model (mlflag==2)
c ---    'mxlktb' = KT:      activate alternative mixed layer model (mlflag==3)
c ---    'mxlpwp' = PWP:     activate mixed layer model (mlflag==4)
c ---    'mxlmy ' = MY:      activate mixed layer model (mlflag==5)
c ---   'mxlgiss' = GISS:    activate mixed layer model (mlflag==6)
c --- 'pensol' = KT/PWP   activate penetrating solar radiation (0=F,1=T)
c --- 'dtrate' = KT:      maximum permitted m.l. detrainment rate (m/day)
c --- 'thkmin' = KT/PWP:  minimum mixed-layer thickness (m)
c --- 'dypflg' = KT/PWP:  diapycnal mixing flag (0=none,1=KPP,2=explicit)
c ---                      always none (0) for MY/GISS and explicit (2) for PWP
c --- 'mixfrq' = KT/PWP:  number of time steps between diapycnal mixing calcs
c --- 'diapyc' = KT/PWP:  diapycnal diffusivity x buoyancy freq. (m**2/s**2)
c --- 'rigr'   = PWP:     critical gradient richardson number
c --- 'ribc'   = PWP:     critical bulk     richardson number
c --- 'rinfty' = KPP:     maximum  gradient richardson number (shear inst.)
c --- 'ricr'   = KPP:     critical bulk     richardson number
c --- 'bldmin' = KPP:     minimum surface boundary layer thickness (m)
c --- 'bldmax' = KPP:     maximum surface boundary layer thickness (m)
c --- 'cekman' = KPP/KT:  scale factor for Ekman depth
c --- 'cmonob' = KPP:     scale factor for Monin-Obukov depth
c --- 'bblkpp' = KPP:     activate bottom boundary layer    (0=F,1=T)
c --- 'shinst' = KPP:     activate shear instability mixing (0=F,1=T)
c --- 'dbdiff' = KPP:     activate double diffusion  mixing (0=F,1=T)
c --- 'nonloc' = KPP:     activate nonlocal b. layer mixing (0=F,1=T)
c --- 'latdiw' = KPROF:   activate lat.depen.int.wav mixing (0=F,1=T)
c --- 'botdiw' = GISS:    activate bot.enhan.int.wav mixing (0=F,1=T)
c --- 'difsmo' = KPROF:   activate horiz smooth diff coeffs (0=F,1=T)
c --- 'difout' = KPROF:   output visc/diff coffs in archive (0=F,1=T)
c --- 'difm0'  = KPP:     max viscosity   due to shear instability (m**2/s)
c --- 'difs0'  = KPP:     max diffusivity due to shear instability (m**2/s)
c --- 'difmiw' = KPP/MY:  background/internal wave viscosity       (m**2/s)
c --- 'difsiw' = KPP/MY:  background/internal wave diffusivity     (m**2/s)
c --- 'dsfmax' = KPP:     salt fingering diffusivity factor        (m**2/s)
c --- 'rrho0'  = KPP:     salt fingering rp=(alpha*delT)/(beta*delS)
c --- 'cs'     = KPP:     value for nonlocal flux term
c --- 'cstar'  = KPP:     value for nonlocal flux term
c --- 'cv'     = KPP:     buoyancy frequency ratio (0.0 to use a fn. of N)
c --- 'c11'    = KPP:     value for turb velocity scale
c --- 'hblflg' = KPP:     b. layer interpolation flag (0=con.,1=lin.,2=quad.)
c --- 'niter'  = KPP:     iterations for semi-implicit soln. (2 recomended)
c
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(mlflag,'mlflag')
c
      if (mlflag.lt.0 .or. mlflag.gt.6) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - mlflag must be between 0 and 6'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if (isopyc .and. mlflag.ne.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - isopycnal mode requires KT mixed layer (mlflag=2)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      mxlkta = mlflag.eq.2 .and. hybrid
      mxlktb = mlflag.eq.3 .and. hybrid
      mxlkrt = mlflag.eq.2 .or. mlflag.eq.3
      call blkinl(pensol,'pensol')
      call blkinr(dtrate,'dtrate','(a6," =",f10.4," m/day")')
      call blkinr(thkmin,'thkmin','(a6," =",f10.4," m")')
      call blkini(dypflg,'dypflg')
      call blkini(mixfrq,'mixfrq')
      call blkinr(diapyc,'diapyc','(a6," =",f10.4," m**2/s**2")')
c
      if (isopyc .and. pensol) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - isopycnal mode not consistent with pensol'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if (dypflg.lt.0 .or. dypflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - dypflg must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if ((mxlkta .or. mxlktb) .and.
     &    dypflg.ne.0 .and. thkriv.gt.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - thkriv>0 not consistent with KT unless dypflg=1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      mxlmy = mlflag.eq.5
c
      if (mxlmy) then
        if (kkmy25.ne.kdm) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     &      'error - kkmy25 (dimensions.h) must be kdm when mlflag==5'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        elseif (dypflg.ne.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     &      'warning - dypflg reset to 0 for M-Y 2.5'
          call flush(lp)
          endif !1st tile
          dypflg = 0
        endif
      endif
c
      mxlgiss = mlflag.eq.6
c
      if (mxlgiss) then
        if (nlgiss.ne.762) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     &      'error - nlgiss (dimensions.h) must be 762 when mlflag==6'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        elseif (dypflg.ne.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     &      'warning - dypflg reset to 0 for GISS'
          call flush(lp)
          endif !1st tile
          dypflg = 0
        endif
      endif
c
      mxlpwp = mlflag.eq.4
      call blkinr(rigc  ,'rigr  ','(a6," =",f10.4," ")')
      call blkinr(ribc  ,'ribc  ','(a6," =",f10.4," ")')
      call blkinr(qrinfy,'rinfty','(a6," =",f10.4," ")')
      qrinfy = 1.0/qrinfy
      call blkinr(ricr  ,'ricr  ','(a6," =",f10.4," ")')
c
      if (mxlpwp .and. dypflg.ne.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'warning - dypflg  reset to 2 for PWP mixed layer'
        call flush(lp)
        endif !1st tile
        dypflg = 2
      endif
c
      if (mxlpwp .and. thkriv.gt.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - thkriv>0 not consistent with PWP mixed layer'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      call blkinr(bldmin,'bldmin','(a6," =",f10.4," m")')
      call blkinr(bldmax,'bldmax','(a6," =",f10.4," m")')
c
      call blkinr(cekman,'cekman','(a6," =",f10.4," ")')
      call blkinr(cmonob,'cmonob','(a6," =",f10.4," ")')
c
      mxlkpp = mlflag.eq.1
      call blkinl(bblkpp,'bblkpp')
      call blkinl(shinst,'shinst')
      call blkinl(dbdiff,'dbdiff')
      call blkinl(nonloc,'nonloc')
      call blkinl(latdiw,'latdiw')
      call blkinl(botdiw,'botdiw')
      call blkinl(difsmo,'difsmo')
      call blkinl(difout,'difout')
      call blkinr(difm0 ,'difm0 ','(a6," =",f10.4," m**2/s")')
      call blkinr(difs0 ,'difs0 ','(a6," =",f10.4," m**2/s")')
      call blkinr(difmiw,'difmiw','(a6," =",f10.4," m**2/s")')
      call blkinr(difsiw,'difsiw','(a6," =",f10.4," m**2/s")')
      call blkinr(dsfmax,'dsfmax','(a6," =",f10.4," m**2/s")')
      call blkinr(rrho0 ,'rrho0 ','(a6," =",f10.4," ")')
      call blkinr(cs    ,'cs    ','(a6," =",f10.4," ")')
      call blkinr(cstar ,'cstar ','(a6," =",f10.4," ")')
      call blkinr(cv    ,'cv    ','(a6," =",f10.4," ")')
      call blkinr(c11   ,'c11   ','(a6," =",f10.4," ")')
      call blkini(hblflg,'hblflg')
      call blkini(niter ,'niter ')
c
      if (hblflg.lt.0 .or. hblflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ 2a /)') 
     &    'error - hblflg must be',
     &    ' 0 (constant) or 1 (linear) or 2 (quadratic)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if (mxlgiss) then
        if (botdiw .and. .not.latdiw) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     &      'error - botdiw requires latdiw'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      endif !mxlgiss
c
      if (mxlkpp) then
c ---   for KPP, diapyc and vertmx are not used
        dypflg = 0
        diapyc = 0.0
        vertmx = 0.0
      endif
      if (mxlmy) then
c ---   for M-Y, diapyc and vertmx are not used
        diapyc = 0.0
        vertmx = 0.0
      endif
      if (mxlgiss) then
c ---   for GISS, diapyc and vertmx are not used
        diapyc = 0.0
        vertmx = 0.0
      endif
      if (mxlpwp) then
c ---   for PWP, vertmx is not used
        vertmx = 0.0
      endif
      if (mxlkta .or. mxlktb) then
c ---   for HYCOM KT, vertmx is not currently used
        vertmx = 0.0
      endif
c
      if     (dypflg.eq.2) then
        if     (max(tofset,sofset).gt.0.0) then
          if     (diapyc.le.0.0) then
            if (mnproc.eq.1) then
            write(lp,'(/ a /)') 
     &        'error - diapyc must be positive if [ts]ofset is'
            call flush(lp)
            endif !1st tile
            call xcstop('(blkdat)')
                   stop '(blkdat)'
          elseif (mixfrq.ne.tsofrq) then
            if (mnproc.eq.1) then
            write(lp,'(/ a /)') 
     &        'error - mixfrq and tsofrq must be equal'
            call flush(lp)
            endif !1st tile
            call xcstop('(blkdat)')
                   stop '(blkdat)'
          endif !error checks
        endif ![ts]ofset>0.0
      endif !dypflg.eq.2
c
c --- weights for time smoothing
      wuv1  = 0.75
      wuv2  = 0.125
      wts1  = 0.875
      wts2  = 0.0625
      wbaro = 0.125
c
c --- 'spcifh' = specific heat of sea water (j/kg/deg)
c --- 'epsil'  = small nonzero number used to prevent division by zero
      spcifh=3990.
      epsil =1.0e-11
c
c --- 'clmflg' = climatology frequency flag (6=bimonthly,12=monthly)
c --- 'lbflag' = lateral barotropic bndy flag (0=none,1=port,2=input)
c --- 'wndflg' = wind stress input flag (0=none,1=u&v-grid,2,3=p-grid)
c ---             (=3 to calculate wind speed from wind stress)
c --- 'ustflg' = ustar forcing   flag          (3=input,1,2=wndspd,4=stress)
c --- 'flxflg' = thermal forcing flag   (0=none,3=net_flux,1,2,4=sst-based)
c --- 'empflg' = E-P     forcing flag   (0=none,3=net_E-P, 1,2,4=sst-based_E)
c --- 'sssflg' = SSS relaxation  flag   (0=none,1=clim)
c --- 'lwflag' = longwave corr.  flag   (0=none,1=clim,2=atmos), sst-based
c --- 'sstflg' = SST relaxation  flag   (0=none,1=clim,2=atmos,3=observed)
c --- 'icmflg' = ice mask        flag   (0=none,1=clim,2=atmos,3=observed)
c --- 'flxoff' = net flux offset flag   (0=F,1=T)
c --- 'flxsmo' = smooth surface fluxes  (0=F,1=T)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(clmflg,'clmflg')
      call blkini(lbflag,'lbflag')
      call blkini(wndflg,'wndflg')
      call blkini(ustflg,'ustflg')
      call blkini(flxflg,'flxflg')
      call blkini(empflg,'empflg')
      call blkini(sssflg,'sssflg')
      call blkini(lwflag,'lwflag')
      call blkini(sstflg,'sstflg')
      call blkini(icmflg,'icmflg')
      call blkinl(flxoff,'flxoff')
      call blkinl(flxsmo,'flxsmo')
c
      if (clmflg.ne.6 .and. clmflg.ne.12) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - clmflg must be 6 or 12'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (lbflag.lt.0 .or. lbflag.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - lbflag must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (lbflag.ne.2 .and. bnstfq.ne.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - bnstfq must be 0.0 unless lbflag is 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (lbflag.eq.2 .and. bnstfq.le.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - bnstfq must be positive when lbflag is 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (wndflg.lt.0 .or. wndflg.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - wndflg must be between 0 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (ustflg.lt.1 .or. ustflg.gt.4) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')
     &    'error - ustflg must be between 1 and 4'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.lt.0 .or. flxflg.gt.4) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - flxflg must be between 0 and 4'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.gt.0 .and.
     &    wndflg.eq.0      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - wndflg must be 1,2,3 when flxflg>0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.eq.3 .and.   
     &    ustflg.eq.4      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')
     &    'error - ustflg must be 1,2,3 when flxflg==3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif    
      if (empflg.lt.0 .or. empflg.gt.4) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - empflg must be between 0 and 4'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.eq.3 .and. (empflg.ne.0 .and. empflg.ne.3)) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - empflg must be 0 or 3 when flxflg==3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (empflg.ne.flxflg .and. empflg.ne.0 .and. empflg.ne.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - empflg must be 0 or 3 or flxflg'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.eq.0 .and.
     &    empflg.ne.0      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - empflg must be 0 when flxflg==0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if (lwflag.lt.0 .or. lwflag.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - lwflag must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (lwflag.gt.0 .and.
     &    flxflg.eq.0      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - flxflg must be >0 when lwflag>0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if (flxflg.eq.0) then
        flxoff = .false. !no thermal forcing
      endif
c
      if (iceflg.eq.0) then
        icmflg = 0  !no ice
      endif
c
      if (icmflg.lt.0 .or. icmflg.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - icmflg must be between 0 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      elseif (icmflg.eq.1 .or. icmflg.eq.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - icmflg 1 and 3 not yet implemented'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if (sstflg.lt.0 .or. sstflg.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - sstflg must be between 0 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (sstflg.gt.1 .and.
     &    yrflag.lt.2      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)') 
     &    'error - yrflag must be >1 (high frequency forcing)',
     &    ' when sstflg>1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if (sssflg.lt.0 .or. sssflg.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - sssflg must be 0 or 1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
c ---  s w i t c h e s    (if set to .true., then...)
c ---  (due to a SGI bug: read in an integer with 0=F,1=T)
c --- windf       use wind stress forcing (wndflg>0)
c --- thermo      use thermodynamic forcing 
c --- pensol      use penetrating solar radiation (input above)
c --- pcipf       use E-P forcing (may be redefined in forfun)
c
c --- relax       activate lateral boundary T/S/p  climatological nudging
c --- trcrlx      activate lateral boundary tracer climatological nudging
c
c --- srelax      activate surface salinity        climatological nudging
c --- trelax      activate surface temperature     climatological nudging
c --- relaxt      input tracer  climatological relaxation fields
c --- relaxf      input T/S/p   climatological relaxation fields
c ---              (note that relaxt implies relaxf)
c --- relaxs      input surface climatological relaxation fields only
c --- priver      rivers as a precipitation bogas
c --- epmass      treat evap-precip as a mass exchange
c
      windf  = wndflg.ne.0
      thermo = flxflg.ne.0
      pensol = pensol .and. thermo
      pcipf  = empflg.ne.0  ! if .true., might later be set .false. by forfun
c
      if (.not.(thermo .or. sstflg.gt.0 .or. srelax)) then
        niter=1
      elseif (mxlkta) then
        if (mnproc.eq.1) then
        write(lp,'(/ a / a /)') 
     &    'error - KT mixed layer needs thermal forcing, i.e.',
     &    '        mlflag=2 needs max(sstflg,sstflg,flxflg)>0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !.not(thermo ...):mxlkta
c
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinl(relax, 'relax ')
      call blkinl(trcrlx,'trcrlx')
      call blkinl(priver,'priver')
      call blkinl(epmass,'epmass')
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
c
      if     (priver .and. .not.thermo) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - priver must be .false. for flxflg=0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if     (epmass .and. .not.thermo) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - epmass must be .false. for flxflg=0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      close(unit=99)  !file='blkdat.input'
c
      srelax = sssflg.eq.1
      trelax = sstflg.eq.1
      relaxt = trcout .and. (trcrlx .or.
     &                       (.not.trcrin .and. iniflg.eq.2))
      relaxf =        relax
     &          .or.  srelax
     &          .or.  trelax
     &          .or.  lwflag.eq.1
     &          .or.  relaxt
      relaxs = .not. relax
     &         .and. relaxf
c
c --- initialize from climatology (update relaxf and relaxs)?
      if     (iniflg.eq.2 .and. yrflag.ne.3) then
        open(unit=99,file='limits')
        read(99,*) day1
        close(unit=99) !file='limits'
        if     (day1.le.0.0) then
          relaxf = .true.
          relaxs = .false.
        endif !initialize from climatology
      endif
c
      if (kkwall.ne.kdm .and. relaxf .and. .not. relaxs) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     &    'error - kkwall (dimensions.h) must be kdm for 3-d clim'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
c --- use 'huge' to initialize array portions that the code should never access
      huge = 2.0**100  ! 2^100, or about 1.2676506e30
c
c --- i/o file names
c
      flnmdep  = 'regional.depth'
      flnmgrd  = 'regional.grid'
      flnmrsi  = 'restart_in'
      flnmrso  = 'restart_out'
      flnmarc  = 'archv.0000_000_00'
      flnmovr  = 'ovrtn_out'
      flnmflx  = 'flxdp_out'
c
c --- i/o directory names
c
      flnmfor  = './'
      flnmforw = './'
      return
      end
      subroutine blkinr(rvar,cvar,cfmt)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      real      rvar
      character cvar*6,cfmt*(*)
c
c     read in one real value
c
      character*6 cvarin
c
      read(99,*) rvar,cvarin
      if (mnproc.eq.1) then
      write(lp,cfmt) cvarin,rvar
      call flush(lp)
      endif !1st tile
c
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkinr)')
               stop
      endif
      return
      end
      subroutine blkini(ivar,cvar)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value
c
      character*6 cvarin
c
      read(99,*) ivar,cvarin
      if (mnproc.eq.1) then
      write(lp,6000) cvarin,ivar
      call flush(lp)
      endif !1st tile
c
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkini - input ',cvarin,
     +                      ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkini)')
               stop
      endif
      return
 6000 format(a6,' =',i6)
      end
      subroutine blkinl(lvar,cvar)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      logical     lvar
      character*6 cvar
c
c     read in one logical value
c     due to a SGI bug for logical I/O: read in an integer 0=F,1=T
c
      character*6 cvarin
      integer     ivar
c
      read(99,*) ivar,cvarin
      lvar = ivar .ne. 0
      if (mnproc.eq.1) then
      write(lp,6000) cvarin,lvar
      call flush(lp)
      endif !1st tile
c
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkinl - input ',cvarin,
     +                      ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkinl)')
               stop
      endif
      return
 6000 format(a6,' =',l6)
      end
c>
c> Revision history
c>
c> Oct. 1999 - added variables to control penetrating solar radiation
c> Oct. 1999 - added switch to select mixed layer model
c> Oct. 1999 - dp00 for hybgen is now set here
c> Dec. 1999 - multiple heat flux transfer coefficients (cts1, cts2, ctl)
c>             replace old coefficient ct
c> Jan. 2000 - changed to subroutine with run-time input
c> May. 2000 - conversion to SI units
c> Nov. 2000 - added kapflg,thflag,hybflg,advflg,wndflg
c> Dec. 2000 - added flxflg
c> Aug. 2001 - added bnstfq,nestfq
c> Aug. 2001 - added mapflg==4 for an f-plane
c> Aug. 2001 - betard and betabl inverted to replace division by multiply
c> May  2002 - map projection via regional.grid file (see geopar.f)
c> May  2002 - vertical coordinate via d[ps]00*
c> May  2002 - diffusion variable names include 2/4 for Laplacian/biharmonic
c> May  2002 - added PWP and MY2.5 mixed layer options
c> May  2002 - added thkmlr, distinct from thkmin
c> Aug  2002 - added ntracr and trcflg to control tracers
c> Nov  2002 - added jerlv0=0 to use kpar-based turbidity
c> Nov  2002 - split thkmlr into thkmls and thkmlt
c> Apr  2003 - added dp00i, vsigma, and priver
c> May  2003 - added bldmin, bldmax, flxsmo, and icmflg
c> Jun  2003 - added locsig and removed thflag=4
c> Oct  2003 - thkdf4 negative now signals spacial variation
c> Nov  2003 - added advtyp
c> Jan  2004 - added latdiw
c> Jan  2004 - added bblkpp
c> Jan  2004 - added hblflg and cv=0.0 option
c> Feb  2004 - added botdiw (and latdiw) for GISS
c> Feb  2004 - added temdfc
c> Mar  2004 - added thkriv and epmass
c> Mar  2004 - added isotop
c> Mar  2005 - added tfrz_0,tfrz_s,ticegr,hicemn,hicemx
c> Mar  2005 - added tsofrq,tofset,sofset
c> Mar  2005 - added empflg
c> Mar  2005 - added ustflg, reordered thermal forcing flags
c> Mar  2005 - added flxoff
c> Apr  2005 - replaced kapflg with kapref
c> Jun  2005 - added hybrlx