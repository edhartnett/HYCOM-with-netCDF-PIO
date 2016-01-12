      subroutine blkdat
      use mod_xc  ! HYCOM communication interface
c
      include 'common_blocks.h'
c
      real      grdlat,grdlon,pntlat,pntlon,reflat,reflon
      integer   k,mlflag,thflag,kapflg
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
c --- pi-related values
      radian=57.2957795
      pi    = 3.1415926536
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
      if (iversn.lt.20 .or. iversn.gt.20) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3,a,i3 /)') 
     .    'error - iversn must be between',20,' and',20
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
c --- s e l e c t  m a p   p r o j e c t i o n
c
c --- 'mapflg' = map flag (0=mercator,2=uniform,3=beta-plane,4=input)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(mapflg,'mapflg')
c
      if     (mapflg.eq.0) then
c ---   square mercator projection
c ---   'idm   ' = longitudinal array size
c ---   'pntlon' = longitudinal reference grid point on pressure grid
c ---   'reflon' = longitude of reference grid point on pressure grid
c ---   'grdlon' = longitudinal grid size (degrees)
c ---   'jdm   ' = latitudinal  array size
c ---   'pntlat' = latitudinal  reference grid point on pressure grid
c ---   'reflat' = latitude of  reference grid point on pressure grid
c ---   'grdlat' = latitudinal  grid size at the equator (degrees)
        call blkini(itest ,'idm   ')
        call blkinr(pntlon,'pntlon','(a6," =",f10.4," ")')
        call blkinr(reflon,'reflon','(a6," =",f10.4," deg E")')
        call blkinr(grdlon,'grdlon','(a6," =",f10.4," degrees")')
        call blkini(jtest ,'jdm   ')
        call blkinr(pntlat,'pntlat','(a6," =",f10.4," ")')
        call blkinr(reflat,'reflat','(a6," =",f10.4," deg N")')
        call blkinr(grdlat,'grdlat','(a6," =",f10.4," degrees")')
c
c ---   'ypivn'  = the j index of the equator on the pressure grid
c ---   'gridn'  = latitudinal  grid size (degrees)
        ypivn = pntlat  ! since reflat==0
        gridn = grdlat
c
        if     (itest.ne.itdm) then
          if (mnproc.eq.1) then
          write(lp,'(/ a,i5 /)') 
     .      'error - expected idm =',itdm
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
        if     (jtest.ne.jtdm) then
          if (mnproc.eq.1) then
          write(lp,'(/ a,i5 /)') 
     .      'error - expected jdm =',jtdm
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
        if     (grdlon.ne.grdlat) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     .      'error - grdlon.ne.grdlat'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
        if     (reflat.ne.0.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     .      'error - reflat must be the equator'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      elseif (mapflg.eq.2) then
c ---   uniform latitude grid
c ---   square mercator projection
c ---   'idm   ' = longitudinal array size
c ---   'pntlon' = longitudinal reference grid point on pressure grid
c ---   'reflon' = longitude of reference grid point on pressure grid
c ---   'grdlon' = longitudinal grid size (degrees)
c ---   'jdm   ' = latitudinal  array size
c ---   'pntlat' = latitudinal  reference grid point on pressure grid
c ---   'reflat' = latitude of  reference grid point on pressure grid
c ---   'grdlat' = latitudinal  grid size at the equator (degrees)
        call blkini(itest ,'idm   ')
        call blkinr(pntlon,'pntlon','(a6," =",f10.4," ")')
        call blkinr(reflon,'reflon','(a6," =",f10.4," deg E")')
        call blkinr(grdlon,'grdlon','(a6," =",f10.4," degrees")')
        call blkini(jtest ,'jdm   ')
        call blkinr(pntlat,'pntlat','(a6," =",f10.4," ")')
        call blkinr(reflat,'reflat','(a6," =",f10.4," deg N")')
        call blkinr(grdlat,'grdlat','(a6," =",f10.4," degrees")')
c
c ---   ypivn = the j index of the equator on the pressure grid
c ---   gridn = latitudinal  grid size (degrees)
c ---   grido = longitudinal grid size (degrees)
        ypivn = pntlat
        gridn = grdlat
        grido = grdlon
c
        if     (itest.ne.itdm) then
          if (mnproc.eq.1) then
          write(lp,'(/ a,i5 /)') 
     .      'error - expected idm =',itdm
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
        if     (jtest.ne.jtdm) then
          if (mnproc.eq.1) then
          write(lp,'(/ a,i5 /)') 
     .      'error - expected jdm =',jtdm
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
        if     (reflat.ne.0.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     .      'error - reflat must be the equator'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      elseif (mapflg.eq.3) then
c ---   square beta-plane grid
c ---   'idm   ' = longitudinal array size
c ---   'pntlon' = longitudinal reference grid point on pressure grid
c ---   'reflon' = longitude of reference grid point on pressure grid
c ---   'grdlon' = longitudinal grid size (degrees)
c ---   'jdm   ' = latitudinal  array size
c ---   'pntlat' = latitudinal  reference grid point on pressure grid
c ---   'reflat' = latitude of  reference grid point on pressure grid
c ---   'grdlat' = latitudinal  grid size at the equator (degrees)
        call blkini(itest ,'idm   ')
        call blkinr(pntlon,'pntlon','(a6," =",f10.4," ")')
        call blkinr(reflon,'reflon','(a6," =",f10.4," deg E")')
        call blkinr(grdlon,'grdlon','(a6," =",f10.4," degrees")')
        call blkini(jtest ,'jdm   ')
        call blkinr(pntlat,'pntlat','(a6," =",f10.4," ")')
        call blkinr(reflat,'reflat','(a6," =",f10.4," deg N")')
        call blkinr(grdlat,'grdlat','(a6," =",f10.4," degrees")')
c
c ---   'ypivn'  = the j index of the equator on the pressure grid
c ---   'gridn'  = latitudinal  grid size (degrees)
        ypivn = pntlat
        gridn = grdlat
c
        if     (itest.ne.itdm) then
          if (mnproc.eq.1) then
          write(lp,'(/ a,i5 /)') 
     .      'error - expected idm =',itdm
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
        if     (jtest.ne.jtdm) then
          if (mnproc.eq.1) then
          write(lp,'(/ a,i5 /)') 
     .      'error - expected jdm =',jtdm
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
        if     (reflat.ne.0.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     .      'error - reflat must be the equator'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      elseif (mapflg.eq.4) then
c ---   input a full grid of latitude and longitude values (see geopar)
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
        call blkini(itest ,'idm   ')
        call blkini(jtest ,'jdm   ')
c
        if     (itest.ne.itdm) then
          if (mnproc.eq.1) then
          write(lp,'(/ a,i5 /)') 
     .      'error - expected idm =',itdm
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
        if     (jtest.ne.jtdm) then
          if (mnproc.eq.1) then
          write(lp,'(/ a,i5 /)') 
     .      'error - expected jdm =',jtdm
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      else
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - unsupported mapflg value'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
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
     .    'error - maximum itest is',itdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (jttest.gt.jtdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i5 /)') 
     .    'error - maximum jtest is',jtdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
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
*     enddo
c
c --- 'kdm   ' = number of layers
c --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'dp00s'  = sigma   spacing minimum thickness (m)
c --- 'dp00'   = z-level spacing minimum thickness (m)
c --- 'dp00x'  = z-level spacing maximum thickness (m)
c --- 'dp00f'  = z-level spacing stretching factor (1.0=const.spacing)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(k,     'kdm   ')
      call blkini(nhybrd,'nhybrd')
      call blkini(nsigma,'nsigma')
      call blkinr(dp00s, 'dp00s ','(a6," =",f10.4," m")')
      call blkinr(dp00,  'dp00  ','(a6," =",f10.4," m")')
      call blkinr(dp00x, 'dp00x ','(a6," =",f10.4," m")')
      call blkinr(dp00f, 'dp00f ','(a6," =",f10.4," ")')
c
c --- isopycnal (MICOM-like) iff nhybrd is 0
      isopyc = nhybrd .eq. 0
      hybrid = .not. isopyc
      if (hybrid .and. nsigma.le.1) then
        nsigma=1
        dp00s=9999.0
      endif
c
      if     (k.ne.kdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3 /)') 
     .    'error - expected kdm =',kdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (nhybrd.gt.kdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3 /)') 
     .    'error - maximum nhybrd is kdm =',kdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (nsigma.gt.nhybrd) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3 /)') 
     .    'error - maximum nsigma is nhybrd =',nhybrd
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (isopyc .and. max(dp00,dp00x).ne.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - must have dp00x==dp00==0.0 for isopycnal case'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (dp00f.eq.1.0 .and. dp00.ne.dp00x) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - must have dp00x==dp00 for dp00f==1.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (dp00.gt.dp00x) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - dp00x must be at least dp00'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
c --- 'saln0'  = initial salinity value (psu), only used for iniflg<2
c --- 'kapflg' = thermobaric compressibility flag (-1=none or thflag)
c --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2, 4=Sigma-4)
c ---            this is a check on the compile-time stmt_funcs.h setup.
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinr(saln0, 'saln0 ','(a6," =",f10.4," psu")')
      call blkini(kapflg,'kapflg')
      call blkini(thflag,'thflag')
c --- tbaric is thermobaric compressibility switch
      if     (kapflg.eq.-1) then
        tbaric=.false.
      elseif (kapflg.eq.thflag) then
        tbaric=.true.
        if     (kapflg.eq.0 .and. pref.ne.0.0    ) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     .      'error - kapflg not consistent with pref'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        elseif (kapflg.eq.2 .and. pref.ne.2000.e4) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     .      'error - kapflg not consistent with pref'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        elseif (kapflg.eq.4 .and. pref.ne.4000.e4) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     .      'error - kapflg not consistent with pref'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      else
        if (mnproc.eq.1) then
        write(lp,'(/ a,i1 /)') 
     .    'error - kapflg must be -1 or ',thflag
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (thflag.eq.0) then
        if     (abs(sig(0.0,0.0)+0.136471).gt.0.001) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     .      'error - thflag not consistent with sig()'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      elseif (thflag.eq.2) then
        if     (abs(sig(0.0,0.0)-9.77093 ).gt.0.001) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     .      'error - thflag not consistent with sig()'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      elseif (thflag.eq.4) then
        if     (abs(sig(0.0,0.0)-19.2362 ).gt.0.001) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') 
     .      'error - thflag not consistent with sig()'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      else
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - thflag must be 0 or 2 or 4'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (thflag.eq.0) then
        sigfmt = '(a6," =",f10.4," sigma-0")'
      elseif (thflag.eq.2) then
        sigfmt = '(a6," =",f10.4," sigma-2")'
      elseif (thflag.eq.4) then
        sigfmt = '(a6," =",f10.4," sigma-4")'
      endif
c
c --- 'thbase' = reference density (sigma units)
      call blkinr(thbase,'thbase',sigfmt)
c
c --- layer densities (sigma units)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      do k=1,kdm
        call blkinr(sigma(k),'sigma ',sigfmt)
c
        if     (k.gt.1) then
          if      (sigma(k).le.sigma(k-1)) then
            if (mnproc.eq.1) then
            write(lp,'(/ a,i3 /)') 
     .        'error - sigma is not stabally stratified'
            call flush(lp)
            endif !1st tile
            call xcstop('(blkdat)')
                   stop '(blkdat)'
          endif
        endif
      enddo
c
c --- 'iniflg' = initial state flag (0=level,1=zonal,2=clim.,3=archv)
c --- 'jerlv0' = initial jerlov water type (1 to 5)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(iniflg,'iniflg')
      call blkini(jerlv0,'jerlv0')
c
      if (iniflg.lt.0 .or. iniflg.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - iniflg must be between 0 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (jerlv0.lt.1 .or. jerlv0.gt.5) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - jerlv0 must be between 1 and 5'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
c --- red and blue light extinction coefficients (pressure units)
c --- for jerlov water types 1 to 5 - fraction of penetrating red light
      betard(1) =  0.35*onem
      betard(2) =  0.6 *onem
      betard(3) =  1.0 *onem
      betard(4) =  1.5 *onem
      betard(5) =  1.4 *onem
      betabl(1) = 23.0 *onem
      betabl(2) = 20.0 *onem
      betabl(3) = 17.0 *onem
      betabl(4) = 14.0 *onem
      betabl(5) =  7.9 *onem
      redfac(1) = 0.58
      redfac(2) = 0.62
      redfac(3) = 0.67
      redfac(4) = 0.77
      redfac(5) = 0.78
c
c --- 'yrflag' = days in year flag (0=360,1=366,2=366Jan1,3=actual)
c --- 'dsurfq' = number of days between model diagnostics at the surface
c --- 'diagfq' = number of days between model diagnostics
c --- 'rstrfq' = number of days between model restart output
c --- 'baclin' = baroclinic time step (seconds), int. divisor of 86400
c --- 'batrop' = barotropic time step (seconds), int. divisor of baclin/2
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(yrflag,'yrflag')
      call blkinr(dsurfq,'dsurfq','(a6," =",f10.4," days")')
      call blkinr(diagfq,'diagfq','(a6," =",f10.4," days")')
      call blkinr(rstrfq,'rstrfq','(a6," =",f10.4," days")')
      call blkinr(baclin,'baclin','(a6," =",f10.4," sec")')
      call blkinr(batrop,'batrop','(a6," =",f10.4," sec")')
c
      if (yrflag.lt.0 .or. yrflag.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - yrflag must be between 0 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (yrflag.le.1) then  ! monthly forcing
        if     (abs(nint(86400.0/baclin)-86400.0/baclin).gt.0.01) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')
     .      'error - baclin not an integer divisor of 86400 secs'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      else  ! high frequency forcing
        if     (abs(nint(21600.0/baclin)-21600.0/baclin).gt.0.01) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')
     .      'error - baclin not an integer divisor of 21600 secs'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      endif
      if     (abs(nint(0.5*baclin/batrop)-
     .                 0.5*baclin/batrop  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')
     .    'error - batrop not an integer divisor of baclin/2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
c --- 'hybflg' = hybrid generator  flag (0=T&S, 1=th&S, 2=th&T)
c --- 'advflg' = thermal advection flag (0=T&S, 1=th&S, 2=th&T)
c --- 'slip'   = +1 for free-slip, -1  for non-slip boundary conditions
c --- 'biharm' = fraction of diffusion that is biharmonic (0.0 to 1.0)
c --- 'viscos' = deformation-dependent viscosity factor (nondimensional)
c --- 'veldff' = diffusion velocity (m/s) for momentum dissipation
c --- 'thkdff' = diffusion velocity (m/s) for thickness diffusion
c --- 'temdff' = diffusion velocity (m/s) for temp/saln diffusion
c --- 'vertmx' = diffusion velocity (m/s) for mom.mixing at mix.layr.base
c ---             (vertmx only used in MICOM-like isopycnal mode)
c --- 'cbar'   = rms flow speed     (m/s) for linear bottom friction
c --- 'cb'     = coefficient of quadratic bottom friction
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(hybflg,'hybflg')
      call blkini(advflg,'advflg')
      call blkinr(slip,  'slip  ',
     .            '(a6," =",f10.4," (-1=no-slip, +1=free-slip)")')
      call blkinr(biharm,'biharm',
     .            '(a6," =",f10.4," (0.0=Laplacian, 1.0=biharmonic)")')
      call blkinr(viscos,'viscos','(a6," =",f10.4," ")')
      call blkinr(veldff,'veldff','(a6," =",f10.4," m/s")')
      call blkinr(thkdff,'thkdff','(a6," =",f10.4," m/s")')
      call blkinr(temdff,'temdff','(a6," =",f10.4," m/s")')
      call blkinr(vertmx,'vertmx','(a6," =",f10.4," m/s")')
      call blkinr(cbar,  'cbar  ','(a6," =",f10.4," m/s")')
      call blkinr(cb,    'cb    ','(a6," =",f10.4," ")')
c
      if (hybflg.lt.0 .or. hybflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - hybflg must be 0 (T&S) or 1 (th&S) or 2 (th&T)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (advflg.lt.0 .or. advflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - advflg must be 0 (T&S) or 1 (th&S) or 2 (th&T)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (advflg.eq.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - advflg==2 (th&T) not yet implemented'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (slip.ne.-1.0 .and. slip.ne.1.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - slip must be -1.0 (no-slip) or +1.0 (free-slip)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (biharm.lt.0.0 .or. biharm.gt.1.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - biharm must be between 0.0 and 1.0'
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
c --- 'thkmin' = minimum mixed-layer thickness (m)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      h1     = tenm
      call blkinr(thkbot,'thkbot','(a6," =",f10.4," m")')
      call blkinr(sigjmp,'sigjmp','(a6," =",f10.4," kg/m**3")')
      call blkinr(tmljmp,'tmljmp','(a6," =",f10.4," degC")')
      call blkinr(thkmin,'thkmin','(a6," =",f10.4," m")')
c
c --- 'iceflg' = ice model flag (0=none,1=energy loan model)
c --- 'mlflag' = mixed layer flag (0=none,1=KPP,2=KTa,3=KTb)
c --- 'mxlkta' = KT:  activate    original mixed layer model (mlflag==2)
c --- 'mxlktb' = KT:  activate alternative mixed layer model (mlflag==3)
c --- 'mxlkrt' = KT:  MICOM or HYCOM Kraus-Turner (mlflag==2,3)
c --- 'pensol' = KT:  activate penetrating solar radiation
c --- 'dypflg' = KT: diapycnal mixing flag (0=none,1=KPP,2=explicit)
c --- 'mixfrq' = KT: number of time steps between diapycnal mixing calcs
c --- 'diapyc' = KT: diapycnal diffusivity x buoyancy freq. (m**2/s**2)
c --- 'dtrate' = KT: maximum permitted m.l. detrainment rate (m/day)
c --- 'mxlkpp' = KPP: activate mixed layer model (mlflag==1)
c --- 'shinst' = KPP: activate shear instability mixing (0=F,1=T)
c --- 'dbdiff' = KPP: activate double diffusion  mixing (0=F,1=T)
c --- 'nonloc' = KPP: activate nonlocal b. layer mixing (0=F,1=T)
c --- 'difsmo' = KPP: activate horiz smooth diff coeffs (0=F,1=T)
c --- 'rinfty' = KPP: value for calculating rshear instability
c --- 'difm0'  = KPP: max viscosity   due to shear instability (m**2/s)
c --- 'difs0'  = KPP: max diffusivity due to shear instability (m**2/s)
c --- 'difmiw' = KPP: background/internal wave viscosity       (m**2/s)
c --- 'difsiw' = KPP: background/internal wave diffusivity     (m**2/s)
c --- 'dsfmax' = KPP: salt fingering diffusivity factor        (m**2/s)
c --- 'rrho0'  = KPP: salt fingering rp=(alpha*delT)/(beta*delS)
c --- 'ricr'   = KPP: critical bulk richardson number
c --- 'cs'     = KPP: value for nonlocal flux term
c --- 'cstar'  = KPP: value for nonlocal flux term
c --- 'cv'     = KPP: value for turb shear contributon to bulk rich. no.
c --- 'c11'    = KPP: value for turb velocity scale
c --- 'niter'  = KPP: iterations for semi-implicit soln. (2 recomended)
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
     .    'error - iceflg must be between 0 and 1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(mlflag,'mlflag')
c
      if (mlflag.lt.0 .or. mlflag.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - mlflag must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if (isopyc .and. mlflag.eq.1) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - isopycnal mode not consistent with KPP mixed layer'
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
      call blkini(dypflg,'dypflg')
      call blkini(mixfrq,'mixfrq')
      call blkinr(diapyc,'diapyc','(a6," =",f10.4," m**2/s**2")')
      call blkinr(dtrate,'dtrate','(a6," =",f10.4," m/day")')
c
      if (isopyc .and. pensol) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - isopycnal mode not consistent with pensol'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      if (dypflg.lt.0 .or. dypflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - dypflg must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
c
      mxlkpp = mlflag.eq.1
      call blkinl(shinst,'shinst')
      call blkinl(dbdiff,'dbdiff')
      call blkinl(nonloc,'nonloc')
      call blkinl(difsmo,'difsmo')
      call blkinr(rinfty,'rinfty','(a6," =",f10.4," ")')
      call blkinr(difm0 ,'difm0 ','(a6," =",f10.4," m**2/s")')
      call blkinr(difs0 ,'difs0 ','(a6," =",f10.4," m**2/s")')
      call blkinr(difmiw,'difmiw','(a6," =",f10.4," m**2/s")')
      call blkinr(difsiw,'difsiw','(a6," =",f10.4," m**2/s")')
      call blkinr(dsfmax,'dsfmax','(a6," =",f10.4," m**2/s")')
      call blkinr(rrho0 ,'rrho0 ','(a6," =",f10.4," ")')
      call blkinr(ricr  ,'ricr  ','(a6," =",f10.4," ")')
      call blkinr(cs    ,'cs    ','(a6," =",f10.4," ")')
      call blkinr(cstar ,'cstar ','(a6," =",f10.4," ")')
      call blkinr(cv    ,'cv    ','(a6," =",f10.4," ")')
      call blkinr(c11   ,'c11   ','(a6," =",f10.4," ")')
      call blkini(niter ,'niter ')
c
      if (mxlkpp) then
c ---   for KPP, diapyc and vertmx are not used
        dypflg = 2
        diapyc = 0.0
        vertmx = 0.0
      endif
      if (mxlkta .or. mxlktb) then
c ---   for HYCOM KT, vertmx is not currently used
        vertmx = 0.0
      endif
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
c --- 'wndflg' = wind stress input flag (0=none,1=on u/v grid,2=on p grid)
c --- 'flxflg' = thermal forcing flag (0=none,1=orig,2=new-flux-calc)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(clmflg,'clmflg')
      call blkini(lbflag,'lbflag')
      call blkini(wndflg,'wndflg')
      call blkini(flxflg,'flxflg')
c
      if (clmflg.ne.6 .and. clmflg.ne.12) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - clmflg must be 6 or 12'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (lbflag.lt.0 .or. lbflag.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - lbflag must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (wndflg.lt.0 .or. wndflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - wndflg must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.lt.0 .or. flxflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 
     .    'error - flxflg must be between 0 and 2'
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
c --- relax       activate lateral boundary nudging
c --- srelax      activate surface salinity nudging
c --- trelax      activate surface temperature nudging
c --- relaxf      input relaxation fields
c
      windf = wndflg.ne.0
      thermo= flxflg.ne.0
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinl(relax, 'relax ')
      call blkinl(srelax,'srelax')
      call blkinl(trelax,'trelax')
      relaxf = relax .or. srelax .or. trelax .or. iniflg.eq.2
      pensol = pensol .and. thermo .and. hybrid
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
c
c --- use 'huge' to initialize array portions that the code should never access
      huge = 2.0**100  ! 2^100, or about 1.2676506e30
c
c --- i/o file names
c
      flnmdep  = 'regional.depth'
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
c
      close(unit=99)
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
