      subroutine forday(dtime,yrflag, iyear,iday,ihour)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c
      real*8  dtime
      integer yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,julian-day,hour).
c
      real*8  dtim1,day
      integer iyr,nleap
c
      if     (yrflag.eq.0) then
c ---   360 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/360.d0) + 1
        iday  =  mod( dtime+15.001d0 ,360.d0) + 1
        ihour = (mod( dtime+15.001d0 ,360.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.1) then
c ---   366 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/366.d0) + 1
        iday  =  mod( dtime+15.001d0 ,366.d0) + 1
        ihour = (mod( dtime+15.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.2) then
c ---   366 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/366.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,366.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
        iyr   = (dtime-1.d0)/365.25d0
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
        day   = dtime - dtim1 + 1.d0
        if     (dtim1.gt.dtime) then
          iyr = iyr - 1
        elseif (day.ge.367.d0) then
          iyr = iyr + 1
        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
c
        iyear =  1901 + iyr
        iday  =  dtime - dtim1 + 1
        ihour = (dtime - dtim1 + 1.d0 - iday)*24.d0
c
      else
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in forday - unsupported yrflag value'
        write(lp,*)
        endif !1st tile
        call xcstop('(forday)')
               stop '(forday)'
      endif
      return
      end
c
c
      subroutine forfuna
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
      include 'common_blocks.h'
c
c --- initialize input of atmospheric forcing fields
c
c --- units of tau_x  are N/sq_m (positive eastwards)
c --- units of tau_y  are N/sq_m (positive northwards)
c --- units of wndspd are m/s
c --- units of airtmp are degC
c --- units of vapmix are kg/kg
c --- units of precip are m/s    (postive into ocean)
c --- units of radflx are w/sq_m (postive into ocean)
c --- units of swflx  are w/sq_m (postive into ocean)
c
c --- tau_x and tau_y are either on u&v grids or both on the p grid,
c --- depending on the value of blkdat input parameter "wndflg".
c --- all other fields, including wndspd, are always on the p grid.
c
c --- I/O and array I/O units 901-908 are reserved for the entire run.
c
c --- all input fields much be defined at all grid points
c
      integer         mreca,mrecr
      common/rdforfi/ mreca,mrecr
      save  /rdforfi/
c
      integer   lgth,mo
      character preambl(5)*79
      real      pcmax,one,oneps
      integer   i,j
c
 103  format (a79)
 102  format (/(a79))
c
      mreca=1
      if     (mnproc.eq.1) then
      write (lp,*) ' now opening forcing fields ...'
      call flush(lp)
      endif !1st tile
c
      do 12 lgth=28,1,-1
        if (flnmfor(lgth:lgth).ne.' ') then
          goto 112
        endif
  12  continue
      call xcstop('(flnmfor)')
             stop '(flnmfor)'
 112  continue
c
      if (windf) then
c
      call zaiopf(flnmfor(1:lgth)//'forcing.tauewd.a', 'old', 901)
      open (unit=901,file=flnmfor(1:lgth)//'forcing.tauewd.b',
     &      status='old', action='read')
       read (901,103) preambl
      if     (mnproc.eq.1) then
      write (lp, 102) preambl
      call flush(lp)
      endif !1st tile
      call rdmonth(util1, 901)
cdiag call prtmsk(ip,util1,util2,idm,idm,jdm,  0.,1000.,
cdiag.     'tau_x (x 1000 N/sq_m) ')
c
      call zaiopf(flnmfor(1:lgth)//'forcing.taunwd.a', 'old', 902)
      open (unit=902,file=flnmfor(1:lgth)//'forcing.taunwd.b',
     &   status='old', action='read')
       read (902,103) preambl
      if     (mnproc.eq.1) then
      write (lp,102) preambl
      call flush(lp)
      endif !1st tile
      call rdmonth(util1, 902)
cdiag call prtmsk(ip,util1,util2,idm,idm,jdm,  0.,1000.,
cdiag.     'tau_y (x 1000 N/sq_m) ')
c
      end if				!  windf = .true.
c
      if (thermo) then
c
      call zaiopf(flnmfor(1:lgth)//'forcing.wndspd.a', 'old', 903)
      open (unit=903,file=flnmfor(1:lgth)//'forcing.wndspd.b',
     &   status='old', action='read')
       read (903,103) preambl
      if     (mnproc.eq.1) then
      write (lp,102) preambl
      call flush(lp)
      endif !1st tile
      call rdmonth(util1, 903)
cdiag call prtmsk(ip,util1,util2,idm,idm,jdm,  0.,10.,
cdiag.     'wind speed  (x 10 m/s)')
c
      call zaiopf(flnmfor(1:lgth)//'forcing.airtmp.a', 'old', 904)
      open (unit=904,file=flnmfor(1:lgth)//'forcing.airtmp.b',
     &   status='old', action='read')
       read (904,103) preambl
      if     (mnproc.eq.1) then
      write (lp,102) preambl
      call flush(lp)
      endif !1st tile
      call rdmonth(util1, 904)
cdiag call prtmsk(ip,util1,util2,idm,idm,jdm,  0.,10.,
cdiag.     'air temperature  (0.1 c)')
c
      call zaiopf(flnmfor(1:lgth)//'forcing.vapmix.a', 'old', 905)
      open (unit=905,file=flnmfor(1:lgth)//'forcing.vapmix.b',
     &   status='old', action='read')
       read (905,103) preambl
      if     (mnproc.eq.1) then
      write (lp,102) preambl
      call flush(lp)
      endif !1st tile
      call rdmonth(util1, 905)
cdiag call prtmsk(ip,util1,util2,idm,idm,jdm,  0.,10000.,
cdiag.     'mixing ratio  (0.1 g/kg)')
c
      call zaiopf(flnmfor(1:lgth)//'forcing.precip.a', 'old', 906)
      open (unit=906,file=flnmfor(1:lgth)//'forcing.precip.b',
     &   status='old', action='read')
       read (906,103) preambl
      if     (mnproc.eq.1) then
      write (lp,102) preambl
      call flush(lp)
      endif !1st tile
      pcmax = -huge
      call rdmonth(util1, 906)
      do j=1,jj
        do i=1,ii
          pcmax = max(pcmax,util1(i,j))
        enddo
      enddo
      call xcmaxr(pcmax)
cdiag call prtmsk(ip,util1,util2,idm,idm,jdm,  0.,86400.*36000.,
cdiag.     'precipitation (cm/year) ')
c
c --- zero fields implies no surface salinity flux
      pcipf = pcmax.ne.0.0
      if     (.not.pcipf) then
        if     (mnproc.eq.1) then
        write (lp,*)
        write (lp,*) '***** no surface salinity flux *****'
        write (lp,*)
        call flush(lp)
        endif !1st tile
      endif
c
      call zaiopf(flnmfor(1:lgth)//'forcing.radflx.a', 'old', 907)
      open (unit=907,file=flnmfor(1:lgth)//'forcing.radflx.b',
     &   status='old', action='read')
       read (907,103) preambl
      if     (mnproc.eq.1) then
      write (lp,102) preambl
      call flush(lp)
      endif !1st tile
      call rdmonth(util1, 907)
cdiag call prtmsk(ip,util1,util2,idm,idm,jdm,  0.,1.0,
cdiag.     'net radiation (w/sq_m)  ')
c
      call zaiopf(flnmfor(1:lgth)//'forcing.shwflx.a', 'old', 908)
      open (unit=908,file=flnmfor(1:lgth)//'forcing.shwflx.b',
     &   status='old', action='read')
       read (908,103) preambl
      if     (mnproc.eq.1) then
      write (lp,102) preambl
      call flush(lp)
      endif !1st tile
      call rdmonth(util1, 908)
cdiag call prtmsk(ip,util1,util2,idm,idm,jdm,  0.,1.0,
cdiag.     'sw radiation (w/sq_m)  ')
c
c --- calculate jerlov water type, which governs the penetration depth of
c --- shortwave radiation. water type could be given a seasonal cycle, or
c --- it could be made dependent on variables such as plankton/chloryphyll
c --- concentration or suspended sediment concentration.
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
c ---     map shallow depths to high jerlov numbers
          jerlov(i,j)=6-max(1,min(5,int(depths(i,j)/15.0)))
          jerlov(i,j)=max(jerlv0,jerlov(i,j))
        enddo
      enddo
c
      end if                    !  thermo = .true.
c
      if     (mnproc.eq.1) then
      write (lp,*) ' ...finished opening forcing fields'
      endif !1st tile
      call xcsync(flush_lp)
c
      return
      end
c
c
      subroutine forfunh(dtime)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
      include 'common_blocks.h'
c
      real*8    dtime
c
c --- high frequency atmospheric forcing field processing.
c
c --- units of tau_x  are N/sq_m (positive eastwards)
c --- units of tau_y  are N/sq_m (positive northwards)
c --- units of wndspd are m/s
c --- units of airtmp are degC
c --- units of vapmix are kg/kg
c --- units of precip are m/s    (postive into ocean)
c --- units of radflx are w/sq_m (postive into ocean)
c --- units of swflx  are w/sq_m (postive into ocean)
c
c --- I/O and array I/O units 901-908 are reserved for the entire run.
c
c --- all input fields much be defined at all grid points
c
      real*8    dtime0,dtime1
      save      dtime0,dtime1
c
      character preambl(5)*79,cline*80
      real      pcmax
      integer   i,ios,iunit,j,lgth,nrec
c
 103  format (a79)
 102  format (/(a79))
c
c --- w0 negative on first call only.
      if     (w0.lt.-1.0) then
c
c ---   initialize forcing fields
c
        if      (.not.windf) then
          if     (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error in forfunh - windf must be .true.'
          write(lp,*)
          endif !1st tile
          call xcstop('(forfunh)')
                 stop '(forfunh)'
        elseif (.not.thermo) then
          if     (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error in forfunh - thermo must be .true.'
          write(lp,*)
          endif !1st tile
          call xcstop('(forfunh)')
                 stop '(forfunh)'
        endif
c
c ---   linear interpolation in time, so slots 3 and 4 are zero.
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
              taux(i,j,3) = 0.0
              taux(i,j,4) = 0.0
              tauy(i,j,3) = 0.0
              tauy(i,j,4) = 0.0
            wndspd(i,j,3) = 0.0
            wndspd(i,j,4) = 0.0
            airtmp(i,j,3) = 0.0
            airtmp(i,j,4) = 0.0
            vapmix(i,j,3) = 0.0
            vapmix(i,j,4) = 0.0
            precip(i,j,3) = 0.0
            precip(i,j,4) = 0.0
            radflx(i,j,3) = 0.0
            radflx(i,j,4) = 0.0
             swflx(i,j,3) = 0.0
             swflx(i,j,4) = 0.0
          enddo
        enddo
c
c ---   open all forcing files.
        if     (mnproc.eq.1) then
        write (lp,*) ' now initializing forcing fields ...'
        call flush(lp)
        endif !1st tile
c
        do 12 lgth=28,1,-1
          if (flnmfor(lgth:lgth).ne.' ') then
            goto 112
          endif
  12    continue
        call xcstop('(flnmfor)')
               stop '(flnmfor)'
 112    continue
c
        call zaiopf(flnmfor(1:lgth)//'forcing.tauewd.a', 'old', 901)
        open (unit=901,file=flnmfor(1:lgth)//'forcing.tauewd.b',
     &        status='old', action='read')
         read (901,103) preambl
        if     (mnproc.eq.1) then
        write (lp, 102) preambl
        call flush(lp)
        endif !1st tile
c
        call zaiopf(flnmfor(1:lgth)//'forcing.taunwd.a', 'old', 902)
        open (unit=902,file=flnmfor(1:lgth)//'forcing.taunwd.b',
     &     status='old', action='read')
         read (902,103) preambl
        if     (mnproc.eq.1) then
        write (lp, 102) preambl
        call flush(lp)
        endif !1st tile
c
        call zaiopf(flnmfor(1:lgth)//'forcing.wndspd.a', 'old', 903)
        open (unit=903,file=flnmfor(1:lgth)//'forcing.wndspd.b',
     &     status='old', action='read')
         read (903,103) preambl
        if     (mnproc.eq.1) then
        write (lp, 102) preambl
        call flush(lp)
        endif !1st tile
c
        call zaiopf(flnmfor(1:lgth)//'forcing.airtmp.a', 'old', 904)
        open (unit=904,file=flnmfor(1:lgth)//'forcing.airtmp.b',
     &     status='old', action='read')
         read (904,103) preambl
        if     (mnproc.eq.1) then
        write (lp, 102) preambl
        call flush(lp)
        endif !1st tile
c
        call zaiopf(flnmfor(1:lgth)//'forcing.vapmix.a', 'old', 905)
        open (unit=905,file=flnmfor(1:lgth)//'forcing.vapmix.b',
     &     status='old', action='read')
         read (905,103) preambl
        if     (mnproc.eq.1) then
        write (lp, 102) preambl
        call flush(lp)
        endif !1st tile
c
        call zaiopf(flnmfor(1:lgth)//'forcing.precip.a', 'old', 906)
        open (unit=906,file=flnmfor(1:lgth)//'forcing.precip.b',
     &     status='old', action='read')
         read (906,103) preambl
        if     (mnproc.eq.1) then
        write (lp, 102) preambl
        call flush(lp)
        endif !1st tile
c
        call zaiopf(flnmfor(1:lgth)//'forcing.radflx.a', 'old', 907)
        open (unit=907,file=flnmfor(1:lgth)//'forcing.radflx.b',
     &     status='old', action='read')
         read (907,103) preambl
        if     (mnproc.eq.1) then
        write (lp, 102) preambl
        call flush(lp)
        endif !1st tile
c
        call zaiopf(flnmfor(1:lgth)//'forcing.shwflx.a', 'old', 908)
        open (unit=908,file=flnmfor(1:lgth)//'forcing.shwflx.b',
     &     status='old', action='read')
         read (908,103) preambl
        if     (mnproc.eq.1) then
        write (lp, 102) preambl
        call flush(lp)
        endif !1st tile
c
c ---   skip ahead to the start time.
        dtime1 = huge
        do 11 nrec=1,9999999
          dtime0 = dtime1
          read(901,'(a)',iostat=ios) cline
          if     (ios.ne.0) then
            if     (mnproc.eq.1) then
              write(lp,*)
              write(lp,*) 'error in forfunh - hit end of input'
              write(lp,*) 'dtime0,dtime1 = ',dtime0,dtime1
              write(lp,*) 'dtime = ',dtime
              write(lp,*)
            endif !1st tile
            call xcstop('(forfunh)')
                   stop '(forfunh)'
          endif
          i = index(cline,'=')
          read (cline(i+1:),*) dtime1
          if     (yrflag.eq.2) then
            if     (nrec.eq.1 .and. abs(dtime1-1096.0d0).gt.0.01) then
c
c ---         climatology must start on wind day 1096.0, 01/01/1904.
              if     (mnproc.eq.1) then
              write(lp,*)
              write(lp,*) 'error in forfunh - forcing climatology',
     &                    ' must start on wind day 1096'
              write(lp,*) 'dtime1 = ',dtime1
              write(lp,*)
              endif !1st tile
              call xcstop('(forfunh)')
                     stop '(forfunh)'
            endif
            dtime1 = (dtime1 - 1096.0d0) + 
     &               366.0d0*int((dtime+0.00001d0)/366.0d0)
            if     (nrec.ne.1 .and. dtime1.lt.dtime0) then
              dtime1 = dtime1 + 366.0d0
            endif
          endif
          if     (dtime0.le.dtime .and. dtime1.gt.dtime) then
            goto 111
          endif
   11   continue
  111   continue
        rewind(unit=901)
        read (901,103) preambl
c
        do iunit= 901,908
          do i= 1,nrec-2
            read (iunit,*)
            call zaiosk(iunit)
          enddo
        enddo
        pcipf  = .true.
        dtime1 = huge
        call rdpall(dtime0,dtime1)
        if     (yrflag.eq.2) then
          dtime1 = (dtime1 - 1096.0d0) + 
     &             366.0d0*int((dtime+0.00001d0)/366.0d0)
        endif
        call rdpall(dtime0,dtime1)
        if     (yrflag.eq.2) then
          dtime1 = (dtime1 - 1096.0d0) + 
     &             366.0d0*int((dtime+0.00001d0)/366.0d0)
          if     (dtime1.lt.dtime0) then
            dtime1 = dtime1 + 366.0d0
          endif
        endif
c
c ---   zero precip field implies no surface salinity flux
        pcmax = -huge
        do j=1,jj
          do i=1,ii
            pcmax = max(pcmax,precip(i,j,1),precip(i,j,2))
          enddo
        enddo
        call xcmaxr(pcmax)
        pcipf = pcmax.ne.0.0
        if     (.not.pcipf) then
          if     (mnproc.eq.1) then
          write (lp,*)
          write (lp,*) '***** no surface salinity flux *****'
          write (lp,*)
          call flush(lp)
          endif !1st tile
        endif
c
c ---   calculate jerlov water type, which governs the penetration depth of
c ---   shortwave radiation. water type could be given a seasonal cycle, or
c ---   it could be made dependent on variables such as plankton/chloryphyll
c ---   concentration or suspended sediment concentration.
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
c ---       map shallow depths to high jerlov numbers
            jerlov(i,j)=6-max(1,min(5,int(depths(i,j)/15.0)))
            jerlov(i,j)=max(jerlv0,jerlov(i,j))
          enddo
        enddo
        if     (mnproc.eq.1) then
        write (lp,*) 
        write (lp,*) ' dtime,dtime0,dtime1 = ',dtime,dtime0,dtime1
        write (lp,*) 
        write (lp,*) ' ...finished initializing forcing fields'
        endif !1st tile
        call xcsync(flush_lp)
      endif  ! initialization
c
      if     (dtime.gt.dtime1) then
c
c ---   get the next set of fields.
*           if     (mnproc.eq.1) then
*           write(lp,*) 'enter rdpall - ',dtime,dtime0,dtime1
*           call flush(lp)
*           endif !1st tile
        call rdpall(dtime0,dtime1)
        if     (yrflag.eq.2) then
          dtime1 = (dtime1 - 1096.0d0) + 
     &             366.0d0*int((dtime+0.00001d0)/366.0d0)
          if     (dtime1.lt.dtime0) then
            dtime1 = dtime1 + 366.0d0
          endif
        endif
*           if     (mnproc.eq.1) then
*           write(lp,*) ' exit rdpall - ',dtime,dtime0,dtime1
*           call flush(lp)
*           endif !1st tile
      endif
c
c --- linear interpolation in time.
      w0 = (dtime1-dtime)/(dtime1-dtime0)
      w1 = 1.0 - w0
*           if     (mnproc.eq.1) then
*           write(lp,*) 'rdpall - dtime,w0,w1 = ',dtime,w0,w1
*           call flush(lp)
*           endif !1st tile
      return
      end
c
c
      subroutine forfunr
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
      include 'common_blocks.h'
c
c --- initialize input of relaxation forcing fields
c
c --- rmu   is a single field specifying 1 / e-folding time (1/s)
c ---        set to zero where there is no boundary relaxation
c --- twall is temperature climatology for all layers
c --- swall is salinity    climatology for all layers
c --- pwall is interface   climatology for all layers (pressure units)
c
c --- I/O array I/O units 910-913 are reserved for the entire run.
c
c --- all input fields much be defined at all grid points
c
      integer         mreca,mrecr
      common/rdforfi/ mreca,mrecr
      save  /rdforfi/
c
      integer   lgth,mo
      character preambl(5)*79
      real      one,oneps
      integer   i,incmon,j,k,nrec
c
 103  format (a79)
 102  format (/(a79))
c
      if (.not.relaxf) then
        return
      endif
c
c --- read fields needed for boundary relaxation
c
      mrecr=1
      if     (mnproc.eq.1) then
      write (lp,*) ' now opening relaxation fields ...'
      call flush(lp)
      endif !1st tile
c
      do 13 lgth=28,1,-1
        if (flnmforw(lgth:lgth).ne.' ') then
          goto 113
        endif
  13  continue
      call xcstop('(flnmforw)')
             stop '(flnmforw)'
 113  continue
c
      call zaiopf(flnmforw(1:lgth)//'relax.rmu.a', 'old', 910)
      open (unit=910,file=flnmforw(1:lgth)//'relax.rmu.b',
     &      status='old', action='read')
       read (910,103) preambl
      if     (mnproc.eq.1) then
      write (lp, 102) preambl
      call flush(lp)
      endif !1st tile
      call rdmonth(rmu, 910)
      close (unit=910)
      call zaiocl(910)
c
      call zaiopf(flnmforw(1:lgth)//'relax.temp.a', 'old', 911)
      open (unit=911,file=flnmforw(1:lgth)//'relax.temp.b',
     &      status='old', action='read')
       read (911,103) preambl
      if     (mnproc.eq.1) then
      write (lp, 102) preambl
      call flush(lp)
      endif !1st tile
      do k=1,kk
        call rdmonth(util1, 911)
      enddo
c
      call zaiopf(flnmforw(1:lgth)//'relax.saln.a', 'old', 912)
      open (unit=912,file=flnmforw(1:lgth)//'relax.saln.b',
     &      status='old', action='read')
       read (912,103) preambl
      if     (mnproc.eq.1) then
      write (lp, 102) preambl
      call flush(lp)
      endif !1st tile
      do k=1,kk
        call rdmonth(util1, 912)
      enddo
c
      call zaiopf(flnmforw(1:lgth)//'relax.intf.a', 'old', 913)
      open (unit=913,file=flnmforw(1:lgth)//'relax.intf.b',
     &      status='old', action='read')
       read (913,103) preambl
      if     (mnproc.eq.1) then
      write (lp, 102) preambl
      call flush(lp)
      endif !1st tile
      do k=1,kk
        call rdmonth(util1, 913)
      enddo
c
      if     (mnproc.eq.1) then
      write (lp,*) ' ...finished opening relaxation fields'
      endif !1st tile
      call xcsync(flush_lp)
c
      return
      end
c
c
      subroutine rdmonth(field, iunit)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
      include 'common_blocks.h'
c
      integer   iunit
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &          field
c
c --- read a single array field from unit iunit.
c
c --- iunit=901-908; atmospheric forcing field
c --- iunit=910;     relaxation time scale field
c --- iunit=911-913; relaxation forcing field
c
      integer   i,layer,mnth
      real      denlay,hmina,hminb,hmaxa,hmaxb
      character cline*80
c
      read ( iunit,'(a)')  cline
      if     (mnproc.eq.1) then
      write (lp,   '(a)')  cline
      endif !1st tile
      i = index(cline,'=')
      if     (iunit.ge.901 .and. iunit.le.908) then
c ---   atmospheric forcing
        read (cline(i+1:),*) mnth,hminb,hmaxb
      elseif (iunit.ge.911 .and. iunit.le.913) then
c ---   relaxation forcing
        read (cline(i+1:),*) mnth,layer,denlay,hminb,hmaxb
      elseif (iunit.eq.910) then
c ---   relaxation time scale
        read (cline(i+1:),*) hminb,hmaxb
      else
        if     (mnproc.eq.1) then
        write(lp,'(a / a,i5)')
     &    'error - iunit must be 901-908 or 910-913',
     &    'iunit =',iunit
        endif !1st tile
        call xcstop('(rdmonth)')
               stop '(rdmonth)'
      endif
c
      call zaiord(field,ip,.false., hmina,hmaxa, iunit)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,i3 / a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    'iunit = ',iunit,
     &    cline,
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        endif !1st tile
        call xcstop('(rdmonth)')
               stop '(rdmonth)'
      endif
c
      return
      end
c
c
      subroutine rdpall(dtime0,dtime1)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      real*8  dtime0,dtime1
c
c --- copy slot 2 into slot 1, and 
c --- read a set of high frequency forcing fields into slot 2.
c --- on exit, dtime0 and dtime1 are the associated times (wind days).
c
      integer k
      real*8  dtime(901:908)
c
      call rdpall1(  taux,dtime(901),901)
      call rdpall1(  tauy,dtime(902),902)
      call rdpall1(wndspd,dtime(903),903)
      call rdpall1(airtmp,dtime(904),904)
      call rdpall1(vapmix,dtime(905),905)
      if     (pcipf) then
        call rdpall1(precip,dtime(906),906)
      else
        dtime(906) = dtime(905)
      endif
      call rdpall1(radflx,dtime(907),907)
      call rdpall1( swflx,dtime(908),908)
c
      dtime0 = dtime1
      dtime1 = dtime(901)
c
c --- check the input times.
      do k= 902,908
        if     (dtime(k).ne.dtime1) then
          if     (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error in rdpall - inconsistent forcing times'
          write(lp,*) 'dtime0,dtime1 = ',dtime0,dtime1
          write(lp,*) 'dtime = ',dtime
          write(lp,*)
          endif !1st tile
          call xcstop('(rdpall)')
                 stop '(rdpall)'
        endif
      enddo
      return
      end
c
c
      subroutine rdpall1(field,dtime,iunit)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
      include 'common_blocks.h'
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     &        field
      real*8  dtime
      integer iunit
c
c --- copy field(:,:,2) into field(:,:,1), and
c --- read a high frequency forcing field into field(:,:,2).
c --- on exit, dtime is the time (wind day) of the forcing field.
c
      integer   i,j
      character cline*80
      real      hmina,hminb,hmaxa,hmaxb,span
c
      read ( iunit,'(a)',end=100) cline
      if     (mnproc.eq.1) then
      write (lp,   '(a)')         cline
      endif !1st tile
      goto 110
  100 continue
        if     (yrflag.eq.2) then
*         if     (mnproc.eq.1) then
*         write(lp,*) 'rdpall1 - rewind unit ',iunit
*         call flush(lp)
*         endif !1st tile
          rewind iunit
          read (iunit,*)
          read (iunit,*)
          read (iunit,*)
          read (iunit,*)
          read (iunit,*)
          read (iunit,'(a)') cline
          if     (mnproc.eq.1) then
          write (lp,  '(a)') cline
          endif !1st tile
          call zaiorw(iunit)
        else
          if     (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error in rdpall1 - end of file'
          write(lp,*) 'iunit = ',iunit
          write(lp,*)
          endif !1st tile
          call xcstop('(rdpall)')
                 stop '(rdpall)'
        endif
  110 continue
c
      i = index(cline,'=')
      read (cline(i+1:),*) dtime,span,hminb,hmaxb
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1-nbdy,jj+nbdy
        do i= 1-nbdy,ii+nbdy
          field(i,j,1) = field(i,j,2)
        enddo
      enddo
c
      call zaiord(field(1-nbdy,1-nbdy,2),ip,.false., hmina,hmaxa, iunit)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,i3 / a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    'iunit = ',iunit,
     &    cline,
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        endif !1st tile
        call xcstop('(rdpall1)')
               stop '(rdpall1)'
      endif
c
c --- wind stress uses the the halo.
c
      if     (iunit.eq.901) then  ! taux
        call xctilr(field(1-nbdy,1-nbdy,2),1,1, nbdy,nbdy, halo_uv)
      elseif (iunit.eq.902) then  ! tauy
        call xctilr(field(1-nbdy,1-nbdy,2),1,1, nbdy,nbdy, halo_vv)
      endif
      return
      end
c
c
      subroutine rdforf(mnth,lslot)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
      include 'common_blocks.h'
c
      integer lslot,mnth
c
c --- read forcing functions for one month.
c
      integer         mreca,mrecr
      common/rdforfi/ mreca,mrecr
      save  /rdforfi/
c
      integer i,irec,iunit,j
c
      if     (mnth.le.mreca) then
c
c ---   rewind all units
c
        if     (windf) then
          do iunit= 901,902
            rewind iunit
            read  (iunit,*)
            read  (iunit,*)
            read  (iunit,*)
            read  (iunit,*)
            read  (iunit,*)
            call zaiorw(iunit)
          enddo
        endif
        if     (thermo) then
          do iunit= 903,908
            rewind iunit
            read  (iunit,*)
            read  (iunit,*)
            read  (iunit,*)
            read  (iunit,*)
            read  (iunit,*)
            call zaiorw(iunit)
          enddo
        endif
        if     (mnproc.eq.1) then
        write(lp,*) 'rdforf: mreca,mnth = ',mreca,mnth,'  (rewind)'
        call flush(lp)
        endif !1st tile
        mreca = 0
      endif
c
c --- skip forward to desired month
c
      do irec= mreca+1,mnth-1
        if     (mnproc.eq.1) then
        write(lp,*) 'rdforf: mreca,mnth = ',mreca,mnth,
     &              '  (skipping ',irec,')'
        call flush(lp)
        endif !1st tile
        if     (windf) then
          do iunit= 901,902
            read (iunit,*)
            call zaiosk(iunit)
          enddo
        endif
        if     (thermo) then
          do iunit= 903,908
            read (iunit,*)
            call zaiosk(iunit)
          enddo
        endif
      enddo
c
c --- read desired month
c
      if     (windf) then
        call rdmonth(taux(1-nbdy,1-nbdy,lslot),901)
        call xctilr( taux(1-nbdy,1-nbdy,lslot),1,1, nbdy,nbdy, halo_uv)
        call rdmonth(tauy(1-nbdy,1-nbdy,lslot),902)
        call xctilr( tauy(1-nbdy,1-nbdy,lslot),1,1, nbdy,nbdy, halo_vv)
      else
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            taux(i,j,lslot) = 0.0
            tauy(i,j,lslot) = 0.0
          enddo
        enddo
      endif
      if     (thermo) then
        call rdmonth(wndspd(1-nbdy,1-nbdy,lslot),903)
        call rdmonth(airtmp(1-nbdy,1-nbdy,lslot),904)
        call rdmonth(vapmix(1-nbdy,1-nbdy,lslot),905)
        call rdmonth(precip(1-nbdy,1-nbdy,lslot),906)
        call rdmonth(radflx(1-nbdy,1-nbdy,lslot),907)
        call rdmonth( swflx(1-nbdy,1-nbdy,lslot),908)
      else
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            wndspd(i,j,lslot) = 0.0
            airtmp(i,j,lslot) = 0.0
            vapmix(i,j,lslot) = 0.0
            precip(i,j,lslot) = 0.0
            radflx(i,j,lslot) = 0.0
             swflx(i,j,lslot) = 0.0
          enddo
        enddo
      endif
c
      mreca = mnth
c
      if     (mnproc.eq.1) then
      write (lp,'(2(a,i3))') ' forcing functions for month',mnth,
     &   ' written into slot',lslot
      call flush(lp)
      endif !1st tile
      return
      end
c
c
      subroutine rdrlax(month,lslot)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      include 'common_blocks.h'
c
      integer lslot,month
c
c --- read relaxation fields for one month,
c --- monthly (clmflg==12) or bi-monthly (clmflg==6) data.
c
      integer         mreca,mrecr
      common/rdforfi/ mreca,mrecr
      save  /rdforfi/
c
      integer irec,iunit,k,mrec,mnth
c
      mnth=mod(month-1,12)+1
*     if     (mnproc.eq.1) then
*     write(lp,*) 'rdrlax - month = ',month,mnth
*     call flush(lp)
*     endif !1st tile
c
      if     (relaxf) then
        if     (clmflg.eq.12) then
          mrec = mnth
        else
          mrec = (mnth+1)/2
        endif
c
        if     (mrec.le.mrecr) then
c
c ---     rewind all units
c
          do iunit= 911,913
            rewind iunit
            read  (iunit,*)
            read  (iunit,*)
            read  (iunit,*)
            read  (iunit,*)
            read  (iunit,*)
            call zaiorw(iunit)
          enddo
          if     (mnproc.eq.1) then
          write(lp,*) 'rdrlax: mrecr,mrec = ',mrecr,mrec,'  (rewind)'
          call flush(lp)
          endif !1st tile
          mrecr = 0
        endif
c
c ---   skip forward to desired month
c
        do irec= mrecr+1,mrec-1
          if     (mnproc.eq.1) then
          write(lp,*) 'rdrlax: mrecr,mrec = ',mrecr,mrec,
     &                '  (skipping ',irec,')'
          call flush(lp)
          endif !1st tile
          do iunit= 911,913
            do k= 1,kk
              read (iunit,*)
              call zaiosk(iunit)
            enddo
          enddo
        enddo
c
c --- read desired month
c
        do k= 1,kk
          call rdmonth(twall(1-nbdy,1-nbdy,k,lslot),911)
        enddo
        do k= 1,kk
          call rdmonth(swall(1-nbdy,1-nbdy,k,lslot),912)
        enddo
        do k= 1,kk
          call rdmonth(pwall(1-nbdy,1-nbdy,k,lslot),913)
        enddo
c
        mrecr = mrec
c
        if     (mnproc.eq.1) then
        write (lp,'(2(a,i3))') ' relaxation fields for month',mnth,
     &     ' written into slot',lslot
        call flush(lp)
        endif !1st tile
      endif
      return
      end

c
c
c> Revision history:
c>
c> Mar. 1995 - added logical variable 'windf'
c> Oct. 1997 - made necessary changes to reduce time dimension from 12 to 4
c> Oct. 1999 - added code to read and store shortwave heat flux used for
c>             penetrating shortwave radiation
c> Jan. 2000 - removed all conversion factors (apply before input)
c> Jan. 2000 - removed biasrd and biaspc      (apply before input)
c> May  2000 - conversion to SI units, positive flux into ocean
c> Aug  2000 - added option for high frequency atmospheric forcing
c> Jan  2001 - Converted from pakk to array input file type
