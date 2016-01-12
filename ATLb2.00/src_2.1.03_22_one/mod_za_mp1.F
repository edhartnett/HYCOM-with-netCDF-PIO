c
c-----------------------------------------------------------------------
c
c     machine dependent I/O routines.
c     message passing version, with all I/O from first processor.
c     contained in module mod_za.
c
c     author:  Alan J. Wallcraft,  NRL.
c
c-----------------------------------------------------------------------
c
      subroutine zagetc(cline,ios, iunit)
      implicit none
c
      character*80, intent(out)   :: cline
      integer,      intent(out)   :: ios
      integer,      intent(in)    :: iunit
c
c**********
c*
c  1) machine specific routine for reading one text line from a file.
c
c  2) The read is performed on the first processor only.
c*
c**********
c
      integer        iline,ibuf
      common/czgetc/ iline(81,0:1),ibuf
      save  /czgetc/
c
      integer i
c
c --- I/O from first processor only
c
      ibuf = mod(ibuf+1,2)
c
      if     (mnproc.eq.1) then
        read(iunit,'(a)',iostat=ios) cline
        do i= 1,80
          iline(i,ibuf) = ichar(cline(i:i))
        enddo
        iline(81,ibuf) = ios
      endif
c
c     broadcast to all other processors
c
      call xcgetc(iline(:,ibuf))
      do i= 1,80
        cline(i:i) = char(iline(i,ibuf))
      enddo
      ios = iline(81,ibuf)  ! iostat value
      return
      end subroutine zagetc

      subroutine zaiopn(cstat, iaunit)
      implicit none
c
      integer,       intent(in)    :: iaunit
      character*(*), intent(in)    :: cstat
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for opening a file for array i/o.
c
c     must call zaiost before first call to zaiopn.
c     see also 'zaiope' and 'zaiopf'.
c
c  2) the filename is taken from the environment variable FORxxxA,
c       where xxx = iaunit, with default fort.xxxa.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     cstat indicates the file type, it can be 'scratch', 'old', or
c      'new'.
c     all i/o to iaunit must be performed by zaiord and zaiowr.
c     the file should be closed using zaiocl.
c*
c**********
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
c --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
      real*4         w
      common/czioxw/ w(n2drec)
      save  /czioxw/
c
      integer   ios,nrecl
      character cfile*256,cenv*7
      character cact*9
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
c     test file state.
c
      if     (iarec(iaunit).ne.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiopn)')
               stop '(zaiopn)'
      endif
c
      iarec(iaunit) = 0
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiopn - iaunit = ',iaunit
*     call flush(lp)
c
c     get filename.
c
      write(cenv,1000) iaunit
      cfile = ' '
      call getenv(cenv,cfile)
      if     (cfile.eq.' ') then
        write(cfile,1100) iaunit
      endif
c
c     open file.
c
      inquire(iolength=nrecl) w
c
      if     (cstat.eq.'OLD' .or.
     &        cstat.eq.'old'     ) then
        cact = 'READ'
      elseif (cstat.eq.'NEW' .or.
     &        cstat.eq.'new'     ) then
        cact = 'WRITE'
      else
        cact = 'READWRITE'
      endif
#if defined(YMP)
      if     (mod(nrecl,16384).eq.0 .and. nrecl.gt.16384*4) then
       call asnunit(iaunit+1000,'-F syscall -N ieee',ios)
      else
        call asnunit(iaunit+1000,'-F cachea:8:16:2 -N ieee',ios)
      endif
      if     (ios.ne.0) then
        write(lp,9050) iaunit
        write(lp,*) 'ios = ',ios
        call flush(lp)
        call xchalt('(zaiopn)')
               stop '(zaiopn)'
      endif
#endif
      if     (cstat.eq.'scratch' .or.
     &        cstat.eq.'SCRATCH'     ) then
        open(unit=iaunit+1000,             
     &       form='unformatted', status='scratch',
     &       access='direct', recl=nrecl, action=cact, iostat=ios)
      else
        open(unit=iaunit+1000, file=cfile, 
     &       form='unformatted', status=cstat,
     &       access='direct', recl=nrecl, action=cact, iostat=ios)
      endif
      if     (ios.ne.0) then
        write(lp,9100) iaunit
        write(lp,*) 'ios = ',ios
        call flush(lp)
        call xchalt('(zaiopn)')
               stop '(zaiopn)'
      endif
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 1000 format('FOR',i3.3,'A')
 1100 format('fort.',i3.3,'a')
 9000 format(/ /10x,'error in zaiopn -  array I/O unit ',
     &   i3,' is not marked as available.'/ /)
#if defined(YMP)
 9050 format(/ /10x,'error in zaiopn -  can''t asnunit ',i3,
     &   ', for array I/O.'/ /)
#endif
 9100 format(/ /10x,'error in zaiopn -  can''t open unit ',i3,
     &   ', for array I/O.'/ /)
      end subroutine zaiopn

      subroutine zaiope(cenv,cstat, iaunit)
      implicit none
c
      integer,       intent(in)    :: iaunit
      character*(*), intent(in)    :: cenv,cstat
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for opening a file for array i/o.
c
c     must call zaiost before first call to zaiope.
c     see also 'zaiopn' and 'zaiopf'.
c
c  2) the filename is taken from environment variable 'cenv'.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     cstat indicates the file type, it can be 'scratch', 'old', or
c      'new'.
c     all i/o to iaunit must be performed by zaiord and zaiowr.
c      arrays passed to these routines must conform to 'h'.
c     the file should be closed using zaiocl.
c*
c**********
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
c --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
      real*4         w
      common/czioxw/ w(n2drec)
      save  /czioxw/
c
      integer   len_trim
      integer   ios,nrecl
      character cfile*256
      character cact*9
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
c     test file state.
c
      if     (iarec(iaunit).ne.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiope)')
               stop '(zaiope)'
      endif
c
      iarec(iaunit) = 0
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiope - iaunit = ',iaunit
*     call flush(lp)
c
c     get filename.
c
      cfile = ' '
      call getenv(cenv,cfile)
      if     (cfile.eq.' ') then
        write(lp,9300) cenv(1:len_trim(cenv))
        write(lp,*) 'iaunit = ',iaunit
        call flush(lp)
        call xchalt('(zaiope)')
               stop '(zaiope)'
      endif
c
c     open file.
c
      inquire(iolength=nrecl) w
c
      if     (cstat.eq.'OLD' .or.
     &        cstat.eq.'old'     ) then
        cact = 'READ'
      elseif (cstat.eq.'NEW' .or.
     &        cstat.eq.'new'     ) then
        cact = 'WRITE'
      else
        cact = 'READWRITE'
      endif
c
#if defined(YMP)
      if     (mod(nrecl,16384).eq.0 .and. nrecl.gt.16384*4) then
       call asnunit(iaunit+1000,'-F syscall -N ieee',ios)
      else
        call asnunit(iaunit+1000,'-F cachea:8:16:2 -N ieee',ios)
      endif
      if     (ios.ne.0) then
        write(lp,9050) iaunit,cfile(1:len_trim(cfile))
        write(lp,*) 'ios = ',ios
        write(lp,*) 'cenv = ',cenv(1:len_trim(cenv))
        call flush(lp)
        call xchalt('(zaiope)')
               stop '(zaiope)'
      endif
#endif
      open(unit=iaunit+1000, file=cfile, 
     &     form='unformatted', status=cstat,
     &     access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(lp,9100) iaunit,cfile(1:len_trim(cfile))
        write(lp,*) 'ios  = ',ios
        write(lp,*) 'cenv = ',cenv(1:len_trim(cenv))
        call flush(lp)
        call xchalt('(zaiope)')
               stop '(zaiope)'
      endif
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiope -  array I/O unit ',
     &   i3,' is not marked as available.'/ /)
#if defined(YMP)
 9050 format(/ /10x,'error in zaiope -  can''t asnunit ',i3,
     &   ', for array I/O.' /
     &   10x,'cfile = ',a/ /)
#endif
 9100 format(/ /10x,'error in zaiope -  can''t open unit ',i3,
     &   ', for array I/O.' /
     &   10x,'cfile = ',a/ /)
 9300 format(/ /10x,'error in zaiope -  environment variable ',a,
     &   ' not defined'/ /)
      end subroutine zaiope

      subroutine zaiopf(cfile,cstat, iaunit)
      implicit none
c
      integer,       intent(in)    :: iaunit
      character*(*), intent(in)    :: cfile,cstat
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for opening a file for array i/o.
c
c     must call zaiost before first call to zaiopf.
c     see also 'zaiopn' and 'zaiope'.
c
c  2) the filename is taken from 'cfile'.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     cstat indicates the file type, it can be 'scratch', 'old', or
c      'new'.
c     all i/o to iaunit must be performed by zaiord and zaiowr.
c      arrays passed to these routines must conform to 'h'.
c     the file should be closed using zaiocl.
c*
c**********
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
c --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
      real*4         w
      common/czioxw/ w(n2drec)
      save  /czioxw/
c
      integer   len_trim
      integer   ios,nrecl
      character cact*9
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
c     test file state.
c
      if     (iarec(iaunit).ne.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiopf)')
               stop '(zaiopf)'
      endif
c
      iarec(iaunit) = 0
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiopf - iaunit = ',iaunit
*     call flush(lp)
c
c     open file.
c
      inquire(iolength=nrecl) w
c
      if     (cstat.eq.'OLD' .or.
     &        cstat.eq.'old'     ) then
        cact = 'READ'
      elseif (cstat.eq.'NEW' .or.
     &        cstat.eq.'new'     ) then
        cact = 'WRITE'
      else
        cact = 'READWRITE'
      endif
c
#if defined(YMP)
      if     (mod(nrecl,16384).eq.0 .and. nrecl.gt.16384*4) then
       call asnunit(iaunit+1000,'-F syscall -N ieee',ios)
      else
        call asnunit(iaunit+1000,'-F cachea:8:16:2 -N ieee',ios)
      endif
      if     (ios.ne.0) then
        write(lp,9050) iaunit,cfile(1:len_trim(cfile))
        write(lp,*) 'ios   = ',ios
        call flush(lp)
        call xchalt('(zaiopf)')
               stop '(zaiopf)'
      endif
#endif
      open(unit=iaunit+1000, file=cfile, 
     &     form='unformatted', status=cstat,
     &     access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(lp,9100) iaunit,cfile(1:len_trim(cfile))
        write(lp,*) 'ios  = ',ios
        call flush(lp)
        call xchalt('(zaiopf)')
               stop '(zaiopf)'
      endif
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiopf -  array I/O unit ',
     &   i3,' is not marked as available.'/ /)
#if defined(YMP)
 9050 format(/ /10x,'error in zaiopf -  can''t asnunit ',i3,
     &   ', for array I/O.' /
     &   10x,'cfile = ',a/ /)
#endif
 9100 format(/ /10x,'error in zaiopf -  can''t open unit ',i3,
     &   ', for array I/O.' /
     &   10x,'cfile = ',a/ /)
      end subroutine zaiopf

      subroutine zaiopi(lopen, iaunit)
      implicit none
c
      logical, intent(out)   :: lopen
      integer, intent(in)    :: iaunit
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) is an array i/o unit open?
c
c  2) must call zaiost before first call to zaiopi.
c*
c**********
c
      lopen = iarec(iaunit).ne.-1
      return
      end subroutine zaiopi

      subroutine zaiost
      implicit none
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for initializing array i/o.
c
c  2) see also zaiopn, zaiord, zaiowr, and zaiocl.
c*
c**********
c
      integer i
c
      do i= 1,999
        iarec(i) = -1
      enddo
#if defined(TIMER)
c
c     initialize timers.
c
      call xctmrn(16,'zaio**')
      call xctmrn(17,'zaiord')
      call xctmrn(18,'zaiowr')
#endif
      return
      end subroutine zaiost

      subroutine zaiocl(iaunit)
      implicit none
c
      integer, intent(in)    :: iaunit
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array i/o file closing.
c
c     must call zaiopn for this array unit before calling zaiocl.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c*
c**********
c
      integer ios
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
      if     (iarec(iaunit).lt.0) then
        write(lp,9000) iaunit
        call xcstop('(zaiocl)')
               stop '(zaiocl)'
      endif
c
      iarec(iaunit) = -1
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiocl - iaunit = ',iaunit
*     call flush(lp)
c
      close(unit=iaunit+1000, status='keep')
#if defined(T3E) || defined(YMP)
      call asnunit(iaunit+1000,'-R',ios)
#endif
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiocl -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
      end subroutine zaiocl

      subroutine zaiofl(iaunit)
      implicit none
c
      integer, intent(in)    :: iaunit
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array i/o buffer flushing.
c
c     must call zaiopn for this array unit before calling zaiocl.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c*
c**********
c
      integer   irlen
      character cfile*256
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
      if     (iarec(iaunit).lt.0) then
        write(lp,9000) iaunit
        call xcstop('(zaiofl)')
               stop '(zaiofl)'
      endif
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
c
      inquire(unit=iaunit+1000, name=cfile, recl=irlen)
      close(  unit=iaunit+1000, status='keep')
      open(   unit=iaunit+1000, file=cfile, form='unformatted', 
     &        access='direct', recl=irlen)
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiofl -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
      end subroutine zaiofl

      subroutine zaioiq(iaunit, irec)
      implicit none
c
      integer, intent(in)    :: iaunit
      integer, intent(out)   :: irec
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array i/o inquiry.
c
c  2) returns the number of records processed, or -1 for a closed file.
c*
c**********
c
      irec = iarec(iaunit)
      return
      end subroutine zaioiq

      subroutine zaiorw(iaunit)
      implicit none
c
      integer, intent(in)    :: iaunit
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array i/o file rewinding.
c
c     must call zaiopn for this array unit before calling zaiocl.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c*
c**********
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
      if     (iarec(iaunit).lt.0) then
        write(lp,9000) iaunit
        call xcstop('(zaiorw)')
               stop '(zaiorw)'
      endif
c
      iarec(iaunit) = 0
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiorw - iaunit,rec = ',iaunit,iarec(iaunit)
*     call flush(lp)
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiorw -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
      end subroutine zaiorw

      subroutine zaiord3(h, l, mask,lmask, hmin,hmax,  iaunit)
      implicit none
c
      logical, intent(in)    :: lmask
      integer, intent(in)    :: l,iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin(l),hmax(l)
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l),
     &         intent(out)   :: h
#else
      real,    intent(out)   :: hmin(l),hmax(l)
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l),
     &         intent(out)   :: h
#endif
c
c**********
c*
c  1) machine specific routine for 3-d array reading.
c
c     must call zaiopn for this array unit before calling zaiord.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zaiopn.
c
c  4) hmin,hmax are returned as the minimum and maximum value in the 
c     array, ignoring array elements set to 2.0**100.  
c     if lmask==.true. the range is calculated only where mask.ne.0,
c     with all other values unchanged in h on exit.  It is then an
c     error if mask.ne.0 anywhere the input is 2.0**100.
c*
c**********
c
c     this version just calls zaiord l times.
c
      integer k
c
      do k= 1,l
        call zaiord(h(1-nbdy,1-nbdy,k), mask,lmask,
     &              hmin(k),hmax(k), iaunit)
      enddo
c
      return
      end subroutine zaiord3

      subroutine zaiord(h, mask,lmask, hmin,hmax,  iaunit)
      implicit none
c
      logical, intent(in)    :: lmask
      integer, intent(in)    :: iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin,hmax
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(out)   :: h
#else
      real,    intent(out)   :: hmin,hmax
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(out)   :: h
#endif
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array reading.
c
c     must call zaiopn for this array unit before calling zaiord.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zaiopn.
c
c  4) hmin,hmax are returned as the minimum and maximum value in the 
c     array, ignoring array elements set to 2.0**100.  
c     if lmask==.true. the range is calculated only where mask.ne.0,
c     with all other values unchanged in h on exit.  It is then an
c     error if mask.ne.0 anywhere the input is 2.0**100.
c*
c**********
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
c --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
      real*4         w
      common/czioxw/ w(n2drec)
      save  /czioxw/
c
      real           atmp,amsk,tmsk
      common/czioxr/ atmp(itdm,jtdm),amsk(itdm,jtdm),
     &               tmsk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      save  /czioxr/
c
      integer   ios, i,j
      real*4    wmin,wmax
      real      rmin(1),rmax(1),vsave
#if defined(TIMER)
c
      call xctmr0(17)
#endif
c
      if     (iarec(iaunit).lt.0) then
        write(lp,9000) iaunit
        call xcstop('(zaiord)')
               stop '(zaiord)'
      endif
c
      iarec(iaunit) = iarec(iaunit) + 1
      rmin(1) =  spval 
      rmax(1) = -spval 
      if     (lmask) then
        tmsk(1:ii,1:jj) = mask(1:ii,1:jj)
        vsave = vland
        vland = 0.0
        call xcaget(amsk,tmsk, 1)
        vland = spval
        call xcaget(atmp,h,    1)  ! to preserve masked values
        vland = vsave
      endif
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiord - iaunit,rec = ',iaunit,iarec(iaunit)
*     write(lp,*) 'zaiord - mask.1,1    = ',amsk(1,1)
*     write(lp,*) 'zaiord - h.1,1       = ',atmp(1,1)
*     call flush(lp)
c
      call zaiordd(w,n2drec, iaunit+1000,iarec(iaunit),ios)
      if     (ios.ne.0) then
        write(lp,9100) iarec(iaunit),iaunit
        write(lp,*) 'ios = ',ios
        call flush(lp)
        call xchalt('(zaiord)')
               stop '(zaiord)'
      endif
      wmin =  spval 
      wmax = -spval 
      if     (lmask) then
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&              REDUCTION(MIN:wmin) REDUCTION(MAX:wmax)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j= 1,jtdm
          do i= 1,itdm
            if     (amsk(i,j).ne.0.0) then
              atmp(i,j) = w(i+(j-1)*itdm)
              wmin = min( wmin, w(i+(j-1)*itdm) )
              wmax = max( wmax, w(i+(j-1)*itdm) )
            endif
          enddo
        enddo
        if     (wmax.eq.spval) then
          write(lp,9200) iarec(iaunit),iaunit
          call flush(lp)
          call xchalt('(zaiord)')
                 stop '(zaiord)'
        endif
      else
!$OMP   PARALLEL DO PRIVATE(j,i)
!$OMP&              REDUCTION(MIN:wmin) REDUCTION(MAX:wmax)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j= 1,jtdm
          do i= 1,itdm
            atmp(i,j) = w(i+(j-1)*itdm)
            if     (w(i+(j-1)*itdm).ne.spval) then
              wmin = min( wmin, w(i+(j-1)*itdm) )
              wmax = max( wmax, w(i+(j-1)*itdm) )
            endif
          enddo
        enddo
      endif
      rmin(1) = wmin
      rmax(1) = wmax
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
      call xcaput(atmp,h, 1)
      call xcminr(rmin)
      call xcmaxr(rmax)
      hmin = rmin(1)
      hmax = rmax(1)
#if defined(TIMER)
c
      call xctmr1(17)
#endif
      return
c
 9000 format(/ /10x,'error in zaiord -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
 9100 format(/ /10x,'error in zaiord -  can''t read record',
     &   i4,' on array I/O unit ',i3,'.'/ /)
 9200 format(/ /10x,'error in zaiord -  record',
     &   i4,' on array I/O unit ',i3,
     &   ' has 2.0**100 outside masked region.'/ /)
      end subroutine zaiord

      subroutine zaiordd(a,n, iunit,irec,ios)
      implicit none
c
      integer, intent(in)    :: n,iunit,irec
      integer, intent(out)   :: ios
      real*4,  intent(out)   :: a(n)
c
c**********
c*
c 1)  direct access read a single record.
c
c 2)  expressed as a subroutine because i/o with 
c     implied do loops can be slow on some machines.
c*
c**********
c
      read(unit=iunit, rec=irec, iostat=ios) a
#if defined(ENDIAN_IO)
      call zaio_endian(a,n)
#endif
      return
      end subroutine zaiordd

      subroutine zaiosk(iaunit)
      implicit none
c
      integer, intent(in)    :: iaunit
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for skipping an array read.
c
c     must call zaiopn for this array unit before calling zaiosk.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zaiopn.
c*
c**********
#if defined(TIMER)
c
      call xctmr0(16)
#endif
c
      if     (iarec(iaunit).lt.0) then
        write(lp,9000) iaunit
        call xcstop('(zaiosk)')
               stop '(zaiosk)'
      endif
c
      iarec(iaunit) = iarec(iaunit) + 1
c
c --- I/O from first processor only
c
      if     (mnproc.eq.1) then
*
*     write(lp,*) 'zaiosk - iaunit,rec = ',iaunit,iarec(iaunit)
*     call flush(lp)
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
#if defined(TIMER)
c
      call xctmr1(16)
#endif
      return
c
 9000 format(/ /10x,'error in zaiosk -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
      end subroutine zaiosk

      subroutine zaiowr3(h, l, mask,lmask, hmin,hmax, iaunit, lreal4)
      implicit none
c
      logical, intent(in)    :: lmask,lreal4
      integer, intent(in)    :: l,iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin(l),hmax(l)
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l),
     &         intent(inout) :: h
#else
      real,    intent(out)   :: hmin(l),hmax(l)
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l),
     &         intent(inout) :: h
#endif
c
c**********
c*
c  1) machine specific routine for 3-d array writing.
c
c     must call zaiopn for this array unit before calling zaiord.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zaiopn.
c
c  4) hmin,hmax are returned as the minimum and maximum value in the array.
c     if lmask==.true. the range is only where mask.ne.0, with all other
c     values output as 2.0**100.
c
c  5) If lreal4==.true. then h is overwritten on exit with real*4 version
c     of the same array.  This is typically used for reproducability on
c     restart.
c*
c**********
c
c     this version just calls zaiowr l times.
c
      integer k
c
      do k= 1,l
        call zaiowr(h(1-nbdy,1-nbdy,k), mask,lmask,
     &              hmin(k),hmax(k), iaunit, lreal4)
      enddo
      return
      end subroutine zaiowr3

      subroutine zaiowr(h, mask,lmask, hmin,hmax,  iaunit, lreal4)
      implicit none
c
      logical, intent(in)    :: lmask,lreal4
      integer, intent(in)    :: iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin,hmax
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(inout) :: h
#else
      real,    intent(out)   :: hmin,hmax
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(inout) :: h
#endif
c
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
c
c**********
c*
c  1) machine specific routine for array writing.
c
c     must call zaiopn for this array unit before calling zaiord.
c
c  2) array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zaiopn.
c
c  4) hmin,hmax are returned as the minimum and maximum value in the array.
c     if lmask==.true. the range is only where mask.ne.0, with all other
c     values output as 2.0**100.
c
c  5) If lreal4==.true. then h is overwritten on exit with real*4 version
c     of the same array.  This is typically used for reproducability on
c     restart.
c*
c**********
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
c --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
      integer    n2drec
      parameter (n2drec=((itdm*jtdm+4095)/4096)*4096)
c
      real*4         w
      common/czioxw/ w(n2drec)
      save  /czioxw/
c
      real           atmp,amsk,tmsk
      common/czioxr/ atmp(itdm,jtdm),amsk(itdm,jtdm),
     &               tmsk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      save  /czioxr/
c
      integer   ios, i,j
      real*4    wmin,wmax
      real      rmin(1),rmax(1),vsave
#if defined(TIMER)
c
      call xctmr0(18)
#endif
c
      if     (iarec(iaunit).lt.0) then
        write(lp,9000) iaunit
        call xcstop('(zaiowr)')
               stop '(zaiowr)'
      endif
c
      iarec(iaunit) = iarec(iaunit) + 1
      rmin(1) =  spval
      rmax(1) = -spval
c
c --- I/O from first processor only
c
      if     (lmask) then
        tmsk(1:ii,1:jj) = mask(1:ii,1:jj)
        vsave = vland
        vland = 0.0
        call xcaget(amsk,tmsk, 1)
        vland = spval
        call xcaget(atmp,h,    1)
        vland = vsave
      else
        vsave = vland
        vland = spval
        call xcaget(atmp,h,    1)
        vland = vsave
      endif
c
      if     (mnproc.eq.1) then
c
      wmin =  spval
      wmax = -spval
      if     (lreal4) then
        if     (lmask) then
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&                REDUCTION(MIN:wmin) REDUCTION(MAX:wmax)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j= 1,jtdm
            do i= 1,itdm
              if     (amsk(i,j).ne.0.0) then
                w(i+(j-1)*itdm) = atmp(i,j)
                wmin = min( wmin, w(i+(j-1)*itdm) )
                wmax = max( wmax, w(i+(j-1)*itdm) )
              else
                w(i+(j-1)*itdm) = spval
              endif
#if defined(REAL4)
! ---         atmp(i,j) = w(i+(j-1)*itdm)  ! h is already real*4
#else
              atmp(i,j) = w(i+(j-1)*itdm)  ! h is not real*4, so update it
#endif
            enddo
          enddo
        else
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&                REDUCTION(MIN:wmin) REDUCTION(MAX:wmax)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j= 1,jtdm
            do i= 1,itdm
              w(i+(j-1)*itdm) = atmp(i,j)
              if     (w(i+(j-1)*itdm).ne.spval) then
                wmin = min( wmin, w(i+(j-1)*itdm) )
                wmax = max( wmax, w(i+(j-1)*itdm) )
              endif
#if defined(REAL4)
! ---         atmp(i,j) = w(i+(j-1)*itdm)  ! h is already real*4
#else
              atmp(i,j) = w(i+(j-1)*itdm)  ! h is not real*4, so update it
#endif
            enddo
          enddo
        endif
      else
        if     (lmask) then
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&                REDUCTION(MIN:wmin) REDUCTION(MAX:wmax)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j= 1,jtdm
            do i= 1,itdm
              if     (amsk(i,j).ne.0.0) then
                w(i+(j-1)*itdm) = atmp(i,j)
                wmin = min( wmin, w(i+(j-1)*itdm) )
                wmax = max( wmax, w(i+(j-1)*itdm) )
              else
                w(i+(j-1)*itdm) = spval
              endif
            enddo
          enddo
        else
!$OMP     PARALLEL DO PRIVATE(j,i)
!$OMP&                REDUCTION(MIN:wmin) REDUCTION(MAX:wmax)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j= 1,jtdm
            do i= 1,itdm
              w(i+(j-1)*itdm) = atmp(i,j)
              if     (w(i+(j-1)*itdm).ne.spval) then
                wmin = min( wmin, w(i+(j-1)*itdm) )
                wmax = max( wmax, w(i+(j-1)*itdm) )
              endif
            enddo
          enddo
        endif
      endif
      do i= itdm*jtdm+1,n2drec
        w(i) = spval
      enddo
      rmin(1) = wmin
      rmax(1) = wmax
      call zaiowrd(w,n2drec, iaunit+1000,iarec(iaunit),ios)
      if     (ios.ne.0) then
        write(lp,9100) iarec(iaunit),iaunit
        call flush(lp)
        call xchalt('(zaiowr)')
               stop '(zaiowr)'
      endif
c
      endif  ! I/O from first processor only
      call xcsync(no_flush)
      if     (lreal4) then
#if defined(REAL4)
! ---   call xcaput(atmp,h, 1)  ! h is already real*4
#else
        call xcaput(atmp,h, 1)  ! h is not real*4, so update it
#endif
      endif
      call xcminr(rmin)
      call xcmaxr(rmax)
      hmin = rmin(1)
      hmax = rmax(1)
#if defined(TIMER)
c
      call xctmr1(18)
#endif
      return
c
 9000 format(/ /10x,'error in zaiowr -  array I/O unit ',
     &   i3,' is not marked as open.'/ /)
 9100 format(/ /10x,'error in zaiowr -  can''t write record',
     &   i4,' on array I/O unit ',i3,'.'/ /)
      end subroutine zaiowr

      subroutine zaiowrd(a,n, iunit,irec,ios)
      implicit none
c
      integer, intent(in)    :: n,iunit,irec
      integer, intent(out)   :: ios
      real*4,  intent(in)    :: a(n)
c
c**********
c*
c 1)  direct access write a single record.
c
c 2)  expressed as a subroutine because i/o with 
c     implied do loops can be slow on some machines.
c*
c**********
c
#if defined(ENDIAN_IO)
      call zaio_endian(a,n)  ! overwrites a
#endif
      write(unit=iunit, rec=irec, iostat=ios) a
      return
      end subroutine zaiowrd
