      module mod_pipe
      use mod_xc  ! HYCOM communication interface
c
c --- HYCOM (named pipe based) debugging interface
c
      logical, save, public  :: lpipe
c
      integer, save, private :: ipunit,lpunit
      logical, save, private :: lslave,ldebug
c
      real,    allocatable, dimension(:,:),
     &         save, private :: field1,field2,tmask,amask

      contains
c
c --- this set of routines facilitates output comparison from two HYCOM
c --- versions running side by side. one model, the 'slave', writes its
c --- output into a named pipe. the other model, the 'master', reads
c --- from the pipe and compares.
c --- differences are recorded in 'PIPE_base.out'.
c
c --- call 'pipe_init' at start of main program.
c
c ---   if the file 'PIPE_MASTER' exists then this is the master,
c ---   if the file 'PIPE_SLAVE'  exists then this is the slave,
c ---   otherwise there is no comparison made.
c
c ---   if the file 'PIPE_DEBUG' exists then debugging printout
c ---   is produced for point itest,jtest.  This works with or
c ---   without a pipe for comparison.
c
c ---   the 'PIPE_MASTER' and 'PIPE_SLAVE' files contain the location
c ---   of an existing named-pipe.  The contents of the 'PIPE_DEBUG'
c ---   file is ignored.
c
c --- call 'pipe_compare'    (from master and slave) anywhere in the code
c ---   to check whether data stored in a single array are identical
c
c --- call 'pipe_compareall' (from master and slave) after major routines
c ---   to check whether data stored in all major arrays are identical
c
      subroutine pipe_init
      implicit none
c
      character*256 cpipe
c
      character*12  cinfo
      integer       irecl
c
      inquire(file='PIPE_MASTER',exist=lpipe)
      inquire(file='PIPE_SLAVE', exist=lslave)
      inquire(file='PIPE_DEBUG', exist=ldebug)
      if     (lpipe .and. lslave) then
        stop 'pipe_init: (master/slave ambiguity)'
      endif
      lpipe = lpipe .or. lslave
c
      if     (lpipe) then
c
c ---   open the pipe and some output files, allocate arrays
c
        allocate( field1(itdm,jtdm) )
        allocate( field2(itdm,jtdm) )
        allocate( tmask( itdm,jtdm) )
        allocate( amask( 1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
c
        ipunit=18
        lpunit=19
c
        if     (mnproc.eq.1) then
          if     (.not.lslave) then
            open (unit=17,file='PIPE_MASTER',status='old',
     &            form='formatted')
            read (     17,'(a)') cpipe
            close(unit=17)
            write(lp,'(a,a)') 'master opening pipe for reading: ',
     &                        cpipe(1:len_trim(cpipe))
            call flush(lp)
            cinfo=' '  !removes spurious compiler warning message
            inquire( iolength=irecl ) cinfo,field1
            call flush(lp)
            open (unit=ipunit,file=cpipe,status='old',
     &            action='read',recl=irecl,
     &            form='unformatted')
            open (unit=lpunit,file='PIPE_base.out',status='unknown')
          else
            open (unit=17,file='PIPE_SLAVE', status='old',
     &            form='formatted')
            read (     17,'(a)') cpipe
            close(unit=17)
            write(lp,'(a,a)') 'slave opening pipe for writing: ',
     &                        cpipe(1:len_trim(cpipe))
            call flush(lp)
            cinfo=' '  !removes spurious compiler warning message
            inquire( iolength=irecl ) cinfo,field1
            open (unit=ipunit,file=cpipe,status='old',
     &            action='write',recl=irecl,
     &            form='unformatted')
            open (unit=lpunit,file='PIPE_test.out',status='unknown')
          endif
        endif  !1st tile only.
        call xcsync(flush_lp)
      endif
c
      return
      end subroutine pipe_init
c
      subroutine pipe_compare(field,mask,what)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: field
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &         intent(in) :: mask
      character*12,
     &         intent(in) :: what
c
c --- call this routine from anywhere in the code (from both versions, of
c --- course) to check whether data stored in 'field' are identical
c
      integer      i,j
      logical      fail
      character*12 which
c
      if (lpipe) then
        do j=1,jj
          do i=1,ii
            amask(i,j) = mask(i,j)
          enddo
        enddo
        call xcaget(tmask,  amask, 1)
        call xcaget(field2, field, 1)
        if     (mnproc.eq.1) then
          if (lslave) then
            write (lpunit,'(2a)') 'writing for comparison: ',what
            write (ipunit) what, field2
          else				!  lslave = .false.
            read  (ipunit) which,field1
            write (lpunit,'(2a)') 'reading for comparison: ',which
            if (what.ne.which) then
              write (lpunit,'(4a)') 'out of sync -- trying to compare ',
     &           what,'  to  ',which
              stop
            endif
c
            fail=.false.
            do j=1,jtdm
              do i=1,itdm
                if (tmask(i,j).gt.0.0 .and.
     &              field2(i,j).ne.field1(i,j)) then
                  write (lpunit,'(a,2i5,1p,2(a,e12.5),4x,a)')
     &              'i,j=',i,j,
     &              '  master:',field2(i,j),
     &              '  error:', field2(i,j)-field1(i,j),what
                  fail=.true.
                endif
              enddo
            enddo
            if (fail) then  ! optional
              stop 'in compare'
            endif
          endif !slave:master
        endif !1st tile
        call xcsync(no_flush) ! wait for 1st tile
      endif !lpipe
      return
      end subroutine pipe_compare

      subroutine pipe_comparall(m,n, cinfo)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer, intent(in) :: m,n
      character*12,
     &         intent(in) :: cinfo
c
c --- write out a standard menu of arrays for testing
c
      character*12 text
      integer i,imax,imin,j,jmax,jmin,k,l,mnp
      real    diso,dmax,dmin,damax,damin
      real*8  tmean,smean,pmean,rmean
      real*8  d1,d2,d3,d4
c
      real*8  tmean0,smean0,rmean0,
     &        tmean1,smean1,rmean1
      save    tmean0,smean0,rmean0,
     &        tmean1,smean1,rmean1
      data    tmean0,smean0,rmean0 / 3*0.0d0 /
c
cdiag if     (mnproc.eq.1) then
cdiag write(lp,'(a,i10)') cinfo,nstep
cdiag call flush(lp)
cdiag endif
c
 103  format (i8,i5,i4,1x,a,a/
     &       (i8,5x,i4,1x,a,a,2f7.3,2f7.3,f8.4,f9.3,f9.2))
c
      if     (ldebug .and. ittest.ne.-1 .and. jttest.ne.-1) then
        if     (min(ittest,jttest).le.0) then
          stop 'comparall: (debug ambiguity)'
        endif
        if     (i0.lt.ittest .and. i0+ii.ge.ittest .and.
     &          j0.lt.jttest .and. j0+jj.ge.jttest      ) then
        if     (.not.mxlkrt) then
          write (lp,103) 
     &       nstep,itest+i0,jtest+j0,cinfo(1:6),
     &       ':   utot   vtot   temp   saln    dens    thkns     dpth',
     &      (nstep,k,          cinfo(1:6),':',
     &          u(itest,jtest,k,n)+ubavg(itest,jtest,n),
     &          v(itest,jtest,k,n)+vbavg(itest,jtest,n),
     &       temp(itest,jtest,k,n),
     &       saln(itest,jtest,k,n),
     &       th3d(itest,jtest,k,n)+thbase,
     &         dp(itest,jtest,k,n)/onem,
     &          p(itest,jtest,k+1 )/onem,
     &       k=1,kk)
        else
c ---     include KT mixed layer values.
          write (lp,103) 
     &       nstep,itest+i0,jtest+j0,cinfo(1:6),
     &       ':   utot   vtot   temp   saln    dens    thkns     dpth',
     &       nstep,0,          cinfo(1:6),':',
     &          0.0,
     &          0.0,
     &       tmix(itest,jtest),     
     &       smix(itest,jtest),     
     &      thmix(itest,jtest)+thbase,
     &     dpmixl(itest,jtest,n)/onem,
     &     dpmixl(itest,jtest,n)/onem,
     &      (nstep,k,          cinfo(1:6),':',
     &          u(itest,jtest,k,n)+ubavg(itest,jtest,n),
     &          v(itest,jtest,k,n)+vbavg(itest,jtest,n),
     &       temp(itest,jtest,k,n),
     &       saln(itest,jtest,k,n),
     &       th3d(itest,jtest,k,n)+thbase,
     &         dp(itest,jtest,k,n)/onem,
     &          p(itest,jtest,k+1 )/onem,
     &       k=1,kk)
        endif
        endif  ! ittest,jttest tile
        call xcsync(flush_lp)
      endif
c
      if     (ldebug .and. ittest.eq.-1) then
c
c ---   printout min/max/iospycnal th3d
c
 104    format (i8,a3,1x,a,a)
 105    format (i8,i3,1x,a,a,2i5,f9.5,f7.3,f9.5,2i5,i7)
        if     (mnproc.eq.1) then
        write(lp,104)
     &       nstep,'  k',cinfo(1:6),
     &       ': imin jmin  denamin deniso  denamax imax jmax mnproc'
        endif
        call xcsync(flush_lp)
        do k= 1,kk
          diso=theta(k)
          dmin= huge
          dmax=-huge
          do j= 1,jj
            do i= 1,ii
              if     (ip(i,j).eq.1) then
                if     (th3d(i,j,k,n).lt.dmin) then
                  dmin=th3d(i,j,k,n)
                  imin=i
                  jmin=j
                endif
                if     (th3d(i,j,k,n).gt.dmax) then
                  dmax=th3d(i,j,k,n)
                  imax=i
                  jmax=j
                endif
              endif
            enddo
          enddo
          damin=dmin
          call xcminr(damin)
          damax=dmax
          call xcmaxr(damax)
          do mnp= 1,ijpr
            if     (mnp.eq.mnproc) then
              if     (dmin.eq.damin .or. dmax.eq.damax) then
                write (lp,105) 
     &            nstep,k,cinfo(1:6),
     &            ':',imin,jmin,dmin-diso,
     &                          diso+thbase,
     &                          dmax-diso,imax,jmax,mnproc
              endif
            endif
            call xcsync(flush_lp)
          enddo
        enddo
        call flush(lp)
      endif
c
      if     (ldebug .and. jttest.eq.-1) then
c
c ---   printout basin-wide means.
c
!$OMP   PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do k=1,kk
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=             dp(i,j,k,n)*scp2(i,j)
                util2(i,j)=temp(i,j,k,n)*dp(i,j,k,n)*scp2(i,j)
                util3(i,j)=saln(i,j,k,n)*dp(i,j,k,n)*scp2(i,j)
                util4(i,j)=th3d(i,j,k,n)*dp(i,j,k,n)*scp2(i,j)
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
c
 106    format (i8,3x,1x,a,a,3f8.4,1p3e10.2)
        if     (mnproc.eq.1) then
        write (lp,106)
     &    nstep,cinfo(1:6),
     &    ': t,s,th',
     &    tmean,smean,rmean+thbase,
     &    tmean-tmean0,smean-smean0,rmean-rmean0
        call flush(lp)
        endif
c
        if     (cinfo(1:6).eq.'ENTER ') then
          tmean1=tmean
          smean1=smean
          rmean1=rmean
        elseif (cinfo(1:6).eq.'tsadvc') then
          if     (mnproc.eq.1) then
          write (lp,106)
     &      nstep,'cn+tsa',
     &      ': t,s,th',
     &      tmean,smean,rmean+thbase,
     &      tmean-tmean1,smean-smean1,rmean-rmean1
          call flush(lp)
          endif
        elseif (cinfo(1:6).eq.'hybgen') then
          if     (mnproc.eq.1) then
          write (lp,106)
     &      nstep,'EXIT  ',
     &      ': t,s,th',
     &      tmean,smean,rmean+thbase,
     &      tmean-tmean1,smean-smean1,rmean-rmean1
          call flush(lp)
          endif
        endif
c
        tmean0=tmean
        smean0=smean
        rmean0=rmean
      endif
c
      if     (lpipe) then
        write (lpunit,'(a,i10)') cinfo,nstep
        text='ubavg(n)    '
        call pipe_compare(ubavg(1-nbdy,1-nbdy,n),iu,text)
        text='vbavg(n)    '
        call pipe_compare(vbavg(1-nbdy,1-nbdy,n),iv,text)
        text='pbavg(n)    '
        call pipe_compare(pbavg(1-nbdy,1-nbdy,n),ip,text)
        text='montg(1)    '
        call pipe_compare(montg(1-nbdy,1-nbdy,1),ip,text)
        do k=1,kk
          write (text,'(a9,i3)') 'u(kn)  k=',k
          call pipe_compare(u(1-nbdy,1-nbdy,k,n),iu,text)
          write (text,'(a9,i3)') 'v(kn)  k=',k
          call pipe_compare(v(1-nbdy,1-nbdy,k,n),iv,text)
          write (text,'(a9,i3)') 'dp(kn) k=',k
          call pipe_compare(dp(1-nbdy,1-nbdy,k,n),ip,text)
          write (text,'(a9,i3)') 'temp(kn) ',k
          call pipe_compare(temp(1-nbdy,1-nbdy,k,n),ip,text)
          write (text,'(a9,i3)') 'saln(kn) ',k
          call pipe_compare(saln(1-nbdy,1-nbdy,k,n),ip,text)
          write (text,'(a9,i3)') 'th3d(kn) ',k
          call pipe_compare(th3d(1-nbdy,1-nbdy,k,n),ip,text)
        enddo
        write (lpunit,'(a,i10,a)') cinfo,nstep,' -- OK'
      endif
c
      return
      end subroutine pipe_comparall
c
      end module mod_pipe
c
c
c> Revision history:
c>
c> Oct  2000 - added PIPE_DEBUG for debugging printout
