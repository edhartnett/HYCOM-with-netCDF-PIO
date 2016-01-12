      subroutine restart_in(nstep0, dtime0)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
c
      integer nstep0
      real*8  dtime0
c
c     read in a restart file on unit 11.
c
      integer   i,ios,j,k,ktr
      character cline*80
c
      call zaiopf(flnmrsi(1:len_trim(flnmrsi))//'.a','old', 11)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=11,file=flnmrsi(1:len_trim(flnmrsi))//'.b',
     &        status='old',action='read',form='formatted')
      endif
      call zagetc(cline,ios, 11)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)')
     &      'I/O error from zagetc, iunit,ios = ',11,ios
        endif !1st tile
        call xcstop('(restart_in)')
               stop '(restart_in)'
      endif
      if     (mnproc.eq.1) then
      write(lp,'(a)') cline(1:len_trim(cline))
      endif !1st tile
      if     (cline(1:9).ne.'RESTART:') then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') 'error in hycom - ice restart failed'
        endif !1st tile
        call xcstop('(restart_in)')
               stop '(restart_in)'
      endif
      call zagetc(cline,ios, 11)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)')
     &      'I/O error from zagetc, iunit,ios = ',11,ios
        endif !1st tile
        call xcstop('(restart_in)')
               stop '(restart_in)'
      endif
      if     (mnproc.eq.1) then
      write(lp,'(a)') cline(1:len_trim(cline))
      call flush(lp)
      endif !1st tile
      i = index(cline,'=')
      read(cline(i+1:),*) nstep0,dtime0
c
      call restart_in3d(u,     2*kdm, iu, 'u       ')
      call restart_in3d(v,     2*kdm, iv, 'v       ')
      call restart_in3d(dp,    2*kdm, ip, 'dp      ')
      call restart_in3d(temp,  2*kdm, ip, 'temp    ')
      call restart_in3d(saln,  2*kdm, ip, 'saln    ')
      call restart_in3d(th3d,  2*kdm, ip, 'th3d    ')
      if (mxlmy) then
        call restart_in3d(q2    ,2*kdm+4, ip, 'q2      ')
        call restart_in3d(q2l   ,2*kdm+4, ip, 'q2l     ')
        call restart_in3d(vctymy,  kdm+2, ip, 'vctymy  ')
        call restart_in3d(difqmy,  kdm+2, ip, 'difqmy  ')
        call restart_in3d(diftmy,  kdm+2, ip, 'diftmy  ')
      endif
      call restart_in3d(ubavg,     3, iu, 'ubavg   ')
      call restart_in3d(vbavg,     3, iv, 'vbavg   ')
      call restart_in3d(pbavg,     3, ip, 'pbavg   ')
      call restart_in3d(pbot,      1, ip, 'pbot    ')
      call restart_in3d(psikk,     1, ip, 'psikk   ')
      call restart_in3d(thkk,      1, ip, 'thkk    ')
      call restart_in3d(dpmixl,    2, ip, 'dpmixl  ')
      if (icegln) then
        call zagetc(cline,ios, 11)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)')
     &        'I/O error from zagetc, iunit,ios = ',11,ios
          endif !1st tile
          call xcstop('(restart_in)')
                 stop '(restart_in)'
        endif
        if     (ios.ne.0 .or. cline(1:8).ne.'temice  ') then
c
c ---     assume this is an addition of ice to the simulation.
c
          if     (mnproc.eq.1) then
          write(lp,'(/ a /)') 'adding ice to the simulation.'
          call flush(lp)
          endif !1st tile
c
          do j= 1,jj
            do i= 1,ii
              temice(i,j) = temp(i,j,1,1)
              covice(i,j) = 0.0
              thkice(i,j) = 0.0
            enddo
          enddo
          if     (trcrin .and. cline(1:8).eq.'tracer  ') then
c ---       reposition file for tracer input
            if     (mnproc.eq.1) then  ! .b file from 1st tile only
              rewind(11)
            endif
            do k= 1,12*kdm+16
              call zagetc(cline,ios, 11)
              if     (ios.ne.0) then
                if     (mnproc.eq.1) then
                  write(lp,'(/ a,i4,i9 /)')
     &              'I/O error from zagetc, iunit,ios = ',11,ios
                endif !1st tile
                call xcstop('(restart_in)')
                       stop '(restart_in)'
              endif
            enddo
          endif
        else
c
c ---     reposition file for ice input
c
          if     (mnproc.eq.1) then  ! .b file from 1st tile only
            rewind(11)
          endif
          do k= 1,12*kdm+16
            call zagetc(cline,ios, 11)
            if     (ios.ne.0) then
              if     (mnproc.eq.1) then
                write(lp,'(/ a,i4,i9 /)')
     &            'I/O error from zagetc, iunit,ios = ',11,ios
              endif !1st tile
              call xcstop('(restart_in)')
                     stop '(restart_in)'
            endif
*           if     (mnproc.eq.1) then
*           write(lp,'(a)') cline
*           endif !1st tile
          enddo
          call restart_in3d(temice,    1, ip, 'temice  ')
          call restart_in3d(covice,    1, ip, 'covice  ')
          call restart_in3d(thkice,    1, ip, 'thkice  ')
        endif
      endif
      if (trcrin) then
        do ktr= 1,ntracr
          call restart_in3d(tracer(1-nbdy,1-nbdy,1,1,ktr),
     &                             2*kdm, ip, 'tracer  ')
        enddo
      endif
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        close (unit=11)
      endif
      call zaiocl(11)
      return
      end subroutine restart_in

      subroutine restart_in3d(field,l, mask, cfield)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      integer   l
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l) ::
     & field
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & mask
      character cfield*8
c
c --- read a single restart 3-d array field from unit 11.
c
      integer   i,ios,layer,level,k
      real      hmina(2*kdm+4),hminb,hmaxa(2*kdm+4),hmaxb
      character cline*80
c
*     if     (mnproc.eq.1) then
*     write(lp,'(a,i3,2x,a)') 'restart_in3d - l,cfield = ',l,cfield
*     call flush(lp)
*     endif !1st tile
      call zaiord3(field,l, mask,.false., hmina,hmaxa, 11)
c
      do k= 1,l
        call zagetc(cline,ios, 11)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)')
     &        'I/O error from zagetc, iunit,ios = ',11,ios
          endif !1st tile
          call xcstop('(restart_in3d)')
                 stop '(restart_in3d)'
        endif
*       if     (mnproc.eq.1) then
*       write (lp,'(a)')  cline(1:len_trim(cline))
*       endif !1st tile
        if     (cline(1:8).ne.cfield) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cline(1:len_trim(cline)),
     &           'error in restart_in3d - expected ',cfield
          endif !1st tile
          call xcstop('(restart_in3d)')
                 stop '(restart_in3d)'
        endif
        i = index(cline,'=')
        read (cline(i+1:),*) layer,level,hminb,hmaxb
        if     (abs(hmina(k)-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa(k)-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,3i3 / a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      'iunit,k,l = ',11,k,l,
     &      cline,
     &      '.a,.b min = ',hmina(k),hminb,hmina(k)-hminb,
     &      '.a,.b max = ',hmaxa(k),hmaxb,hmaxa(k)-hmaxb
          endif !1st tile
          call xcstop('(restart_in3d)')
                 stop '(restart_in3d)'
        endif
      enddo
c
      return
      end subroutine restart_in3d

      subroutine restart_out(nstepx, dtimex)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
c
      integer nstepx
      real*8  dtimex
c
c     write out in a restart file on unit 12 (and a flux file on unit 25).
c
      logical   lopen
      integer   i,j,k,ktr,l
      real      xmin(2*kdm+4),xmax(2*kdm+4)
      character cline*80
c
      call zaiopi(lopen, 12)
      if (.not.lopen) then
        call zaiopf(flnmrso(1:len_trim(flnmrso))//'.a','new', 12)
        if     (mnproc.eq.1) then
          open (unit=12,file=flnmrso(1:len_trim(flnmrso))//'.b',
     &          status='new',action='write',form='formatted')
          write(lp,'(a)') ' creating a new restart file'
          call flush(lp)
        endif !1st tile
      else
        call zaiorw(12)
        if     (mnproc.eq.1) then
          rewind(unit=12)
          write(lp,'(a)') ' over-writing any previous restart'
          call flush(lp)
        endif !1st tile
      endif
c
      if     (mnproc.eq.1) then
      write(12,'(a,3i6)') 'RESTART: iexpt,iversn,yrflag = ',
     &                              iexpt,iversn,yrflag
      write(cline,*)                nstepx,dtimex
      write(12,'(a,a)')   'RESTART: nstep,dtime = ',
     &                              cline(1:len_trim(cline))
      call flush(12)
      endif !1st tile
c
      call zaiowr3(u,      2*kdm, iu,.false., xmin,xmax, 12, .true.)
      call xctilr( u,    1,2*kdm, nbdy,nbdy, halo_uv)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'u       ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(v,      2*kdm, iv,.false., xmin,xmax, 12, .true.)
      call xctilr( v,    1,2*kdm, nbdy,nbdy, halo_vv)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'v       ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(dp,     2*kdm, ip,.false., xmin,xmax, 12, .true.)
      call xctilr( dp,   1,2*kdm, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'dp      ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(temp,   2*kdm, ip,.false., xmin,xmax, 12, .true.)
      call xctilr( temp, 1,2*kdm, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'temp    ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(saln,   2*kdm, ip,.false., xmin,xmax, 12, .true.)
      call xctilr( saln, 1,2*kdm, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'saln    ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(th3d,   2*kdm, ip,.false., xmin,xmax, 12, .true.)
      call xctilr( th3d, 1,2*kdm, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(12,4100) 'th3d    ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(12)
      endif !1st tile
c
      if (mxlmy) then
        call zaiowr3(q2, 2*kdm+4, ip,.false., xmin,xmax, 12, .true.)
        if     (mnproc.eq.1) then
        do l= 0,1
          do k= 1,kdm+2
            write(12,4100) 'q2      '
     &                    ,k,l+1,xmin(k+l*(kdm+2)),xmax(k+l*(kdm+2))
          enddo
        enddo
        call flush(12)
        endif !1st tile
        call zaiowr3(q2l, 2*kdm+4, ip,.false., xmin,xmax, 12, .true.)
        if     (mnproc.eq.1) then
        do l= 0,1
          do k= 1,kdm+2
            write(12,4100) 'q2l     '
     &                    ,k,l+1,xmin(k+l*(kdm+2)),xmax(k+l*(kdm+2))
          enddo
        enddo
        call flush(12)
        endif !1st tile
        call zaiowr3(vctymy, kdm+2, ip,.false., xmin,xmax, 12, .true.)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 1,kdm+2
            write(12,4100) 'vctymy  ',k,l,xmin(k),xmax(k)
          enddo
        enddo
        call flush(12)
        endif !1st tile
        call zaiowr3(difqmy, kdm+2, ip,.false., xmin,xmax, 12, .true.)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 1,kdm+2
            write(12,4100) 'difqmy  ',k,l,xmin(k),xmax(k)
          enddo
        enddo
        call flush(12)
        endif !1st tile
        call zaiowr3(diftmy, kdm+2, ip,.false., xmin,xmax, 12, .true.)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 1,kdm+2
            write(12,4100) 'diftmy  ',k,l,xmin(k),xmax(k)
          enddo
        enddo
        call flush(12)
        endif !1st tile
      endif
c
      call zaiowr3(ubavg,      3, iu,.false., xmin,xmax, 12, .true.)
      call xctilr( ubavg,    1,3, nbdy,nbdy, halo_uv)
      if     (mnproc.eq.1) then
      do l= 1,3
        do k= 0,0
          write(12,4100) 'ubavg   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(vbavg,      3, iv,.false., xmin,xmax, 12, .true.)
      call xctilr( vbavg,    1,3, nbdy,nbdy, halo_vv)
      if     (mnproc.eq.1) then
      do l= 1,3
        do k= 0,0
          write(12,4100) 'vbavg   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(pbavg,      3, ip,.false., xmin,xmax, 12, .true.)
      call xctilr( pbavg,    1,3, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 1,3
        do k= 0,0
          write(12,4100) 'pbavg   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(pbot,       1, ip,.false., xmin,xmax, 12, .true.)
      call xctilr( pbot,     1,1, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 1,1
        do k= 0,0
          write(12,4100) 'pbot    ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(psikk,      1, ip,.false., xmin,xmax, 12, .true.)
      call xctilr( psikk,    1,1, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 1,1
        do k= 0,0
          write(12,4100) 'psikk   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(thkk,       1, ip,.false., xmin,xmax, 12, .true.)
      call xctilr( thkk,     1,1, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 1,1
        do k= 0,0
          write(12,4100) 'thkk    ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      call zaiowr3(dpmixl,     2, ip,.false., xmin,xmax, 12, .true.)
      call xctilr( dpmixl,   1,2, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 1,2
        do k= 0,0
          write(12,4100) 'dpmixl  ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(12)
      endif !1st tile
      if (icegln) then
        call zaiowr3(temice,     1, ip,.false., xmin,xmax, 12, .true.)
        call xctilr( temice,   1,1, nbdy,nbdy, halo_ps)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(12,4100) 'temice  ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(12)
        endif !1st tile
        call zaiowr3(covice,     1, ip,.false., xmin,xmax, 12, .true.)
        call xctilr( covice,   1,1, nbdy,nbdy, halo_ps)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(12,4100) 'covice  ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(12)
        endif !1st tile
        call zaiowr3(thkice,     1, ip,.false., xmin,xmax, 12, .true.)
        call xctilr( thkice,   1,1, nbdy,nbdy, halo_ps)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(12,4100) 'thkice  ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(12)
        endif !1st tile
      endif
      if (trcout) then
        do ktr= 1,ntracr
          call zaiowr3(tracer(1-nbdy,1-nbdy,1,1,ktr),  2*kdm,
     &                 ip,.false., xmin,xmax, 12, .true.)
          call xctilr( tracer(1-nbdy,1-nbdy,1,1,ktr),1,2*kdm,
     &                  nbdy,nbdy, halo_ps)
          if     (mnproc.eq.1) then
          do l= 0,1
            do k= 1,kdm
              write(12,4100) 'tracer  ',k,l+1,xmin(k+l*kdm),
     &                                        xmax(k+l*kdm)
            enddo
          enddo
          call flush(12)
          endif !1st tile
        enddo !ktr
      endif !trcout
      call zaiofl(12)
      if     (mnproc.eq.1) then
      call flush( 12)
      write(lp,'(a,f11.2)') ' restart created at model day',dtimex
      call flush(12)
      endif !1st tile
      call xcsync(flush_lp)
c
c --- output to flux file
c
      if (.FALSE.) then ! turn on/off flux output
        call zaiopi(lopen, 25)
        if (.not.lopen) then
          call zaiopf(flnmflx(1:len_trim(flnmflx))//'.a','new', 25)
          if     (mnproc.eq.1) then
          open (unit=25,file=flnmflx(1:len_trim(flnmflx))//'.b',
     &          status='new',action='write',form='formatted')
          write(25,'(a,3i6)') 'FLUXES: iexpt,iversn,yrflag = ',
     &                                 iexpt,iversn,yrflag
          call flush(25)
          endif !1st tile
        endif
c
        if     (mnproc.eq.1) then
        write(cline,*)               nstepx,dtimex
        write(25,'(a,a)')   'FLUXES: nstep,dtime = ',
     &                               cline(1:len_trim(cline))
        call flush(25)
        endif !1st tile
c
        call zaiowr3(dpav,     kdm, ip,.true.,  xmin,xmax, 25, .false.)
        if     (mnproc.eq.1) then
        do l= 0,0
          do k= 1,kdm
            write(25,4100) 'dpav    ',k,l,  xmin(k+l*kdm),xmax(k+l*kdm)
          enddo
        enddo
        call flush(25)
        endif !1st tile
        call zaiowr3(uflxav,   kdm, iu,.true.,  xmin,xmax, 25, .false.)
        if     (mnproc.eq.1) then
        do l= 0,0
          do k= 1,kdm
            write(25,4100) 'uflxav  ',k,l,  xmin(k+l*kdm),xmax(k+l*kdm)
          enddo
        enddo
        call flush(25)
        endif !1st tile
        call zaiowr3(vflxav,   kdm, iv,.true.,  xmin,xmax, 25, .false.)
        if     (mnproc.eq.1) then
        do l= 0,0
          do k= 1,kdm
            write(25,4100) 'vflxav  ',k,l,  xmin(k+l*kdm),xmax(k+l*kdm)
          enddo
        enddo
        call flush( 25)
        endif !1st tile
        call zaiofl(25)
      endif !flux output
c
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do k=1,kk
          do i=1-nbdy,ii+nbdy
              dpav(i,j,k)=0.
            uflxav(i,j,k)=0.
            vflxav(i,j,k)=0.
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
      return
 4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
      end subroutine restart_out
