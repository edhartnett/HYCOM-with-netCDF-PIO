      subroutine archiv(n, kkout, iyear,iday,ihour, intvl)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
c
      integer   n, kkout, iyear,iday,ihour
      character intvl*3
c
c --- write an archive file.
c
      integer   k,l,nop
      real      coord,xmin,xmax
c
      do l=1,40
        if (flnmarc(l:l).eq.'.') then
          goto 5
        endif
      enddo
      if     (mnproc.eq.1) then
      write (lp,*) 'need decimal point in flnmarc'
      endif
      call xcstop('(flnmarc)')
             stop '(flnmarc)'
 5    continue
c
      if     ((kkout.eq.1 .and. dsurfq.ge.1.0/24.0) .or.
     &        (kkout.gt.1 .and. diagfq.ge.1.0/24.0)     ) then
c ---   indicate the archive date
        write(flnmarc(l+1:l+11),'(i4.4,a1,i3.3,a1,i2.2)') 
     &   iyear,'_',iday,'_',ihour
        l=l+11
      else
c ---   indicate the archive time step
        write(flnmarc(l+1:l+11),'(i6.6,5x)') nstep
        l=l+6
      endif
      nop=13
      call zaiopf(flnmarc(1:l)//'.a', 'new', nop)
      if     (mnproc.eq.1) then
      open (unit=nop,file=flnmarc(1:l)//'.b',status='new')
      write(nop,116) ctitle,iversn,iexpt,yrflag,itdm,jtdm
      call flush(nop)
      endif !1st tile
 116  format (a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & 'field       time step  model day',
     & '  k  dens        min              max')
c
c --- surface fields
c
      coord=0.
c
      call zaiowr(montg(1-nbdy,1-nbdy,1),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'montg1  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(util1,ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'srfhgt  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
c
      call zaiowr(surflx,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'surflx  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(salflx,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'salflx  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
c
      call zaiowr(dpbl,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'bl_dpth ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(dpmixl(1-nbdy,1-nbdy,n),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'mix_dpth',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(tmix,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'tmix    ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(smix,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'smix    ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(thmix,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'thmix   ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(umix,iu,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'umix    ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(vmix,iv,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'vmix    ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      if(icegln) then
        call zaiowr(covice,ip,.true., xmin,xmax, nop, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'covice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
        call zaiowr(thkice,ip,.true., xmin,xmax, nop, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'thkice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
        call zaiowr(temice,ip,.true., xmin,xmax, nop, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'temice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
      endif
c
c --- depth averaged fields
c
      call zaiowr(ubavg(1-nbdy,1-nbdy,n),iu,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'u_btrop ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(vbavg(1-nbdy,1-nbdy,n),iv,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'v_btrop ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
c
c --- layer loop.
c
      do 75 k=1,kkout
      coord=sigma(k)
      call zaiowr(u(1-nbdy,1-nbdy,k,n),iu,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'u-vel.  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(v(1-nbdy,1-nbdy,k,n),iv,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'v-vel.  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(dp(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'thknss  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(temp(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'temp    ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(saln(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'salin   ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(th3d(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'density ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      if(trcout) then
        call zaiowr(tracer(1-nbdy,1-nbdy,k),ip,.true.,
     &              xmin,xmax, nop, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'tracer  ',nstep,time,k,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
      endif
 75   continue
c
 117  format (a8,' =',i11,f11.2,i3,f7.3,1p2e16.7)
c
c --- output time-averaged mass fluxes
c
      if     (kkout.eq.kk) then
        do k=1,kk
          coord=sigma(k)
c
          call zaiowr(diaflx(1-nbdy,1-nbdy,k),ip,.true.,
     &                xmin,xmax, nop, .false.)
          if     (mnproc.eq.1) then
          write (nop,118) 'diafx',intvl,nstep,time,k,coord,xmin,xmax
          call flush(nop)
          endif !1st tile
 118      format (a5,a3,' =',i11,f11.2,i3,f7.3,1p2e16.7)
        enddo
      endif
c
      close (unit=nop)
      call zaiocl(nop)
c
cccc
cccc --- output to line printer
cccc
ccc      call prtmsk(ip,util1,util3,idm,ii,jj,0.,1./(thref*onecm),
ccc     .     'sea surface height (cm)')
ccc      if(mxlkpp) call prtmsk(ip,dpbl,util3,idm,ii,jj,0.,1./onem,
ccc     .     'turb. b.l. depth (m)')
ccc      call prtmsk(ip,dpmixl,util3,idm,ii,jj,0.,1./onem,
ccc     .     'mixed layer depth (m)')
ccc      call prtmsk(ip,tmix,util3,idm,ii,jj,0.,10.,
ccc     .     'mix.layer temp. (.1 deg)')
ccc      call prtmsk(ip,smix,util3,idm,ii,jj,35.,100.,
ccc     .     'mx.lay. salin. (.01 mil)')
cccc!$OMP PARALLEL DO PRIVATE(j,l,i)
ccc!$OMP&         SCHEDULE(STATIC,jblk)
ccc      do j=1-margin,jj+margin
ccc        do l=1,isu(j)
ccc          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
ccc            util1(i,j)=umix(i,j)+ubavg(i,j,n)
ccc          enddo
ccc        enddo
ccc        do l=1,isv(j)
ccc          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
ccc            util2(i,j)=vmix(i,j)+vbavg(i,j,n)
ccc          enddo
ccc        enddo
ccc      enddo
ccc!$OMP END PARALLEL DO
ccc      call prtmsk(iu(2,1),util1(2,1),util3,idm,ii-2,jj,0.,1000.,
ccc     .     'mix.layer u vel. (mm/s)')
ccc      call prtmsk(iv(1,2),util2(1,2),util3,idm,ii,jj-2,0.,1000.,
ccc     .     'mix.layer v vel. (mm/s)')
ccc      call prtmsk(iu(2,1),ubavg(2,1,n),util3,idm,ii-2,jj,0.,1000.,
ccc     .     'barotrop. u vel. (mm/s)')
ccc      call prtmsk(iv(2,1),vbavg(1,2,n),util3,idm,ii,jj-2,0.,1000.,
ccc     .     'barotrop. v vel. (mm/s)')
      return
      end subroutine archiv
