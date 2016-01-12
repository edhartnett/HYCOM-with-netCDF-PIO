      program hycom
      use mod_xc    ! HYCOM communication interface
      use mod_za    ! HYCOM I/O interface
      use mod_pipe  ! HYCOM debugging interface
c
c --- -----------------------------------------
c --- MICOM-based HYbrid Coordinate Ocean Model
c ---               H Y C O M
c ---           v e r s i o n  2.0
c --- -----------------------------------------
c
      implicit none
c
      include 'common_blocks.h'
c
      integer i,ios,j,k,l,m,n,k1n
      real*8  d1,d2,d3,d4,ddsurf,ddiagf,drstrf,
     &        dske,dskea,dsms,dsmsa,dsmt,dsmta,dsum,dsuma,
     &        dtime,dtime0,dbimon,dmonth,dyear,dyear0
      real    sms,smt,sum,smin,smax,
     &        coord,day1,day2,x,x1,time0,timav,cold,utotp,vtotp
      integer nstep0,nod,
     &        lt,ma0,ma1,ma2,ma3,mr0,mr1,mr2,mr3,mnth,
     &        iday,ihour,iyear
*     integer index
      logical linit,diagsv,hisurf,histry,restrt
      character intvl*3
c
      real*8     dsmall,dsmall2,days1,days6
      parameter (dsmall=0.0001d0,dsmall2=0.0002d0)
      parameter (days1=1.d0,days6=6.d0)
c
      character  charac(50)*1
      integer    nflip
c
      data charac/'1','2','3','4','5','6','7','8','9','0',
     &            '1','2','3','4','5','6','7','8','9','0',
     &            '1','2','3','4','5','6','7','8','9','0',
     &            '1','2','3','4','5','6','7','8','9','0',
     &            '1','2','3','4','5','6','7','8','9','0'/
      data nflip/0/
c
      include 'stmt_fns.h'
c
c --- initialize SPMD processsing
      call xcspmd
c
c --- initialize timer names.
c
      call xctmrn(40,'cnuity')
      call xctmrn(41,'tsadvc')
      call xctmrn(42,'momtum')
      call xctmrn(43,'barotp')
      call xctmrn(44,'thermf')
      call xctmrn(45,'ic****')
      call xctmrn(46,'mx****')
      call xctmrn(47,'conv**')
      call xctmrn(48,'diapf*')
      call xctmrn(49,'hybgen')
      call xctmrn(50,'archiv')
      call xctmrn(51,'restrt')
      call xctmrn(52,'overtn')
c
c --- machine-specific initialization
      call machine
c
c --- initialize array i/o.
      call zaiost
c
c --- initiate named-pipe comparison utility
      call pipe_init
c
c --- initialize common variables
c
      call blkdat
      if (dsurfq.ge.1.0) then
        ddsurf = dsurfq
      else
        ddsurf = (baclin/86400.0d0)*
     &           max(1,nint((86400.0d0*dsurfq)/baclin))
      endif
      if (diagfq.ge.1.0) then
        ddiagf = diagfq
      else
        ddiagf = (baclin/86400.0d0)*
     &           max(1,nint((86400.0d0*diagfq)/baclin))
      endif
      if     (rstrfq.eq.0.0) then
        drstrf = rstrfq
      elseif (rstrfq.ge.1.0) then
        drstrf = rstrfq
      else
        drstrf = (baclin/86400.0d0)*
     &           max(1,nint((86400.0d0*rstrfq)/baclin))
      endif
      if (mnproc.eq.1) then
      write(lp,*)
      write(lp,*) 'ddsurf = ',ddsurf,nint((86400.0d0*ddsurf)/baclin)
      write(lp,*) 'ddiagf = ',ddiagf,nint((86400.0d0*ddiagf)/baclin)
      write(lp,*) 'drstrf = ',drstrf,nint((86400.0d0*drstrf)/baclin)
      write(lp,*)
      write (lp,101) thkdff,temdff,veldff,viscos,biharm,diapyc,vertmx
 101  format (' turb. flux parameters:',1p/
     &  ' thkdff,temdff        =',2e9.2/
     &  ' veldff,viscos,biharm =',3e9.2/
     &  ' diapyc,vertmx        =',2e9.2)
      endif !1st tile

c
c --- days in year.
c
      if     (yrflag.eq.0) then
c ---   360 days, starting Jan 16
        dmonth =  30.0d0
        dbimon =  60.0d0
        dyear  = 360.0d0
        dyear0 =   0.0d0
      elseif (yrflag.eq.1) then
c ---   366 days, starting Jan 16
        dmonth =  30.5d0
        dbimon =  61.0d0
        dyear  = 366.0d0
        dyear0 =   0.0d0
      elseif (yrflag.eq.2) then
c ---   366 days, starting Jan 1
c ---   also implies high frequency atmospheric forcing
        dmonth =  30.5d0
        dbimon =  61.0d0
        dyear  = 366.0d0
        dyear0 = -15.0d0+dyear
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
c ---   also implies high frequency atmospheric forcing
        dyear  = 365.25d0
        dmonth = dyear/12.d0
        dbimon = dyear/ 6.d0
        dyear0 = -15.0d0+dyear
      else
        if (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in hycom - unsupported yrflag value'
        write(lp,*)
        call flush(lp)
        endif !1st tile
        call xcstop('(hycom)')
               stop
      endif
c
c --- 'lstep' = number of barotropic time steps per baroclinic time step.
c --- lstep   m u s t   be even.
c
      lstep=baclin/batrop
      lstep=2*((lstep+1)/2)
      dlt=baclin/lstep
      if (mnproc.eq.1) then
      write (lp,'(i4,'' barotropic steps per baroclinic time step'')')
     &  lstep
      endif !1st tile
c
c --- set up parameters defining the geographic environment
c
      call geopar
c
c --- set up forcing functions
c
      if (yrflag.lt.2) then
        call forfuna  ! monthly atmospheric forcing
        call forfunr
      else
        call forfunr
      endif
      watcum=0.
      empcum=0.
c
c --- set minimum salinity for each isopycnic layer
      do k=2,kk
        cold=-3.0
        salmin(k)=(sigma(k)-c1-cold*(c2+cold*(c4+c6*cold)))/
     &      (c3+cold*(c5+c7*cold))
      enddo
c
c --- layer specific volume is defined as (1-theta)*thref
c --- subtract constant 'thbase' from theta to reduce roundoff errors
c
      do k=1,kk
        theta(k)=sigma(k)-thbase
      enddo
c
c --- model is to be integrated from time step 'nstep1' to 'nstep2'
      open(unit=99,file='limits')
      read(99,*) day1,day2,trcrin,trcout
      close(unit=99)
c --- non-positive day1 indicates a new initialization, or
c --- the start of a yrflag==3 case.
      linit =day1.le.0.0 
      day1  =abs(day1)
c
      dtime=day1
      nstep1=nint(dtime*(86400.0d0/baclin))
      dtime=nstep1/(86400.0d0/baclin)
      day1 =dtime
c
      dtime=day2
      nstep2=nint(dtime*(86400.0d0/baclin))
      dtime=nstep2/(86400.0d0/baclin)
      day2 =dtime
c
c --- initialize kpp mixing
      call inikpp
c
      if (linit .and. yrflag.lt.3) then
c
c ---   set up initial conditions
c
        nstep0=nstep1
        dtime0=nstep0/(86400.0d0/baclin)
        time0=dtime0
        delt1=baclin
        if     (clmflg.eq.12) then
          mnth=   1.+nint(mod(dtime0+dyear0,dyear)/dmonth)
        elseif (clmflg.eq.6) then
          mnth=2*(1.+nint(mod(dtime0+dyear0,dyear)/dbimon))-1
        endif
        call inicon(mnth)
c
c ---   output to archive file
c
        m=mod(nstep0  ,2)+1
        n=mod(nstep0+1,2)+1
        nstep=nstep0
        time=dtime0
        call forday(dtime0,yrflag, iyear,iday,ihour)
c
!$OMP   PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=montg(i,j,1)+thref*pbavg(i,j,m)
            enddo
          enddo
          if (isopyc .or. mxlkrt) then
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                umix(i,j)=u(i,j,1,n)
              enddo
            enddo
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                vmix(i,j)=v(i,j,1,n)
              enddo
            enddo
          endif
          do k= 1,kk
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                diaflx(i,j,k)=0.0
              enddo
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
c
        if (mnproc.eq.1) then
        write (intvl,'(i3.3)') 0
        endif !1st tile
        call archiv(n, kk, iyear,iday,ihour, intvl)
c
      else
c
c ---   start from restart file
c
        call restart_in(nstep0,dtime0)
c
        if     (linit) then
c
c ---     start a new calendar-day (yrflag==3) case.
          if (mnproc.eq.1) then
          time0=dtime0
          write (lp,110) time0,day1,nstep1
          endif !1st tile
 110      format (9x,'restart file for day',f8.1,
     &                ', is treated as day',f8.1,'  step',i9)
          nstep0=nstep1
          dtime0=nstep0/(86400.0d0/baclin)
          time0=dtime0
          delt1=baclin+baclin
        else
          nstep0=nint(dtime0*(86400.0d0/baclin))
          dtime0=nstep0/(86400.0d0/baclin)
          time0=dtime0
          delt1=baclin+baclin
          if (mnproc.eq.1) then
          write (lp,111) nstep0,nstep1
          endif !1st tile
 111      format (9x,'time step in restart file -',i9,5x,'wanted -',i9)
          if (nstep0.ne.nstep1) then
            call xcstop('(hycom)')
                   stop
          endif
        endif
      endif
c
c --- set barotp.pot.vort. and layer thickness (incl.bottom pressure) at
c --- u,v points
c
      call dpthuv
c
      call xctilr(dp(     1-nbdy,1-nbdy,1,1),1,kk, nbdy,nbdy, halo_ps)
      call xctilr(dp(     1-nbdy,1-nbdy,1,2),1,kk, nbdy,nbdy, halo_ps)
      call xctilr(dpmixl( 1-nbdy,1-nbdy,  1),1, 2, nbdy,nbdy, halo_ps)
      call xctilr(thkk(   1-nbdy,1-nbdy    ),1, 1, nbdy,nbdy, halo_ps)
      call xctilr(psikk(  1-nbdy,1-nbdy    ),1, 1, nbdy,nbdy, halo_ps)
c
      margin = nbdy
c
      do n=1,2
      m=mod(n,2)+1
c
!$OMP   PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            if     (n.eq.mod(nstep+1,2)+1) then
              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                dpbl(i,j)=dpmixl(i,j,n)
              enddo
            endif
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              p(i,j,1)=0.0
              do k=1,kk
                p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
              enddo
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
c
        call dpudpv(p,depthu,depthv,dpu(1-nbdy,1-nbdy,1,n),
     &                              dpv(1-nbdy,1-nbdy,1,n))
c
        if (.false.) then
c
c ---     ISOPYC TO HYBRID RESTART ONLY
          nstep=nstep1
          call hybgen(m,n)
c
!$OMP     PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do l=1,isp(j)
              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                p(i,j,1)=0.0
                do k=1,kk
                  p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
                enddo
              enddo
            enddo
          enddo
!$OMP     END PARALLEL DO
          call dpudpv(p,depthu,depthv,dpu(1-nbdy,1-nbdy,1,n),
     &                                dpv(1-nbdy,1-nbdy,1,n))
             call pipe_comparall(m,n, 'hybgen, step')
        endif  !isopyc to hybrid restart only
c
      enddo  !n=1,2
c
      nod=14 
      nstep=nstep1
      if (mnproc.eq.1) then
      write (lp,'(/2(a,f8.1),2(a,i9),a/)') 'model starts at day',
     &   time0,', goes to day',time0+day2-day1,'   (steps',nstep1,
     &   ' --',nstep2,')'
      open (unit=nod,file='summary_out',status='unknown')
      write(nod,'(/2(a,f8.1),2(a,i9),a/)') 'model starts at day',
     &   time0,', goes to day',time0+day2-day1,'   (steps',nstep1,
     &   ' --',nstep2,')'
      endif !1st tile
c
      timav=time0
      m=mod(nstep  ,2)+1
      n=mod(nstep+1,2)+1
c
           call pipe_comparall(m,n, 'restrt, step')
c
      if (yrflag.lt.2) then
c
c ---   read in forcing fields for 4 consecutive months
        ma1=1.+mod(dtime0+dyear0,dyear)/dmonth
        ma0=mod(ma1+10,12)+1
        ma2=mod(ma1,   12)+1
        ma3=mod(ma2,   12)+1
        l0=1
        l1=2
        l2=3
        l3=4
        call rdforf(ma0,l0)
        call rdforf(ma1,l1)
        call rdforf(ma2,l2)
        call rdforf(ma3,l3)
      else
c
c ---   initial day of high frequency atmospheric forcing.
c ---   only two fields are used (linear interpolation in time).
        l0=1
        l1=2
        l2=3
        l3=4
        w0=-99.9
        w1=-99.0
        w2=0.0
        w3=0.0
        call forfunh(dtime0)
      endif
c
      if     (clmflg.eq.12) then
c ---   read in relaxation fields for 4 consecutive months
        mr1=1.+mod(dtime0+dyear0,dyear)/dmonth
        mr0=mod(mr1+10,12)+1
        mr2=mod(mr1,   12)+1
        mr3=mod(mr2,   12)+1
        lr0=1
        lr1=2
        lr2=3
        lr3=4
        call rdrlax(mr0,lr0)
        call rdrlax(mr1,lr1)
        call rdrlax(mr2,lr2)
        call rdrlax(mr3,lr3)
      elseif (clmflg.eq.6) then
c ---   read in relaxation fields for 4 consecutive bi-months
        mr1=1.+mod(dtime0+dyear0,dyear)/dbimon
        mr0=mod(mr1+4,6)+1
        mr2=mod(mr1,  6)+1
        mr3=mod(mr2,  6)+1
        lr0=1
        lr1=2
        lr2=3
        lr3=4
        call rdrlax(2*mr0-1,lr0)
        call rdrlax(2*mr1-1,lr1)
        call rdrlax(2*mr2-1,lr2)
        call rdrlax(2*mr3-1,lr3)
      else
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 'error in hycom - unsupported clmflg value'
        call flush(lp)
        endif !1st tile
        call xcstop('(hycom)')
               stop '(hycom)'
      endif
c
c --- report initialization time.
c
      call xctmrp
c
c --- ---------------------
c --- main loop starts here
c --- ---------------------
c
c --- letter 'm' refers to mid-time level (example: dp(i,j,k,m) )
c --- letter 'n' refers to old and new time level
c
 15   continue
      m=mod(nstep  ,2)+1
      n=mod(nstep+1,2)+1
c
      nstep=nstep+1
      dtime=dtime0+(nstep-nstep0)/(86400.0d0/baclin)
      time =dtime
      hisurf=mod(dtime+dsmall,ddsurf).lt.dsmall2 
      histry=mod(dtime+dsmall,ddiagf).lt.dsmall2 
      if     (rstrfq.eq.0.0) then
        restrt=.false.  ! for benchmark cases only
      else
        restrt=mod(dtime+dsmall,drstrf).lt.dsmall2 .or.
     &         nstep.ge.nstep2
      endif
      diagno=mod(dtime+dsmall,ddiagf).lt.dsmall2 .or.
     &       restrt
      if (yrflag.lt.2) then
c
c ---   set weights for quasi-hermite time interpolation for
c ---   monthly atmospheric forcing fields
        x=1.+mod(dtime+dyear0,dyear)/dmonth
c ---   keep quadruplet of forcing functions centered on model time
        if (int(x).ne.ma1) then
          ma1=x
          ma0=mod(ma1+10,12)+1
          ma2=mod(ma1,   12)+1
          ma3=mod(ma2,   12)+1
          lt=l0
          l0=l1
          l1=l2
          l2=l3
          l3=lt
c ---     newest set of forcing functions overwrites set no longer needed
          call rdforf(ma3,l3)
        endif
        x=mod(x,1.)
        x1=1.-x
        w1=x1*(1.+x *(1.-1.5*x ))
        w2=x *(1.+x1*(1.-1.5*x1))
        w0=-.5*x *x1*x1
        w3=-.5*x1*x *x
cdiag   if (mnproc.eq.1) then
cdiag   write (lp,'(i9,'' atmos time interpolation: months'',4i3,
cdiag.    '',  weights '',4f6.3)') nstep,l0,l1,l2,l3,w0,w1,w2,w3
cdiag   endif !1st tile
      else
c
c ---   set weights and fields for high frequency atmospheric forcing.
c ---   only two fields are used (linear interpolation in time).
        call forfunh(dtime)
      endif
c
c --- set weights for quasi-hermite time interpolation for temperature,
c --- salinity and pressure relaxation fields.
      if     (clmflg.eq.12) then
c ---   monthly fields.
        x=1.+mod(dtime+dyear0,dyear)/dmonth
        if (int(x).ne.mr1) then
          mr1=x
          mr0=mod(mr1+10,12)+1
          mr2=mod(mr1,   12)+1
          mr3=mod(mr2,   12)+1
          lt =lr0
          lr0=lr1
          lr1=lr2
          lr2=lr3
          lr3=lt
          call rdrlax(mr3,lr3)
        endif
        x=mod(x,1.)
        x1=1.-x
        wr1=x1*(1.+x *(1.-1.5*x ))
        wr2=x *(1.+x1*(1.-1.5*x1))
        wr0=-.5*x *x1*x1
        wr3=-.5*x1*x *x
      elseif (clmflg.eq.6) then
c ---   bi-monthly fields.
        x=1.+mod(dtime+dyear0,dyear)/dbimon
        if (int(x).ne.mr1) then
          mr1=x
          mr0=mod(mr1+4,6)+1
          mr2=mod(mr1,  6)+1
          mr3=mod(mr2,  6)+1
          lt =lr0
          lr0=lr1
          lr1=lr2
          lr2=lr3
          lr3=lt
          call rdrlax(2*mr3-1,lr3)
        endif
        x=mod(x,1.)
        x1=1.-x
        wr1=x1*(1.+x *(1.-1.5*x ))
        wr2=x *(1.+x1*(1.-1.5*x1))
        wr0=-.5*x *x1*x1
        wr3=-.5*x1*x *x
      endif
cdiag if (mnproc.eq.1) then
cdiag write (lp,'(i9,'' relax time interpolation: months'',4i3,
cdiag.  '',  weights '',4f6.3)') nstep,lr0,lr1,lr2,lr3,wr0,wr1,wr2,wr3
cdiag endif !1st tile
c
         call pipe_comparall(m,n, 'ENTER , step')
         call xctmr0(40)
      call cnuity(m,n)
         call xctmr1(40)
         call pipe_comparall(m,n, 'cnuity, step')
         call xctmr0(41)
      call tsadvc(m,n)
         call xctmr1(41)
         call pipe_comparall(m,n, 'tsadvc, step')
         call xctmr0(42)
      call momtum(m,n)
         call xctmr1(42)
         call pipe_comparall(m,n, 'momtum, step')
         call xctmr0(43)
      call barotp(m,n)
         call xctmr1(43)
         call pipe_comparall(m,n, 'barotp, step')
         call xctmr0(44)
      call thermf(m,n)
         call xctmr1(44)
         call pipe_comparall(m,n, 'thermf, step')
      if (icegln) then
            call xctmr0(45)
         call icloan(m,n)
            call xctmr1(45)
            call pipe_comparall(m,n, 'icloan, step')
      else
            call xctmr0(45)
            call xctmr1(45)
      endif
      if (hybrid) then
         diagsv = diagno
         diagno = diagno .or. nstep.eq.nstep0+1 .or.
     &            histry .or. hisurf .or.
     &            mod(dtime+dsmall,days6).lt.dsmall2
         if (mxlkpp) then
               call xctmr0(46)
            call mxkpp(m,n)
               call xctmr1(46)
               call pipe_comparall(m,n, 'mxkpp , step')
         elseif (mxlkta) then
               call xctmr0(46)
            call mxkrta(m,n)
               call xctmr1(46)
               call pipe_comparall(m,n, 'mxkrta, step')
         elseif (mxlktb) then
               call xctmr0(46)
            call mxkrtb(m,n)
               call xctmr1(46)
               call pipe_comparall(m,n, 'mxkrtb, step')
         else
               call xctmr0(46)
               call xctmr1(46)
         endif
         diagno = diagsv
         if (mxlkpp) then
               call xctmr0(47)
               call xctmr1(47)
               call xctmr0(48)
               call xctmr1(48)
         else
               call xctmr0(47)
            call convch(m,n)
               call xctmr1(47)
               call pipe_comparall(m,n, 'convch, step')
           if (dypflg.eq.1) then
                 call xctmr0(48)
              call diapf1(m,n)
                 call xctmr1(48)
                 call pipe_comparall(m,n, 'diapf1, step')
           elseif (dypflg.eq.2) then
                 call xctmr0(48)
              call diapf2(m,n)
                 call xctmr1(48)
                 call pipe_comparall(m,n, 'diapf2, step')
           else
                 call xctmr0(48)
                 call xctmr1(48)
           endif
        endif
           call xctmr0(49)
        call hybgen(m,n)
           call xctmr1(49)
           call pipe_comparall(m,n, 'hybgen, step')
      else  ! isopyc
            call xctmr0(46)
         call mxkrtm(m,n)
            call xctmr1(46)
            call pipe_comparall(m,n, 'mxkrtm, step')
            call xctmr0(47)
         call convcm(m,n)
            call xctmr1(47)
            call pipe_comparall(m,n, 'convcm, step')
            call xctmr0(48)
         call diapf3(m,n)
            call xctmr1(48)
            call pipe_comparall(m,n, 'diapf3, step')
            call xctmr0(49)
            call xctmr1(49)
      endif
c
c ---------------------------------------------------------------------------
c
c --- output and diagnostic calculations
c
c ---------------------------------------------------------------------------
c
      if (diagno .or. histry .or. hisurf .or.
     &    nstep.eq.nstep0+1 .or.
     &    mod(dtime+dsmall,days1).lt.dsmall2) then  ! at least daily
c
        call forday(dtime,yrflag, iyear,iday,ihour)
c
c ---   diagnose mean sea surface height

        smin= huge
        smax=-huge
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&              REDUCTION(MIN:smin) REDUCTION(MAX:smax)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c ---         compute sea surface height
              util1(i,j)=montg(i,j,1)+thref*pbavg(i,j,m)
              util3(i,j)=util1(i,j)  *scp2(i,j)
              util4(i,j)=montg(i,j,1)*scp2(i,j)
              smin=min(smin,util1(i,j))
              smax=max(smax,util1(i,j))
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xcminr(smin)
        call xcmaxr(smax)
        call xcsum( dsum, util3,ip)
        call xcsum( dsmt, util4,ip)
        sum=dsum
        smt=dsmt
        if (mnproc.eq.1) then
        write (lp,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &              '' mean sea srf.hgt. (mm):'',f8.2,
     &              ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &  nstep,iyear,iday,
     &  sum/(area*thref*onemm),smin/(thref*onemm),smax/(thref*onemm)
*       write (lp,'(i9,'' ('',i4,''/'',i3.3,'')'',
*    .              '' mean Montg. pot.  (mm):'',f8.2)')
*    .  nstep,iyear,iday,
*    .  smt/(area*thref*onemm)
        call flush(lp)
        write(nod,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &              '' mean sea srf.hgt. (mm):'',f8.2,
     &              ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &  nstep,iyear,iday,
     &  sum/(area*thref*onemm),smin/(thref*onemm),smax/(thref*onemm)
*       write(nod,'(i9,'' ('',i4,''/'',i3.3,'')'',
*    .              '' mean Montg. pot.  (mm):'',f8.2)')
*    .  nstep,iyear,iday,
*    .  smt/(area*thref*onemm)
        call flush(nod)
        endif !1st tile
c ---   NaN detection.
        if     (smin.ne.smin .or. sum.ne.sum .or.
     &          smax.ne.smax .or. smt.ne.smt      ) then
          if (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error - NaN detected'
          write(lp,*)
          call flush(lp)
          endif !1st tile
          call xcstop('(hycom)')
                 stop
        endif
*     else
*       if (mnproc.eq.1) then
*       write (lp,'('' time step ='',i9)') nstep
*       call flush(lp)
*       endif !1st tile
      endif
c
c --- diagnose heat flux, ice, layer thickness and temperature,
c --- mean temperature and mean kinetic energy
c --- note that mixed-layer fields must be switched on in mxkpp
      if (diagno .or.
     &    nstep.eq.nstep0+1 .or.
     &    mod(dtime+dsmall,days6).lt.dsmall2) then  ! at least every 6 days
c
        call forday(dtime,yrflag, iyear,iday,ihour)
c
        if (thermo) then
!$OMP     PARALLEL DO PRIVATE(j,l,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=surflx(i,j)*scp2(i,j)
                if (mxlkpp .or. (mxlkrt .and. pensol)) then
                  util2(i,j)=sswflx(i,j)*scp2(i,j)*
     &                      (1.-redfac(jerlov(i,j)))*
     &                      exp(-dpmixl(i,j,n)/betabl(jerlov(i,j)))
                endif
              enddo
            enddo
          enddo
!$OMP     END PARALLEL DO
          call xcsum(dsum, util1,ip)
          call xcsum(dsmt, util2,ip)
          sum= dsum/area
          smt=(dsum-dsmt)/area
          if (mnproc.eq.1) then
          write (lp,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &        '' mean heat flux (w/m/m):'',f8.2,
     &                           ''  ml:'',f8.2)')
     &      nstep,iyear,iday,
     &      sum,smt
          call flush(lp)
          write(nod,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &        '' mean heat flux (w/m/m):'',f8.2,
     &                           ''  ml:'',f8.2)')
     &      nstep,iyear,iday,
     &      sum,smt
          call flush(nod)
          endif !1st tile
        endif
c
        if (icegln) then
!$OMP     PARALLEL DO PRIVATE(j,l,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                if     (covice(i,j).ne.0.0) then
                  util1(i,j)=                        scp2(i,j)
                  util2(i,j)=covice(i,j)*            scp2(i,j)
                  util3(i,j)=            thkice(i,j)*scp2(i,j)
                  util4(i,j)=covice(i,j)*temice(i,j)*scp2(i,j)
                else
                  util1(i,j)=0.0
                  util2(i,j)=0.0
                  util3(i,j)=0.0
                  util4(i,j)=0.0
                endif
              enddo
            enddo
          enddo
!$OMP     END PARALLEL DO
          call xcsum(d1, util1,ip)
          call xcsum(d2, util2,ip)
          call xcsum(d3, util3,ip)
          call xcsum(d4, util4,ip)
          if     (d2.gt.0.0d0) then
            sum=d3/d1
            smt=d4/d2
            sms=d2/area * 100.0
          else
            sum=0.0 
            smt=0.0 
            sms=0.0 
          endif
          if (mnproc.eq.1) then
          write (lp,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &              '' mean      ice thk. (m):'',f8.2,
     &                               ''  temp:'',f6.2,
     &                               ''  pcen:'',f7.3)')
     &      nstep,iyear,iday,
     &      sum,smt,sms
          call flush(lp)
          write(nod,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &              '' mean      ice thk. (m):'',f8.2,
     &                               ''  temp:'',f6.2,
     &                               ''  pcen:'',f7.3)')
     &      nstep,iyear,iday,
     &      sum,smt,sms
          call flush(nod)
          endif !1st tile
        endif  ! icegln
c
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=dpmixl(i,j,n)*scp2(i,j)
              util2(i,j)=dpmixl(i,j,n)*scp2(i,j)*tmix(i,j)
              util3(i,j)=dpmixl(i,j,n)*scp2(i,j)*smix(i,j)
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xcsum(dsum, util1,ip)
        call xcsum(dsmt, util2,ip)
        call xcsum(dsms, util3,ip)
        if     (dsum.ne.0.0d0) then
          sum=dsum/(area*onem)
          smt=dsmt/dsum
          sms=dsms/dsum
        else
          sum=0.0
          smt=0.0
          sms=0.0
        endif
        if (mnproc.eq.1) then
        write (lp,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &              '' mean mixlayer thk. (m):'',f8.2,
     &                               ''  temp:'',f6.2,
     &                               ''  saln:'',f7.3)')
     &      nstep,iyear,iday,
     &      sum,smt,sms
        call flush(lp)
        write(nod,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &              '' mean mixlayer thk. (m):'',f8.2,
     &                               ''  temp:'',f6.2,
     &                               ''  saln:'',f7.3)')
     &      nstep,iyear,iday,
     &      sum,smt,sms
        call flush(nod)
        endif !1st tile
c
        if     (relaxf) then
c
c ---     mean surface climatology.
c
!$OMP     PARALLEL DO PRIVATE(j,l,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=scp2(i,j)*
     &                (twall(i,j,1,lr0)*wr0+twall(i,j,1,lr1)*wr1
     &                +twall(i,j,1,lr2)*wr2+twall(i,j,1,lr3)*wr3)
                util2(i,j)=scp2(i,j)*
     &                (swall(i,j,1,lr0)*wr0+swall(i,j,1,lr1)*wr1
     &                +swall(i,j,1,lr2)*wr2+swall(i,j,1,lr3)*wr3)
              enddo
            enddo
          enddo
!$OMP     END PARALLEL DO
          call xcsum(dsmt, util1,ip)
          call xcsum(dsms, util2,ip)
          smt=dsmt/area
          sms=dsms/area
          if (mnproc.eq.1) then
          write (lp,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &                '' mean surfclim thk. (m):'',f8.2,
     &                                 ''  temp:'',f6.2,
     &                                 ''  saln:'',f7.3)')
     &        nstep,iyear,iday,
     &        thkmin,smt,sms
          call flush(lp)
          write(nod,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &                '' mean surfclim thk. (m):'',f8.2,
     &                                 ''  temp:'',f6.2,
     &                                 ''  saln:'',f7.3)')
     &        nstep,iyear,iday,
     &        thkmin,smt,sms
          call flush(nod)
          endif !1st tile
        endif
c
        call xctilr(u(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
        call xctilr(v(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_vv)
        call xctilr(ubavg(  1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_uv)
        call xctilr(vbavg(  1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_vv)
c
        dsuma=0.0d0
        dsmta=0.0d0
        dsmsa=0.0d0
        dskea=0.0d0
        do k=1,kk
!$OMP     PARALLEL DO PRIVATE(j,l,i,utotp,vtotp)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                utotp=0.5*( u(i,  j,k,n)+ubavg(i,  j,n) +
     &                      u(i+1,j,k,n)+ubavg(i+1,j,n)  )
                vtotp=0.5*( v(i,j,  k,n)+vbavg(i,j,  n) +
     &                      v(i,j+1,k,n)+vbavg(i,j+1,n)  )
                util1(i,j)=dp(i,j,k,n)*scp2(i,j)
                util2(i,j)=dp(i,j,k,n)*scp2(i,j)*temp(i,j,k,n)
                util3(i,j)=dp(i,j,k,n)*scp2(i,j)*saln(i,j,k,n)
                util4(i,j)=dp(i,j,k,n)*scp2(i,j)*
     &                    0.5*(1000.0+th3d(i,j,k,n)+thbase)*
     &                        (utotp**2+vtotp**2)
              enddo
            enddo
          enddo
!$OMP     END PARALLEL DO
          call xcsum(dsum, util1,ip)
          call xcsum(dsmt, util2,ip)
          call xcsum(dsms, util3,ip)
          call xcsum(dske, util4,ip)
          dsuma=dsuma+dsum
          dsmta=dsmta+dsmt
          dsmsa=dsmsa+dsms
          dskea=dskea+dske
          if     (dsum.ne.0.0d0) then
            sum=dsum/(area*onem)
            smt=dsmt/dsum
            sms=dsms/dsum
          else
            sum=0.0
            smt=0.0
            sms=0.0
          endif
          if (mnproc.eq.1) then
          write (lp,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &                '' mean layer '',i2,'' thk. (m):'',f8.2,
     &                                       ''  temp:'',f6.2,
     &                                       ''  saln:'',f7.3)')
     &        nstep,iyear,iday,
     &        k,sum,smt,sms
          call flush(lp)
          write(nod,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &                '' mean layer '',i2,'' thk. (m):'',f8.2,
     &                                       ''  temp:'',f6.2,
     &                                       ''  saln:'',f7.3)')
     &        nstep,iyear,iday,
     &        k,sum,smt,sms
          call flush(nod)
          endif !1st tile
        enddo
        sum=dskea/(area*onem)
        smt=dsmta/dsuma
        sms=dsmsa/dsuma
        if (mnproc.eq.1) then
        write (lp,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &                '' region-wide mean K.E. :'',f8.2,
     &                                 ''  temp:'',f6.2,
     &                                 ''  saln:'',f7.3)')
     &      nstep,iyear,iday,
     &        sum,smt,sms
        call flush(lp)
        write(nod,'(i9,'' ('',i4,''/'',i3.3,'')'',
     &                '' region-wide mean K.E. :'',f8.2,
     &                                 ''  temp:'',f6.2,
     &                                 ''  saln:'',f7.3)')
     &      nstep,iyear,iday,
     &        sum,smt,sms
        call flush(nod)
        endif !1st tile
      endif
c
c --- diagnose meridional overturning and heat flux
****  if     (mod(dtime+dsmall,dmonth).lt.dsmall2) then
****    call xctmr0(52)
****    call overtn(dtime,dyear)
****    call xctmr1(52)
****  elseif (nstep.ge.nstep2) then
****    call xctmr0(52)
****    call overtn(dtime,dyear)
****    call xctmr1(52)
****  endif
c
      if (histry .or.hisurf) then
        call xctmr0(50)
c
c ---   output to archive file
c
        call forday(dtime,yrflag, iyear,iday,ihour)
c
!$OMP   PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j)=montg(i,j,1)+thref*pbavg(i,j,m)
            enddo
          enddo
          if (isopyc .or. mxlkrt) then
            do l=1,isu(j)
              do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
                umix(i,j)=u(i,j,1,n)
              enddo
            enddo
            do l=1,isv(j)
              do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
                vmix(i,j)=v(i,j,1,n)
              enddo
            enddo
          endif
          if (histry) then
            do k= 1,kk
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c ---             convert diapycnal thickness changes into 
c ---             actual interface fluxes
                  if (k.gt.1) then
                    diaflx(i,j,k)=diaflx(i,j,k)/(2.*onem) +
     &                            diaflx(i,j,k-1)
                  else
                    diaflx(i,j,k)=diaflx(i,j,k)/(2.*onem)
                  endif
                enddo
              enddo
            enddo
          endif
        enddo
!$OMP   END PARALLEL DO
c
        if (mnproc.eq.1) then
        write (intvl,'(i3.3)') int(dtime-timav+dsmall)
        endif !1st tile
        if (hisurf .and. .not. histry) then
          call archiv(n, 1,  iyear,iday,ihour, intvl)
        else
          call archiv(n, kk, iyear,iday,ihour, intvl)
        endif
c
        if (histry) then
!$OMP     PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1,jj
            do k= 1,kk
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  diaflx(i,j,k)=0.0
                enddo
              enddo
            enddo
          enddo
!$OMP     END PARALLEL DO
c
          timav=time
        endif
        call xctmr1(50)
      endif  ! histry.or.hisurf
c
      if (restrt) then
        call xctmr0(51)
c
c ---   output to restart and flux statitics files
c
        call forday(dtime,yrflag, iyear,iday,ihour)
        if (mnproc.eq.1) then
        write (lp,100) nstep,iyear,iday,ihour
        call flush(lp)
        endif !1st tile
 100    format (' time step',i9,
     &          9x,'y e a r',i6.4,
     &          9x,'d a y',  i5.3,
     &          9x,'h o u r',i4.2)
c
        if (mxlkpp) then
!$OMP     PARALLEL DO PRIVATE(j,l,i)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                dpmixl(i,j,1) = dpbl(i,j)
                dpmixl(i,j,2) = dpbl(i,j)
              enddo
            enddo
          enddo
!$OMP     END PARALLEL DO
        endif
c
        call restart_out(nstep,dtime)
c
        call xctmr1(51)
      endif  ! restrt
c
      if     (histry .or. hisurf .or. restrt) then
        if (mnproc.eq.1) then
        write (lp,105) nstep,dtime
        call flush(lp)
        endif !1st tile
 105    format (' step',i9,' day',f9.2,' -- archiving completed --')
      endif
c
      if (nstep.ge.nstep2) then
c
c ---   end of the run.
        if (mnproc.eq.1) then
        write(nod,'(a)') 'normal stop'
        call flush(nod)
        endif !1st tile
        call xcstop('(normal)')
               stop '(normal)'
      endif
      delt1=baclin+baclin
      go to 15
      end program hycom
c
c
c> Revision history:
c>
c> May  1997 - removed statement "theta(1)=-thbase" after loop 14
c> June 1997 - added loop 60 to fix bug in vertical summation of -diaflx-
c> Oct. 1999 - option for krt or kpp mixed layer model - convec and diapfl
c>             not called for kpp mixing model
c> Oct. 1999 - dpbl (boundary layer thickness) is output in addition to
c>             dpmixl when the kpp mixing model is selected
c> May  2000 - conversion to SI units
c> Aug. 2000 - added isopycnic (MICOM) vertical coordinate option
c> Oct. 2000 - added option for high frequency atmospheric forcing
c> Nov. 2000 - archive time stamp is either time step or YYYY_DDD_HH
