      subroutine thermf(m,n)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- ---------------
c --- thermal forcing
c --- ---------------
c
      integer i,j,k,ktr,nm,l
      real    pwl,q
      real*8  t1mean,s1mean,tmean,smean,pmean,rmean,
     &        rareac,runsec,secpyr
      real*8  d1,d2,d3,d4
c
      real    pwij(kk+1),trwij(kk,ntracr),
     &        prij(kk+1),trcij(kk,ntracr)
c
      real*8  tmean0,smean0,rmean0
      save    tmean0,smean0,rmean0
c
      include 'stmt_fns.h'
c
      margin = 0  ! no horizontal derivatives
c
!$OMP PARALLEL DO PRIVATE(j,k,l,i,ktr)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do k=1,kk
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
c --- ----------------------------
c --- thermal forcing at nestwalls
c --- ----------------------------
c
      if (nestfq.ne.0.0 .and. nstep.gt.2) then
c
!$OMP PARALLEL DO PRIVATE(j,i,k,pwl)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (ip(i,j).eq.1 .and. rmunp(i,j).ne.0.0) then
            k=1
            saln(i,j,k,n)=saln(i,j,k,n)+delt1*rmunp(i,j)*
     &         ((snest(i,j,k,ln0)*wn0+snest(i,j,k,ln1)*wn1)
     &          - saln(i,j,k,n))
            temp(i,j,k,n)=temp(i,j,k,n)+delt1*rmunp(i,j)*
     &         ((tnest(i,j,k,ln0)*wn0+tnest(i,j,k,ln1)*wn1)
     &          - temp(i,j,k,n))
            th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
c
            if     (hybrid) then
              do k=kk,2,-1
                pwl=pnest(i,j,k,ln0)*wn0+pnest(i,j,k,ln1)*wn1
                if     (pwl.gt.p(i,j,kk+1)-tencm) then
                  pwl=p(i,j,kk+1)
                endif
                p(i,j,k)=min(p(i,j,k+1),
     &                       p(i,j,k)+delt1*rmunp(i,j)*(pwl-p(i,j,k)))
                dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
c
                if     (pwl.lt.p(i,j,kk+1)) then
                  saln(i,j,k,n)=saln(i,j,k,n)+delt1*rmunp(i,j)*
     &                 ((snest(i,j,k,ln0)*wn0+snest(i,j,k,ln1)*wn1)
     &                  - saln(i,j,k,n))
                  if     (k.le.nhybrd) then
                    temp(i,j,k,n)=temp(i,j,k,n)+delt1*rmunp(i,j)*
     &                 ((tnest(i,j,k,ln0)*wn0+tnest(i,j,k,ln1)*wn1)
     &                  - temp(i,j,k,n))
                    th3d(i,j,k,n)=sig(temp(i,j,k,n),
     &                                saln(i,j,k,n))-thbase
                  else
                    th3d(i,j,k,n)=       theta(i,j,k)
                    temp(i,j,k,n)=tofsig(theta(i,j,k)+thbase,
     &                                   saln(i,j,k,n))
                  endif
                endif
              enddo  !k
              dp(i,j,1,n)=p(i,j,2)-p(i,j,1)
            else  ! isopyc
              do k=kk,2,-1
                saln(i,j,k,n)=saln(i,j,k,n)+delt1*rmunp(i,j)*
     &             ((snest(i,j,k,ln0)*wn0+snest(i,j,k,ln1)*wn1)
     &              - saln(i,j,k,n))
                temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
                if (k.ge.3) then
                  pwl=pnest(i,j,k,ln0)*wn0+pnest(i,j,k,ln1)*wn1
                  pwl=max(p(i,j,2),pwl)
                  if     (pwl.gt.p(i,j,kk+1)-tencm) then
                    pwl=p(i,j,kk+1)
                  endif
                  p(i,j,k)=min(p(i,j,k+1),
     &                         p(i,j,k)+delt1*rmunp(i,j)*(pwl-p(i,j,k)))
                endif
                dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
              enddo  !k
            endif  ! hybrid:isopyc
c
c ---       minimal tracer support (non-negative in buffer zone).
            do ktr= 1,ntracr
              tracer(i,j,k,n,ktr)=max(tracer(i,j,k,n,ktr),0.0)
            enddo
          endif  !ip.eq.1 .and. rmunp.ne.0.0
c
          if (iu(i,j).eq.1 .and.
     &        max(rmunv(i,j),rmunv(i-1,j)).ne.0.0) then
            do k= 1,kk
              pwl=u(i,j,k,n)
              u(i,j,k,n)=u(i,j,k,n)+delt1*max(rmunv(i,j),rmunv(i-1,j))*
     &           ((unest(i,j,k,ln0)*wn0+unest(i,j,k,ln1)*wn1)
     &               - u(i,j,k,n))
            enddo  !k
          endif  !iu.eq.1 .and. rmunv.ne.0.0
c
          if (iv(i,j).eq.1 .and.
     &        max(rmunv(i,j),rmunv(i,j-1)).ne.0.0) then
            do k= 1,kk
              pwl=v(i,j,k,n)
              v(i,j,k,n)=v(i,j,k,n)+delt1*max(rmunv(i,j),rmunv(i,j-1))*
     &           ((vnest(i,j,k,ln0)*wn0+vnest(i,j,k,ln1)*wn1)
     &               - v(i,j,k,n))
            enddo  !k
          endif  !iu.eq.1 .and. rmunv.ne.0.0
        enddo  !i
      enddo  !j
!$OMP END PARALLEL DO
c
      endif  !  nestfq.ne.0.0
c
c --- ----------------------------
c --- thermal forcing at sidewalls
c --- ----------------------------
c
      if (relax .and. nstep.gt.2) then
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,pwl)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 53 j=1-margin,jj+margin
      do 53 l=1,isp(j)
      do 53 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
        if (rmu(i,j).ne.0.0) then
          k=1
          saln(i,j,k,n)=saln(i,j,k,n)+delt1*rmu(i,j)*
     &       (( swall(i,j,k,lc0)*wc0+swall(i,j,k,lc1)*wc1
     &         +swall(i,j,k,lc2)*wc2+swall(i,j,k,lc3)*wc3)
     &        - saln(i,j,k,n))
          if     (lwflag.eq.2 .or. sstflg.gt.2   .or.
     &            icmflg.eq.2 .or. ticegr.eq.0.0     ) then
c ---       use seatmp, since it is the best available SST
            temp(i,j,k,n)=temp(i,j,k,n)+delt1*rmu(i,j)*
     &         (( seatmp(i,j,l0)*w0+seatmp(i,j,l1)*w1
     &           +seatmp(i,j,l2)*w2+seatmp(i,j,l3)*w3)
     &          - temp(i,j,k,n))
          else
            temp(i,j,k,n)=temp(i,j,k,n)+delt1*rmu(i,j)*
     &         (( twall(i,j,k,lc0)*wc0+twall(i,j,k,lc1)*wc1
     &           +twall(i,j,k,lc2)*wc2+twall(i,j,k,lc3)*wc3)
     &          - temp(i,j,k,n))
          endif
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
c
          if     (hybrid) then
            do k=kk,2,-1
              pwl=pwall(i,j,k,lc0)*wc0+pwall(i,j,k,lc1)*wc1
     &           +pwall(i,j,k,lc2)*wc2+pwall(i,j,k,lc3)*wc3
              if     (pwl.gt.p(i,j,kk+1)-tencm) then
                pwl=p(i,j,kk+1)
              endif
              p(i,j,k)=min(p(i,j,k+1),
     &                     p(i,j,k)+delt1*rmu(i,j)*(pwl-p(i,j,k)))
              dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
c
              if     (pwl.lt.p(i,j,kk+1)) then
                saln(i,j,k,n)=saln(i,j,k,n)+delt1*rmu(i,j)*
     &             ((swall(i,j,k,lc0)*wc0+swall(i,j,k,lc1)*wc1
     &              +swall(i,j,k,lc2)*wc2+swall(i,j,k,lc3)*wc3)
     &              - saln(i,j,k,n))
                if     (k.le.nhybrd) then
                  temp(i,j,k,n)=temp(i,j,k,n)+delt1*rmu(i,j)*
     &               ((twall(i,j,k,lc0)*wc0+twall(i,j,k,lc1)*wc1
     &                +twall(i,j,k,lc2)*wc2+twall(i,j,k,lc3)*wc3)
     &                - temp(i,j,k,n))
                  th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
                else
                  th3d(i,j,k,n)=       theta(i,j,k)
                  temp(i,j,k,n)=tofsig(theta(i,j,k)+thbase,
     &                                 saln(i,j,k,n))
                endif !hybrid:else
              endif !pwl.lt.p(i,j,kk+1)
            enddo !k
            dp(i,j,1,n)=p(i,j,2)-p(i,j,1)
          else  ! isopyc
            do k=kk,2,-1
              saln(i,j,k,n)=saln(i,j,k,n)+delt1*rmu(i,j)*
     &           ((swall(i,j,k,lc0)*wc0+swall(i,j,k,lc1)*wc1
     &            +swall(i,j,k,lc2)*wc2+swall(i,j,k,lc3)*wc3)
     &            - saln(i,j,k,n))
              temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
              if (k.ge.3) then
                pwl=pwall(i,j,k,lc0)*wc0+pwall(i,j,k,lc1)*wc1
     &             +pwall(i,j,k,lc2)*wc2+pwall(i,j,k,lc3)*wc3
                pwl=max(p(i,j,2),pwl)
                if     (pwl.gt.p(i,j,kk+1)-tencm) then
                  pwl=p(i,j,kk+1)
                endif
                p(i,j,k)=min(p(i,j,k+1),
     &                       p(i,j,k)+delt1*rmu(i,j)*(pwl-p(i,j,k)))
              endif !k.ge.3
              dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
            enddo !k
          endif !hybrid:isopyc
        endif !rmu(i,j).ne.0.0
 53   continue
!$OMP END PARALLEL DO
c
      endif  !  relax = .true.
c
c --- ----------------------------
c --- tracer forcing at sidewalls
c --- ----------------------------
c
      if (trcrlx .and. nstep.gt.2) then
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,ktr,pwij,trwij,prij,trcij)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              if     (rmutra(i,j).ne.0.0) then !at least one mask is non-zero
                prij(1)=0.0
                do k=1,kk
                  prij(k+1) =  prij(k)+dp(i,j,k,n)
                  pwij(k)   =  pwall(i,j,k,lc0)*wc0
     &                        +pwall(i,j,k,lc1)*wc1
     &                        +pwall(i,j,k,lc2)*wc2
     &                        +pwall(i,j,k,lc3)*wc3
                  do ktr= 1,ntracr
                    trwij(k,ktr) =  trwall(i,j,k,lc0,ktr)*wc0
     &                             +trwall(i,j,k,lc1,ktr)*wc1
     &                             +trwall(i,j,k,lc2,ktr)*wc2
     &                             +trwall(i,j,k,lc3,ktr)*wc3
                  enddo !ktr
                enddo !k
                pwij(kk+1)=prij(kk+1)
*               call plctrc(trwij,pwij,kk,ntracr,
*    &                      trcij,prij,kk        )
                call plmtrc(trwij,pwij,kk,ntracr,
     &                      trcij,prij,kk        )
                do ktr= 1,ntracr
                  if     (rmutr(i,j,ktr).ne.0.0) then
                    do k=1,kk
                      tracer(i,j,k,n,ktr) = tracer(i,j,k,n,ktr)+
     &                delt1*rmutr(i,j,ktr)*(trcij(k,ktr)-
     &                                      tracer(i,j,k,n,ktr))
                    enddo !k
                  endif !rmutr.ktr.ne.0.0
                enddo !ktr
              endif !rmutra.ne.0.0
            enddo !i
          enddo !l
        enddo !j
!$OMP   END PARALLEL DO
c
      endif  !  trcrlx = .true.
c
c --- --------------------------------
c --- thermal forcing of ocean surface
c --- --------------------------------
c
      if (thermo .or. sstflg.gt.0 .or. srelax) then
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call thermfj(m,n, j)
      enddo
!$OMP END PARALLEL DO
      call xcsum(d1, util1,ip)
      call xcsum(d2, util2,ip)
      call xcsum(d3, util3,ip)
      call xcsum(d4, util4,ip)
      watcum=watcum+d1
      empcum=empcum+d2
      t1mean=d3
      s1mean=d4
c
c --- smooth surface fluxes?
c
      if     (flxsmo) then
        call psmooth_ice(surflx, 0)
        call psmooth_ice(salflx, 0)
      endif
c
      if (nstep.eq.nstep1+1 .or. diagno) then
!$OMP   PARALLEL DO PRIVATE(j,k,nm,l,i)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          k=1
            nm=n
            if (nstep.eq.nstep1+1) nm=m
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=               dp(i,j,k,nm)*scp2(i,j)
                util2(i,j)=temp(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
                util3(i,j)=saln(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
                util4(i,j)=th3d(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
              enddo
            enddo
          do k=2,kk
            nm=n
            if (nstep.eq.nstep1+1) nm=m
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                util1(i,j)=util1(i,j)+    dp(i,j,k,nm)*scp2(i,j)
                util2(i,j)=util2(i,j)+
     &                     temp(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
                util3(i,j)=util3(i,j)+
     &                     saln(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
                util4(i,j)=util4(i,j)+
     &                     th3d(i,j,k,nm)*dp(i,j,k,nm)*scp2(i,j)
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
        if     (mnproc.eq.1) then
        write (lp,'(i9,a,3f9.3)') 
     &    nstep,' mean basin temp, saln, dens ',
     &    tmean,smean,rmean+thbase
        endif !1st tile
        if     (nstep.eq.nstep1+1) then
c
c ---     save initial basin means.
          tmean0=tmean
          smean0=smean
          rmean0=rmean
        else
c
c ---     diagnostic printout of fluxes.
          rareac=1.0/(area*(nstep-nstep1))
          runsec=   baclin*(nstep-nstep1)
          if      (yrflag.eq.0) then
            secpyr=360.00d0*86400.0d0
          elseif (yrflag.lt.3) then
            secpyr=366.00d0*86400.0d0
          elseif (yrflag.ge.3) then
            secpyr=365.25d0*86400.0d0
          endif
          if     (mnproc.eq.1) then
          write (lp,'(i9,a,2f9.3)') 
     &     nstep,' mean surface temp and saln  ',
     &     t1mean/area,s1mean/area
          write (lp,'(i9,a,2f9.3,a)') 
     &     nstep,' energy residual (atmos,tot) ',
     &     watcum*rareac,
     &     (tmean-tmean0)*(spcifh*avgbot*qthref)/runsec,
     &    ' (w/m**2)'
c ---     note that empcum is now salflx cum.
          write (lp,'(i9,a,2f9.3,a)')
     &     nstep,'  e - p residual (atmos,tot) ',
     &     empcum*(thref/saln0)*rareac*100.0*secpyr,
     &     (smean-smean0)/(saln0*runsec)*avgbot*100.0*secpyr,
     &    ' (cm/year)'
          write (lp,'(i9,a,2f9.3)') 
     &     nstep,' temp drift per century      ',
     &     (watcum*rareac/(spcifh*avgbot*qthref))*(secpyr*100.0d0),
     &     (tmean-tmean0)*(secpyr*100.0d0)/runsec
          write (lp,'(i9,a,2f9.3)') 
     &     nstep,' saln drift per century      ',
     &     (empcum*rareac/(       avgbot*qthref))*(secpyr*100.0d0),
     &     (smean-smean0)*(secpyr*100.0d0)/runsec
          write (lp,'(i9,a,9x,f9.3)') 
     &     nstep,' dens drift per century      ',
     &     (rmean-rmean0)*(secpyr*100.0d0)/runsec
          endif !1st tile
          call xcsync(flush_lp)
        endif
      endif
c
      endif   !  thermo .or.  sstflg.gt.0 .or. srelax
c
      return
      end subroutine thermf
c
      subroutine thermfj(m,n, j)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, j
c
      real    radfl,swfl,sstrlx,wind,airt,vpmx,prcp,xtau,ytau,
     &        evap,emnp,snsibl,dsgdt,tmn,smn,rmus,rmut
      real    cd0,clh,cl0,cl1,csh,
     &        rair,slat,ssen,tdif,tsur,wsph,
     &        tamts,q,qva,va
      integer i,l
c
c --- 'ustrmn' = minimum ustar
c --- 'cormn4' = 4 times minimum coriolis magnitude
c --- 'csubp'  = specific heat of air at constant pressure (j/kg/deg)
c --- 'evaplh' = latent heat of evaporation (j/kg)
C --- 'csice'  = ice-air sensible exchange coefficient
c
      real       ustrmn,cormn4,csubp,evaplh,csice
      parameter (ustrmn=1.0e-5, 
     &           cormn4=4.0e-5,  ! corio(4N) is about 1.e-5
     &           csubp =1005.7,
     &           evaplh=2.47e6,
     &           csice =0.0006)
c
c --- parameters primarily for flxflg=1 (ustflg=1)
c --- 'airdns' = air density at sea level (kg/m**3)
c --- 'cd'     = drag coefficient
c --- 'ctl'    = thermal transfer coefficient (latent)
c --- 'cts1'   = thermal transfer coefficient (sensible, stable)
c --- 'cts2'   = thermal transfer coefficient (sensible, unstable)
c
      real       airdns,cd,ctl,cts1,cts2
      parameter (airdns=1.2)
      parameter (cd  =0.0013, ctl =0.0012,
     &           cts1=0.0012, cts2=0.0012)
c
c --- parameters primarily for flxflg=2 (ustflg=2)
c --- 'pairc'  = air pressure (mb) * 100
c --- 'rgas'   = gas constant (j/kg/k)
c --- 'tzero'  = celsius to kelvin temperature offset
c --- 'clmin'  = minimum allowed cl
c --- 'clmax'  = maximum allowed cl
c --- 'wsmin'  = minimum allowed wind speed (for cl and cd)
c --- 'wsmax'  = maximum allowed wind speed (for cl and cd)
c
      real       pairc,rgas,tzero,clmin,clmax,wsmin,wsmax
      parameter (pairc=1013.0*100.0,
     &           rgas =287.1,   tzero=273.16, 
     &           clmin=0.0003,  clmax=0.002,
     &           wsmin=3.5,     wsmax=27.5)
c
c --- parameters primarily for flxflg=4
c --- 'lvtc'   = include a virtual temperature correction
c --- 'vamin'  = minimum allowed wind speed (for cl)
c --- 'vamax'  = maximum allowed wind speed (for cl)
c --- 'tdmin'  = minimum allowed Ta-Ts      (for cl)
c --- 'tdmax'  = maximum allowed Ta-Ts      (for cl)
c
c --- 'as0_??' =  stable Ta-Ts  polynominal coefficients, va<=5m/s
c --- 'as5_??' =  stable Ta-Ts  polynominal coefficients, va>=5m/s
c --- 'au0_??' = unstable Ta-Ts polynominal coefficients, va<=5m/s
c --- 'au5_??' = unstable Ta-Ts polynominal coefficients, va>=5m/s
c --- 'an0_??' =  neutral Ta-Ts polynominal coefficients, va<=5m/s
c --- 'an5_??' =  neutral Ta-Ts polynominal coefficients, va>=5m/s
c --- 'ap0_??' =    +0.75 Ta-Ts polynominal coefficients, va<=5m/s
c --- 'ap5_??' =    +0.75 Ta-Ts polynominal coefficients, va>=5m/s
c --- 'am0_??' =    -0.75 Ta-Ts polynominal coefficients, va<=5m/s
c --- 'am5_??' =    -0.75 Ta-Ts polynominal coefficients, va>=5m/s
c
      logical, parameter :: lvtc =.true.
      real,    parameter :: vamin= 1.2, vamax=40.0
      real,    parameter :: tdmin=-8.0, tdmax= 7.0
c
      real, parameter ::
     &  as0_00=-2.925e-4,   as0_10= 7.272e-5,  as0_20=-6.948e-6,
     &  as0_01= 5.498e-4,   as0_11=-1.740e-4,  as0_21= 1.637e-5,
     &  as0_02=-5.544e-5,   as0_12= 2.489e-5,  as0_22=-2.618e-6
      real, parameter ::
     &  as5_00= 1.023e-3,   as5_10=-2.672e-6,  as5_20= 1.546e-6,
     &  as5_01= 9.657e-6,   as5_11= 2.103e-4,  as5_21=-6.228e-5,
     &  as5_02=-2.281e-8,   as5_12=-5.329e-3,  as5_22= 5.094e-4
      real, parameter ::
     &  au0_00= 2.077e-3,   au0_10=-2.899e-4,  au0_20=-1.954e-5,
     &  au0_01=-3.933e-4,   au0_11= 7.350e-5,  au0_21= 5.483e-6,
     &  au0_02= 3.971e-5,   au0_12=-6.267e-6,  au0_22=-4.867e-7
      real, parameter ::
     &  au5_00= 1.074e-3,   au5_10= 6.912e-6,  au5_20= 1.849e-7,
     &  au5_01= 5.579e-6,   au5_11=-2.244e-4,  au5_21=-2.167e-6,
     &  au5_02= 5.263e-8,   au5_12=-1.027e-3,  au5_22=-1.010e-4
      real, parameter ::
     &  an0_00= 1.14086e-3, an5_00= 1.073e-3,
     &  an0_01=-3.120e-6,   an5_01= 5.531e-6,
     &  an0_02=-9.300e-7,   an5_02= 5.433e-8
      real, parameter ::
     &  ap0_00= as0_00 + as0_10*0.75 + as0_20*0.75**2,
     &  ap0_01= as0_01 + as0_11*0.75 + as0_21*0.75**2,
     &  ap0_02= as0_02 + as0_12*0.75 + as0_22*0.75**2
      real, parameter ::
     &  ap5_00= as5_00 + as5_10*0.75 + as5_20*0.75**2,
     &  ap5_01= as5_01,
     &  ap5_02= as5_02,
     &  ap5_11=          as5_11*0.75 + as5_21*0.75**2,
     &  ap5_12=          as5_12*0.75 + as5_22*0.75**2
      real, parameter ::
     &  am0_00= au0_00 - au0_10*0.75 + au0_20*0.75**2,
     &  am0_01= au0_01 - au0_11*0.75 + au0_21*0.75**2,
     &  am0_02= au0_02 - au0_12*0.75 + au0_22*0.75**2
      real, parameter ::
     &  am5_00= au5_00 - au5_10*0.75 + au5_20*0.75**2,
     &  am5_01= au5_01,
     &  am5_02= au5_02,
     &  am5_11=        - au5_11*0.75 + au5_21*0.75**2,
     &  am5_12=        - au5_12*0.75 + au5_22*0.75**2
c
      real qsatur
      include 'stmt_fns.h'
c
c --- saturation mixing ratio (kg/kg), from a polynominal approximation
c --- for saturation vapor pressure (lowe, j.appl.met., 16, 100-103, 1976)
      qsatur(t)=.622e-3*(6.107799961e+00+t*(4.436518521e-01
     &               +t*(1.428945805e-02+t*(2.650648471e-04
     &               +t*(3.031240396e-06+t*(2.034080948e-08
     &               +t* 6.136820929e-11))))))
c
c --- salinity relaxation coefficient
      rmus=1./(30.0*86400.0)  !1/30 days
c
c --- temperature relaxation coefficient
      rmut=1./(30.0*86400.0)  !1/30 days
c
c --- ------------------------------------------------------
c --- thermal forcing of ocean surface (positive into ocean)
c --- ------------------------------------------------------
c
!DIR$ NO STREAM
      do l=1,isp(j)
c
!DIR$ CONCURRENT
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      if     (flxflg.gt.0) then
c ---   wind = wind speed (m/s)
        wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
     &      +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3
c ---   radfl= net       radiative thermal flux (w/m**2) positive into ocean
        radfl=radflx(i,j,l0)*w0+radflx(i,j,l1)*w1
     &       +radflx(i,j,l2)*w2+radflx(i,j,l3)*w3
c ---   swfl = shortwave radiative thermal flux (w/m**2) positive into ocean
        swfl =swflx (i,j,l0)*w0+swflx (i,j,l1)*w1
     &       +swflx (i,j,l2)*w2+swflx (i,j,l3)*w3
        if     (lwflag.gt.0) then
c ---     longwave correction to radfl (Qsw+Qlw).
          if     (lwflag.eq.1) then !from climatology (ignore ice)
            tsur = temp(i,j,1,n)
            tdif = tsur -
     &             ( twall(i,j,1,lc0)*wc0+twall(i,j,1,lc1)*wc1
     &              +twall(i,j,1,lc2)*wc2+twall(i,j,1,lc3)*wc3)
          else !w.r.t. atmospheric model's sst (allow for ice)
            tsur = (1.-covice(i,j))*temp(  i,j,1,n) +
     &                 covice(i,j) *temice(i,j)      
            tdif = tsur -
     &             ( surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1
     &              +surtmp(i,j,l2)*w2+surtmp(i,j,l3)*w3)
          endif
          !correction is blackbody radiation from tdif at tsur
          radfl = radfl - (4.506+0.0554*tsur) * tdif
          !count the correction as a relaxation term
          sstflx(i,j) = - (4.506+0.0554*tsur) * tdif
        else
          sstflx(i,j) = 0.0
        endif
        if     (pcipf) then
c ---     prcp = precipitation (m/sec; positive into ocean)
c ---     note that if empflg==3, this is actually P-E
          prcp=precip(i,j,l0)*w0+precip(i,j,l1)*w1
     &        +precip(i,j,l2)*w2+precip(i,j,l3)*w3
        endif
        if     (flxflg.ne.3) then
c ---     airt = air temperature (C)
          airt=airtmp(i,j,l0)*w0+airtmp(i,j,l1)*w1
     &        +airtmp(i,j,l2)*w2+airtmp(i,j,l3)*w3
c ---     vpmx = water vapor mixing ratio (kg/kg)
          vpmx=vapmix(i,j,l0)*w0+vapmix(i,j,l1)*w1
     &        +vapmix(i,j,l2)*w2+vapmix(i,j,l3)*w3
        endif
c ---   ustar = U* (sqrt(N.m/kg))                 
        if     (ustflg.eq.3) then !ustar from input
          ustar(i,j)=ustara(i,j,l0)*w0+ustara(i,j,l1)*w1
     &              +ustara(i,j,l2)*w2+ustara(i,j,l3)*w3
        elseif (ustflg.eq.1) then !ustar from wndspd, constant cd
          ustar(i,j)=sqrt(thref*cd*airdns)*wind
        elseif (ustflg.eq.2) then !ustar from wndspd, variable cd
          wsph = min( wsmax, max( wsmin, wind ) )
          cd0  = 0.862e-3 + 0.088e-3 * wsph - 0.00089e-3 * wsph**2
          rair = pairc / (rgas * ( tzero + airt ))
          ustar(i,j)=sqrt(thref*cd0*rair)*wind
        elseif (ustflg.eq.4) then !ustar from wind stress
          if     (wndflg.eq.2) then ! tau on p grid
            xtau=( taux(i,j,l0)*w0+taux(i,j,l1)*w1
     &            +taux(i,j,l2)*w2+taux(i,j,l3)*w3)
            ytau=( tauy(i,j,l0)*w0+tauy(i,j,l1)*w1
     &            +tauy(i,j,l2)*w2+tauy(i,j,l3)*w3)
          else ! tau on u&v grids
            xtau=( (taux(i,j,l0)+taux(i+1,j,l0))*w0
     &            +(taux(i,j,l1)+taux(i+1,j,l1))*w1
     &            +(taux(i,j,l2)+taux(i+1,j,l2))*w2
     &            +(taux(i,j,l3)+taux(i+1,j,l3))*w3)*0.5
            ytau=( (tauy(i,j,l0)+tauy(i,j+1,l0))*w0
     &            +(tauy(i,j,l1)+tauy(i,j+1,l1))*w1
     &            +(tauy(i,j,l2)+tauy(i,j+1,l2))*w2
     &            +(tauy(i,j,l3)+tauy(i,j+1,l3))*w3)*0.5
          endif
          ustar(i,j)=sqrt(thref*sqrt(xtau**2+ytau**2))
        endif !ustflg
        ustar( i,j)=max(ustrmn,ustar(i,j))
        hekman(i,j)=ustar(i,j)*(cekman*4.0)/
     &               max( cormn4,
     &                    abs(corio(i,j  ))+abs(corio(i+1,j  ))+
     &                    abs(corio(i,j+1))+abs(corio(i+1,j+1)) )
      else !flxlfg==0, i.e. no flux
        swfl=0.0
        ustar( i,j)=0.0
        hekman(i,j)=0.0
      endif !flxflg
c
      if     (flxflg.eq.1) then
c
c ---   MICOM bulk air-sea flux parameterization
c ---   (constant Cl and constant stable/unstable Cs)
c
        if (temp(i,j,1,n).lt.airt) then
          csh=cts1  !stable
        else
          csh=cts2  !unstable
        endif
c ---   evap   = evaporation (w/m**2) into atmos.
c ---   snsibl = sensible heat flux   into atmos.
        if (.not.icegln) then
          evap  =ctl*airdns*evaplh*wind*
     &           max(0.,0.97*qsatur(temp(i,j,1,n))-vpmx)
          snsibl=csh*airdns*csubp*wind*(temp(i,j,1,n)-airt)
        else  ! ice
          evap  =(1.-covice(i,j))*
     &           ctl*airdns*evaplh*wind*
     &           max(0.,0.97*qsatur(temp(i,j,1,n))-vpmx)
          snsibl=airdns*csubp*wind*
     &           ( (1.-covice(i,j))*csh  *(temp(i,j,1,n)-airt) +
     &                 covice(i,j) *csice*(temice(i,j)  -airt)  )
        endif
c ---   surflx = thermal energy flux (w/m**2) into ocean
        surflx(i,j)=radfl - snsibl - evap
      elseif (flxflg.eq.2) then
c
c ---    Cl (and Cs) depend on wind speed and Ta-Ts.
c ---    Kara, A. B., P. A. Rochford, and H. E. Hurlburt, 2002:
c ---    Air-sea flux estimates and the 1997-1998 ENSO event.
c ---    Bound.-Layer Meteor., 103, 439-458.
c
        rair = pairc / (rgas * ( tzero + airt ))
        slat = evaplh*rair
        ssen = csubp *rair
c
        tdif = temp(i,j,1,n) - airt
        wsph = min( wsmax, max( wsmin, wind ) )
        cl0  =  0.885e-3 + 0.0748e-3 * wsph - 0.00143e-3 * wsph**2
        cl1  = -0.113e-4 + 4.89e-4   / wsph
        clh  = min( clmax, max( clmin, cl0 + cl1 * tdif ) )
        csh  = 0.9554*clh
c
c ---   evap   = evaporation         (w/m**2) into atmos.
c ---   snsibl = sensible heat flux  (w/m**2) into atmos.
c ---   surflx = thermal energy flux (w/m**2) into ocean
        if (.not.icegln) then
          evap   = slat*clh*wind*(0.97*qsatur(temp(i,j,1,n))-vpmx)
          snsibl = ssen*csh*wind* tdif
        else  ! ice
          evap   = (1.-covice(i,j))*
     &             slat*clh*wind*(0.97*qsatur(temp(i,j,1,n))-vpmx)
          snsibl = ssen*wind*
     &             ( (1.-covice(i,j))*csh  * tdif +
     &                   covice(i,j) *csice*(temice(i,j)-airt) )
        endif
        surflx(i,j) = radfl - snsibl - evap
c
cdiag   if     (i.eq.itest.and.j.eq.jtest) then
cdiag     write(lp,'(i9,2i5,a,4f8.5)')
cdiag.    nstep,i0+i,j0+j,' cl0,cl,cs,cd    = ',cl0,clh,csh,cd0
cdiag     write(lp,'(i9,2i5,a,2f8.2,f8.5)')
cdiag.    nstep,i0+i,j0+j,' wsph,tdif,ustar = ',wsph,tdif,ustar(i,j)
cdiag     call flush(lp)
cdiag   endif
      elseif (flxflg.eq.4) then
c
c ---   Similar to flxflg.eq.2, but with Cl based on an approximation
c ---   to values from the COARE 3.0 algorithm (Fairall et al. 2003), 
c ---   for Cl over the global ocean in the range 1m/s <= Va <= 40m/s
c ---   and -8degC <= Ta-Ts <= 7degC, that is quadratic in Ta-Ts and
c ---   quadratic in either Va or 1/Va (Kara and Hurlburt, 2004).
c
c ---   Fairall, C. W., E. F. Bradley, J. E. Hare, A. A. Grachev, and J. B.
c ---   Edson, 2003:  Bulk parameterization of air-sea fluxes:  Updates 
c ---   and verification for the COARE algorithm.  J. Climate, 16, 571-591.
c
c ---   Kara, K. Birol, H.E. Hurlburt, 2004: A Note on the Stability-Dependent
c ---   Exchange Coefficients of Air-Sea Fluxes for Use in General Circulation
c ---   Models. submitted to J. Atmos. Oceanic Technol. on January 14, 2004. 
c ---   preprint at http://hycom.rsmas.miami.edu/publications.html
c
        rair = pairc / (rgas * ( tzero + airt ))
        slat = evaplh*rair
        ssen = csubp *rair
c
        tdif  = temp(i,j,1,n) - airt
        if     (lvtc) then !include a virtual temperature correction
          tamts = -tdif - 0.61*(airt+tzero)*(qsatur(airt)-vpmx)
          tamts = min( tdmax, max( tdmin, tamts ) )
        else
          tamts = min( tdmax, max( tdmin, -tdif ) )
        endif !lvtc:else
        va    = min( vamax, max( vamin,  wind ) )
        if     (va.le.5.0) then
          if     (tamts.gt. 0.75) then !stable
            clh =   (as0_00 + as0_01* va + as0_02* va**2)
     &            + (as0_10 + as0_11* va + as0_12* va**2)*tamts
     &            + (as0_20 + as0_21* va + as0_22* va**2)*tamts**2
          elseif (tamts.lt.-0.75) then !unstable
            clh =   (au0_00 + au0_01* va + au0_02* va**2)
     &            + (au0_10 + au0_11* va + au0_12* va**2)*tamts
     &            + (au0_20 + au0_21* va + au0_22* va**2)*tamts**2
          elseif (tamts.ge.-0.098)  then
            q = (tamts-0.75)/0.848  !linear between  0.75 and -0.098
            q = q**2  !favor  0.75
            clh = (1.0-q)*(ap0_00 + ap0_01* va + ap0_02* va**2)
     &               + q *(an0_00 + an0_01* va + an0_02* va**2)
          else
            q = (tamts+0.75)/0.652  !linear between -0.75 and -0.098
            q = q**2  !favor -0.75
            clh = (1.0-q)*(am0_00 + am0_01* va + am0_02* va**2)
     &               + q *(an0_00 + an0_01* va + an0_02* va**2)
          endif !tamts
        else !va>5
          qva = 1.0/va
          if     (tamts.gt. 0.75) then !stable
            clh =   (as5_00 + as5_01* va + as5_02* va**2)
     &            + (as5_10 + as5_11*qva + as5_12*qva**2)*tamts
     &            + (as5_20 + as5_21*qva + as5_22*qva**2)*tamts**2
          elseif (tamts.lt.-0.75) then !unstable
            clh =   (au5_00 + au5_01* va + au5_02* va**2)
     &            + (au5_10 + au5_11*qva + au5_12*qva**2)*tamts
     &            + (au5_20 + au5_21*qva + au5_22*qva**2)*tamts**2
          elseif (tamts.ge.-0.098)  then
            q = (tamts-0.75)/0.848  !linear between  0.75 and -0.098
            q = q**2  !favor  0.75
            clh = (1.0-q)*(ap5_00 + ap5_01* va + ap5_02* va**2
     &                            + ap5_11*qva + ap5_12*qva**2)
     &               + q *(an5_00 + an5_01* va + an5_02* va**2)
          else
            q = (tamts+0.75)/0.652  !linear between -0.75 and -0.098
            q = q**2  !favor -0.75
            clh = (1.0-q)*(am5_00 + am5_01* va + am5_02* va**2
     &                            + am5_11*qva + am5_12*qva**2)
     &               + q *(an5_00 + an5_01* va + an5_02* va**2)
          endif !tamts
        endif !va
        csh  = 0.9554*clh
c
c ---   evap   = evaporation         (w/m**2) into atmos.
c ---   snsibl = sensible heat flux  (w/m**2) into atmos.
c ---   surflx = thermal energy flux (w/m**2) into ocean
        if (.not.icegln) then
          evap   = slat*clh*wind*(0.97*qsatur(temp(i,j,1,n))-vpmx)
          snsibl = ssen*csh*wind* tdif
        else  ! ice
          evap   = (1.-covice(i,j))*
     &             slat*clh*wind*(0.97*qsatur(temp(i,j,1,n))-vpmx)
          snsibl = ssen*wind*
     &             ( (1.-covice(i,j))*csh  * tdif +
     &                   covice(i,j) *csice*(temice(i,j)-airt) )
        endif
        surflx(i,j) = radfl - snsibl - evap
c
cdiag   if     (i.eq.itest.and.j.eq.jtest) then
cdiag     write(lp,'(i9,2i5,a,3f8.5)')
cdiag.    nstep,i0+i,j0+j,' cl,cs,cd    = ',clh,csh,cd0
cdiag     write(lp,'(i9,2i5,a,2f8.2,f8.5)')
cdiag.    nstep,i0+i,j0+j,' va,tamst,ustar = ',va,tamts,ustar(i,j)
cdiag     call flush(lp)
cdiag   endif
      elseif (flxflg.eq.3) then
c
c ---   input radiation flux is the net flux.
c
        evap=0.0
        surflx(i,j)=radfl
      else  ! no flux
        evap=0.0
        surflx(i,j)=0.0  
      endif  ! flxflg
c
c --- add a time-invarient net heat flux offset
      if     (flxoff) then
        surflx(i,j)=surflx(i,j)+offlux(i,j)
      endif
c
c --- relax to surface temperature
c --- use a reference relaxation thickness (min. mixed layer depth)
c --- in shallow water, thkmlt is replaced by the total depth
c --- actual e-folding time is (dpmixl(i,j,n)/(thkmlt*onem))/rmut
c --- in shallow water this is (dpmixl(i,j,n)/p(i,j,kk+1)  )/rmut
      if     (sstflg.eq.1) then !climatological sst
        sstrlx=
     &   (rmut*spcifh*min(p(i,j,kk+1),thkmlt*onem)/g)*
     &   ( ( twall(i,j,1,lc0)*wc0+twall(i,j,1,lc1)*wc1
     &      +twall(i,j,1,lc2)*wc2+twall(i,j,1,lc3)*wc3) -
     &     temp(i,j,1,n) )
        surflx(i,j)=surflx(i,j)+sstrlx
        sstflx(i,j)=sstflx(i,j)+sstrlx
      elseif (sstflg.gt.1) then !synoptic sst
        sstrlx=
     &   (rmut*spcifh*min(p(i,j,kk+1),thkmlt*onem)/g)*
     &   ( ( seatmp(i,j,l0)*w0+seatmp(i,j,l1)*w1
     &      +seatmp(i,j,l2)*w2+seatmp(i,j,l3)*w3) -
     &     temp(i,j,1,n) )
        surflx(i,j)=surflx(i,j)+sstrlx
        sstflx(i,j)=sstflx(i,j)+sstrlx
      endif
c --- sswflx = shortwave radiative energy flux (w/m**2) into ocean
      if (.not.icegln) then
        sswflx(i,j)=swfl
      else
        sswflx(i,j)=swfl*(1.-covice(i,j))
      endif
c --- emnp = evaporation minus precipitation (m/sec) into atmos.
      if     (.not.pcipf) then
        prcp = 0.0
        emnp = 0.0   !no E-P
      elseif (empflg.eq.3) then
        emnp = -prcp !input prcp is P-E
      else
        emnp = evap*thref/evaplh - prcp  !input prcp is P
      endif
c --- allow for rivers as a precipitation bogas
      if     (priver) then
        emnp=emnp - ( rivers(i,j,lr0)*wr0+rivers(i,j,lr1)*wr1
     &               +rivers(i,j,lr2)*wr2+rivers(i,j,lr3)*wr3)
      endif
c --- salflx = salt flux (10**-3 kg/m**2/sec) into ocean
      salflx(i,j)=emnp*(saln(i,j,1,n)*qthref)
c --- relax to surface salinity
      if     (srelax) then
c ---   use a reference relaxation thickness (min. mixed layer depth)
c ---   in shallow water, thkmls is replaced by the total depth
c ---   actual e-folding time is (dpmixl(i,j,n)/(thkmls*onem))/rmus
c ---   in shallow water this is (dpmixl(i,j,n)/p(i,j,kk+1)  )/rmus
        sssflx(i,j)=
     &     (rmus*min(p(i,j,kk+1),thkmls*onem)/g)*
     &     ( ( swall(i,j,1,lc0)*wc0+swall(i,j,1,lc1)*wc1
     &        +swall(i,j,1,lc2)*wc2+swall(i,j,1,lc3)*wc3) -
     &       saln(i,j,1,n) )
        salflx(i,j)=salflx(i,j)+sssflx(i,j)
      else
        sssflx(i,j)=0.0
      endif !srelax
c
      if     (epmass) then
c ---   change total water depth by the water exchanged with the atmosphere
        pbavg(i,j,1)=pbavg(i,j,1)-
     &               onem*0.5*delt1*salflx(i,j)/(saln(i,j,1,n)*qthref)
        pbavg(i,j,2)=pbavg(i,j,2)-
     &               onem*0.5*delt1*salflx(i,j)/(saln(i,j,1,n)*qthref)
        pbavg(i,j,3)=pbavg(i,j,3)-
     &               onem*0.5*delt1*salflx(i,j)/(saln(i,j,1,n)*qthref)
      endif !epmass
c
c --- n o t e : t/s changes due to surflx/salflx are computed in mxlayr routine
c
      util1(i,j)=  surflx(i,j)*scp2(i,j)
      util2(i,j)=  salflx(i,j)*scp2(i,j)
      util3(i,j)=temp(i,j,1,n)*scp2(i,j)
      util4(i,j)=saln(i,j,1,n)*scp2(i,j)
c
cdiag if     (i.eq.itest.and.j.eq.jtest) then
cdiag   write (lp,100) 
cdiag.  nstep,i0+i,j0+j,
cdiag.  '   radfl       swfl        wind        airt         sst',
cdiag.  radfl,swfl,wind,airt,temp(i,j,1,n),
cdiag.  nstep,i0+i,j0+j,
cdiag.  '   vpmx       snsibl       evap        prcp        emnp',
cdiag.  vpmx,snsibl,evap,prcp,emnp,
cdiag.  nstep,i0+i,j0+j,
cdiag.  '  surflx      sswflx      salflx      hekman       ustar',
cdiag.  surflx(i,j),sswflx(i,j),salflx(i,j),hekman(i,j),ustar(i,j)
cdiag   call flush(lp)
cdiag endif
 100  format(i9,2i5,a/18x,1p5e12.4/
     &       i9,2i5,a/18x,1p5e12.4/
     &       i9,2i5,a/18x,1p5e12.4)
      enddo !i
      enddo !l
      return
      end subroutine thermfj
c
c
c> Revision history:
c>
c> Oct. 1999 - surface flux calculations modified for kpp mixed layer model,
c>             including penetrating solar radiation based on jerlov water type
c> Apr. 2000 - conversion to SI units
c> Oct  2000 - added thermfj to simplify OpenMP logic
c> Dec  2000 - modified fluxes when ice is present
c> Dec  2000 - added Kara bulk air-sea flux parameterization (flxflg=2)
c> May  2002 - buoyfl now calculated in mixed layer routine
c> Aug  2002 - added nested velocity relaxation
c> Nov  2002 - separate sss and sst relaxation time scales (thkml[st])
c> Nov  2002 - save sssflx and sstflx for diagnostics
c> Mar  2003 - longwave radiation correction for model vs "longwave" SST
c> May  2003 - use seatmp in place of twall.1, when available
c> Mar  2003 - add option to smooth surface fluxes
c> Mar  2004 - added epmass for treating E-P as a mass exchange
c> Mar  2005 - limit thkml[st] to no more than the actual depth
c> Mar  2005 - added empflg
c> Mar  2005 - replaced qsatur with 97% of qsatur in evap calculation
c> Mar  2005 - added ustflg
c> Mar  2005 - added flxoff
c> Apr  2005 - add a virtual temperature correction to Ta-Ts for flxflg=4.
