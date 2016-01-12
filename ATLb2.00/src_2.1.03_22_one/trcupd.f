      subroutine initrc(mnth)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
      implicit none
c
      include 'common_blocks.h'
c
      integer mnth
c
c --- --------------------------
c --- initializatize all tracers
c --- --------------------------
c
      logical    lpipe_initrc
      parameter (lpipe_initrc=.false.)
c
      character ptxt*12,cformat*99
      integer   i,ibio,nbio,j,k,ktr,l
      real      bio_n,bio_p,zk
c
c --- expand trcflg to allow for number of biology fields.
c
      nbio = 0
      ibio = 0
      do ktr= 1,ntracr+1
        if     (ktr.ne.ntracr+1 .and.
     &          trcflg(min(ktr,ntracr)).eq.9) then
          if     (ibio.eq.0) then !start biology
            ibio = ktr
          endif
        elseif (ibio.ne.0) then !end biology
          nbio = ktr-ibio
          if     (nbio.eq.3) then
c ---       Franks NPZ.
            trcflg(ibio)   =  903
            trcflg(ibio+1) = -903
            trcflg(ibio+2) = -903
            ibio = 0
          elseif (nbio.eq.3) then
c ---       Two Franks NPZ.
            trcflg(ibio)   =  903
            trcflg(ibio+1) = -903
            trcflg(ibio+2) = -903
            trcflg(ibio+3) =  903
            trcflg(ibio+4) = -903
            trcflg(ibio+5) = -903
            ibio = 0
          elseif (nbio.eq.4) then
c ---       Lima/Idrisi NPZD.
            trcflg(ibio)   =  904
            trcflg(ibio+1) = -904
            trcflg(ibio+2) = -904
            trcflg(ibio+3) = -904
            ibio = 0
          elseif (nbio.eq.7) then
c ---       Lima/Idrisi NPZD and Franks NPZ.
            trcflg(ibio)   =  904
            trcflg(ibio+1) = -904
            trcflg(ibio+2) = -904
            trcflg(ibio+3) = -904
            trcflg(ibio+4) =  903
            trcflg(ibio+5) = -903
            trcflg(ibio+6) = -903
            ibio = 0
          elseif (nbio.eq.8) then
c ---       Two Lima/Idrisi NPZD.
            trcflg(ibio)   =  904
            trcflg(ibio+1) = -904
            trcflg(ibio+2) = -904
            trcflg(ibio+3) = -904
            trcflg(ibio+4) =  904
            trcflg(ibio+5) = -904
            trcflg(ibio+6) = -904
            trcflg(ibio+7) = -904
            ibio = 0
          elseif (nbio.eq.9) then
c ---       Chai 9-component.
*           trcflg(ibio)   =  909
*           trcflg(ibio+1) = -909
*           trcflg(ibio+2) = -909
*           trcflg(ibio+3) = -909
*           trcflg(ibio+4) = -909
*           trcflg(ibio+5) = -909
*           trcflg(ibio+6) = -909
*           trcflg(ibio+7) = -909
*           trcflg(ibio+8) = -909
*           ibio = 0
c ---       not yet implemented
            if (mnproc.eq.1) then
            write(lp,'(/ 3a /)')
     &        'error - trcflg=9 (standard biology) configured',
     &        ' with 9 consecutive tracers, but Chai scheme is',
     &        ' not yet implemented'
            call flush(lp)
            endif !1st tile
            call xcstop('(trcini)')
                   stop '(trcini)'
          else
c ---       unknown standard biology.
            if (mnproc.eq.1) then
            write(lp,'(/ 2a,i3 /)')
     &        'error - trcflg=9 (standard biology) expects',
     &        ' 3/4/6/7/8 consecutive tracers but have',nbio
*    &        ' 3/4/6/7/8/9 consecutive tracers but have',nbio
            call flush(lp)
            endif !1st tile
            call xcstop('(trcini)')
                   stop '(trcini)'
          endif
        endif
      enddo
c
      if (ntracr.gt.0) then
        if (mnproc.eq.1) then
        write(lp,*)
        do k= 1,ntracr
          write(lp,'(a,i3,i6)') 'initrc: k,trcflg =',k,trcflg(k)
        enddo
        write(lp,*)
        endif !1st tile
      endif !ntracr.gt.0
c
      if     (nbio.gt.0) then
c
c ---   input bio-tracer parameters.
c ---   note that multiple sets of bio-tracers are allowed,
c ---   each is read from tracer.input in tracer order.
c
        open(unit=99,file='tracer.input')
        do ktr= 1,ntracr
          if     (trcflg(ktr).eq.903) then
c ---       NPZ
            call trcupd_903(1,2, -ktr)
          elseif (trcflg(ktr).eq.904) then
c ---       NPZD
            call trcupd_904(1,2, -ktr)
*         elseif (trcflg(ktr).eq.909) then
* ---       Chai 9-component.
*           call trcupd_909(1,2, -ktr)
          endif
        enddo
        close(unit=99)
      endif
c
      if     (.not.trcout .or. trcrin) then
        return  ! no tracer or tracer from restart
      endif
c
      margin = 0
c
      if     (iniflg.eq.2) then  ! use climatology
        call rdrlax(mnth,1)
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,ktr)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              p(i,j,kk+1)=depths(i,j)*onem
              p(i,j,   1)=0.0
              do k=1,kk
                p(i,j,k+1)=p(i,j,k)+dp(i,j,k,1)
                if     (p(i,j,k).lt.p(i,j,kk+1)-tencm) then
                  do ktr= 1,ntracr
                    tracer(i,j,k,1,ktr)=trwall(i,j,k,1,ktr)
                    tracer(i,j,k,2,ktr)=trwall(i,j,k,1,ktr)
                  enddo !ktr
                else
                  do ktr= 1,ntracr
                    tracer(i,j,k,1,ktr)=tracer(i,j,k-1,1,ktr)
                    tracer(i,j,k,2,ktr)=tracer(i,j,k-1,1,ktr)
                  enddo !ktr
                endif
              enddo !k
            enddo !i
          enddo !l
        enddo !j
      else ! analytic inititalization
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,ktr)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              p(i,j,1)=0.0
              do k=1,kk
                p(i,j,k+1)=p(i,j,k)+dp(i,j,k,1)
                do ktr= 1,ntracr
                  if     (trcflg(ktr).eq.0) then !100% in the mixed layer
                    if     (p(i,j,k).le.dpmixl(i,j,1)) then
                      tracer(i,j,k,1,ktr)=10.0
                      tracer(i,j,k,2,ktr)=10.0
                    else
                      tracer(i,j,k,1,ktr)=0.0
                      tracer(i,j,k,2,ktr)=0.0
                    endif
                  elseif (trcflg(ktr).eq.1) then !20 below euphotic zone
                    if     (p(i,j,k)*betabl(jerlv0).lt.4.0) then
                      tracer(i,j,k,1,ktr)=0.0
                      tracer(i,j,k,2,ktr)=0.0
                    else
                      tracer(i,j,k,1,ktr)=20.0  ! mg/m^3
                      tracer(i,j,k,2,ktr)=20.0  ! mg/m^3
                    endif
                  elseif (trcflg(ktr).eq.2) then !temperature
                    tracer(i,j,k,1,ktr)=temp(i,j,k,1)
                    tracer(i,j,k,2,ktr)=temp(i,j,k,1)
                  elseif (trcflg(ktr).eq.904 .or.
     &                    trcflg(ktr).eq.903     ) then !NPZD or NPZ
                    zk = 0.5*(p(i,j,k+1)+p(i,j,k))*qonem
                    if     (zk.le.300.0) then
                      ! 0.1 at 300m, 1.0 at 100m, 2.025 at 0m
                      bio_p = 0.1 + (300.0-zk)**2 * (0.9/200.0**2)
                    elseif (zk.le.900.0) then
                      ! 0.1 at 300m, 0.0 at 900m
                      bio_p = (900.0-zk) * 0.1/600.0
                    else
                      bio_p = 0.0
                    endif
                    if     (temp(i,j,k,1).lt. 6.0) then
                      bio_n = 37.0
                    elseif (temp(i,j,k,1).gt.27.0) then
                      bio_n =  0.0
                    else
*                     bio_n = (27.0-temp(i,j,k,1)) * 37.0/21.0
                      bio_n = 39.3116-1.335*temp(i,j,k,1)
                    endif
                    tracer(i,j,k,1,ktr  )=bio_n  !N
                    tracer(i,j,k,2,ktr  )=bio_n
                    tracer(i,j,k,1,ktr+1)=bio_p  !P
                    tracer(i,j,k,2,ktr+1)=bio_p
                    tracer(i,j,k,1,ktr+2)=bio_p  !Z=P
                    tracer(i,j,k,2,ktr+2)=bio_p
                    if     (trcflg(ktr).eq.904) then
                      tracer(i,j,k,1,ktr+3)=bio_p + 1.0  !D=P+1
                      tracer(i,j,k,2,ktr+3)=bio_p + 1.0
                    endif
                  endif !trcflg
                enddo !ktr
              enddo !k
            enddo !i
          enddo !l
        enddo !j
      endif !iniflg.eq.2:else
c
      if     (lpipe .and. lpipe_initrc) then
         do ktr= 1,ntracr
           do k= 1,kk
             write (ptxt,'(a4,i2.2,a3,i3)') 'trc.',ktr,' k=',k
             call pipe_compare_sym1(tracer(1-nbdy,1-nbdy,k,1,ktr),
     &                              ip,ptxt)
           enddo !k
         enddo !ktr
       endif !lpipe.and.lpipe_initrc
c
      if     (itest.gt.0 .and. jtest.gt.0) then
         write(cformat,'(a,i2,a,i2,a)')
     &     '(i9,2i5,a,',ntracr,
     &     'a / (23x,i3,2f8.2,', ntracr,'f8.4))'
         write (lp,cformat)
     &     nstep,i0+itest,j0+jtest,
     &     '  istate:  thkns    dpth',
     &     ('  tracer',ktr=1,ntracr),
     &     (k,
     &      dp(itest,jtest,k,1)*qonem,
     &      (p(itest,jtest,k+1)+p(itest,jtest,k))*0.5*qonem,
     &      (tracer(itest,jtest,k,1,ktr),ktr=1,ntracr),
     &      k=1,kk)
         write(lp,'(23x,a,8x,f8.2)') 'bot',depths(itest,jtest)
      endif !test tile
      call xcsync(flush_lp)
c
      return
      end

      subroutine trcupd(m,n)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- -----------------------------------------------------------
c --- tracer-specific operations (side-wall relaxation in thermf)
c --- -----------------------------------------------------------
c
      integer i,j,jrlv,k,ktr,l
      real    pijk,pijkp,q
c
      margin = 0  ! no horizontal derivatives
c
      do ktr= 1,ntracr
        if     (trcflg(ktr).eq.0) then
c ---     tracer always 10.0 at surface
!$OMP     PARALLEL DO PRIVATE(j,k,l,i,ktr)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do l=1,isp(j)
              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                tracer(i,j,1,n,ktr) = 10.0
              enddo !i
            enddo !l
          enddo !j
        elseif (trcflg(ktr).eq.1) then
c ---     psudo-silicate, half-life of 30 days in euphotic zone
          q = 1.0-delt1/(30.0*86400.0)
!$OMP     PARALLEL DO PRIVATE(j,k,l,i,ktr,jrlv,pijk,pijkp)
!$OMP&             SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do l=1,isp(j)
              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                jrlv=jerlov(i,j)
                pijkp=0.0
                do k=1,kk
                  pijk  = pijkp
                  pijkp = pijk+dp(i,j,k,n)
                  if     (0.5*(pijk+pijkp)*betabl(jrlv).lt.4.0) then
                    tracer(i,j,k,n,ktr) = q*tracer(i,j,k,n,ktr)
                  else
                    exit  !too deep
                  endif
                enddo
              enddo !i
            enddo !l
          enddo !j
        elseif (trcflg(ktr).eq.2) then
c ---     temperature-like (do nothing, heat flux forcing in mixed layer)
        elseif (trcflg(ktr).eq.903) then
c ---     NPZ
          call trcupd_903(m,n, ktr)
        elseif (trcflg(ktr).eq.904) then
c ---     NPZD
          call trcupd_904(m,n, ktr)
*       elseif (trcflg(ktr).eq.909) then
* ---     Chai 9-component.
*         call trcupd_909(m,n, ktr)
        endif
      enddo !ktr
      return
      end subroutine trcupd

      subroutine trcupd_903(m,n, ibio)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n,ibio
c
c --- -------------------------------------------------
c --- tracer-specific operations for Franks NPZ biology
c --- -------------------------------------------------
c
      real,    save, dimension(mxtrcr) ::
     & bup,   ! maximum growth  rate of phytoplankton (1/d).
     & bgz,   ! maximum grazing rate of zooplankton   (1/d).
     & bdp,   ! senescence (death) rate of phytoplankton (1/d).
     & bdz,   ! death rate of zooplankton (1/d).
     & buk,   ! = half-saturation coefficient for phytoplankton (mg/m^3)
     & asim,  ! assimilation efficiency of zooplankton.
     & glam   ! Ivlev parameter for grazing efficiency of zooplankton.
c
      integer i,j,jrlv,k,l
      real    bm_n,bm_p,bm_z,bn_n,bn_p,bn_z,bu_n,bu_p,bu_z,
     &        uptake,grazin,pdeath,zdeath,
     &        pijk,pijkp,par
c
      if     (ibio.lt.0) then !initialize only
c
c ---   read from tracer_NN.input:
c ---   'biotyp' = type (90X=std.bio,X=3,4,9) must be 903
c ---   'bup   ' = maximum growth  rate of phytoplankton (1/d).
c ---   'bgz   ' = maximum grazing rate of zooplankton   (1/d).
c ---   'bdp   ' = senescence (death) rate of phytoplankton (1/d).
c ---   'bdz   ' = death rate of zooplankton (1/d).
c ---   'buk   ' = half-saturation coefficient for phytoplankton (mg/m^3)
c ---   'asim  ' = assimilation efficiency of zooplankton.
c ---   'glam  ' = Ivlev parameter for grazing efficiency of zooplankton.
c
        i = -ibio
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3,a,i3,a)')
     &    'Franks NPZ parameters for tracers',i,' to',i+2,':'
        endif !1st tile
c
        call blkini(k, 'biotyp')
        if     (k.ne.903) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')
     &        'error - biotyp must be 903'
          call flush(lp)
          endif !1st tile
          call xcstop('(trcini)')
                 stop '(trcini)'
        endif !biotyp.ne.903
c
        call blkinr(bup(   i), 'bup   ','(a6," =",f10.4," 1/d")')
        call blkinr(bgz(   i), 'bgz   ','(a6," =",f10.4," 1/d")')
        call blkinr(bdp(   i), 'bdp   ','(a6," =",f10.4," 1/d")')
        call blkinr(bdz(   i), 'bdz   ','(a6," =",f10.4," 1/d")')
        call blkinr(buk(   i), 'buk   ','(a6," =",f10.4," mg/m^3")')
        call blkinr(asim(  i), 'asim  ','(a6," =",f10.4," ")')
        call blkinr(glam(  i), 'glam  ','(a6," =",f10.4," ")')
c
        if (mnproc.eq.1) then
        write(lp,*)
        endif !1st tile
        return
      endif !ibio.lt.0
c
c --- leapfrog time step.
c
      margin = 0  ! no horizontal derivatives
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,jrlv,pijk,pijkp,par,
!$OMP&                    bm_n,bm_p,bm_z,bn_n,bn_p,bn_z,
!$OMP&                    bu_n,bu_p,bu_z,
!$OMP&                    uptake,grazin,pdeath,zdeath)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            jrlv=jerlov(i,j)
            pijkp=0.0
            do k=1,kk
              pijk  = pijkp
              pijkp = pijk+dp(i,j,k,n)
              par   = exp(-0.5*(pijk+pijkp)*betabl(jrlv))
c
              bm_n = tracer(i,j,k,m,ibio)
              bm_p = tracer(i,j,k,m,ibio+1)
              bm_z = tracer(i,j,k,m,ibio+2)
              bn_n = tracer(i,j,k,n,ibio)
              bn_p = tracer(i,j,k,n,ibio+1)
              bn_z = tracer(i,j,k,n,ibio+2)
c
              uptake = bup(ibio)*bm_p*bm_n*par/(buk(ibio)+bm_n)
              grazin = bgz(ibio)*bm_z*(1.0-exp(-glam(ibio)*bm_p))
              pdeath = bdp(ibio)*bm_p
              zdeath = bdz(ibio)*bm_z
              ! limit negative terms to 10% of total per single time step
              grazin = min(grazin,bn_p*0.2*86400.0/delt1)
              uptake = min(uptake,bn_n*0.2*86400.0/delt1)
c
              bu_p =                 -grazin       +uptake-pdeath
              bu_z =      asim(ibio) *grazin-zdeath
              bu_n = (1.0-asim(ibio))*grazin+zdeath-uptake+pdeath
c
              tracer(i,j,k,n,ibio)   = bn_n + delt1/86400.0 * bu_n
              tracer(i,j,k,n,ibio+1) = bn_p + delt1/86400.0 * bu_p
              tracer(i,j,k,n,ibio+2) = bn_z + delt1/86400.0 * bu_z
            enddo
          enddo !i
        enddo !l
      enddo !j
      return
      end subroutine trcupd_903

      subroutine trcupd_904(m,n, ibio)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n,ibio
c
c --- -------------------------------------------------------
c --- tracer-specific operations for Lima/Idrisi NPZD biology
c --- -------------------------------------------------------
c
      real,    save, dimension(mxtrcr) ::
     &  pp,   ! zoopl: preference term for phytoplankton
     &  pz,   ! zoopl: preference term for zooplankton
     &  pd,   ! zoopl: preference term for detritus
     &  aa,   ! zoopl: assimilation efficiency
     &  am,   ! zoopl: metabolic    efficiency
     &  fkz,  ! zoopl: half-saturation coefficient (mg/m^3)
     &  gmax, ! zoopl: maximum growth rate (1/day)
     &  zmor  ! zoopl: mortality (1/day)
c
      real,    save, dimension(mxtrcr) ::
*    &  ik,   ! phyto: light absorption efficiency scalar (einst/m^2/h)
     &  fkp,  ! phyto: half-saturation coefficient (mg/m^3)
     &  pmax, ! phyto: maximum growth rate (1/day)
     &  psen  ! phyto: senescence (1/day)
c
      real,    save, dimension(mxtrcr) ::
     &  remn  ! detri: remineralization (1/day)
c
      integer, save, dimension(mxtrcr) ::
     & spcflg ! tmpfn: species type (0=none,1=cold-water,2=warm-water)
c
      real, parameter ::  ! temperature function for cold-water species
     &                    ! thornton and lessem (1978)
     &  theta1 = 16.0,    ! dependence on lower  optimum temperature curve
     &  theta2 =  9.0,    ! dependence on higher optimum temperature curve
     &  theta3 = 11.0,    ! maximum temperature (upper tolerance level)
     &  q10l   =  2.0,    ! the metabolic q10 for temperature response
     &  xk1    =  0.5,    ! scalar constant
     &  xk2    =  0.98,   ! scalar constant
     &  xk3    =  0.01,   ! scalar constant
     &  xk4    =  0.01    ! scalar constant
c
      real, parameter ::  ! temperature function for warm-water species
     &  tmax   = 27.0,    ! Tfunc: maximum tolerated temperature
     &  topt   = 25.0,    ! Tfunc: optimum temperature
     &  q10w   =  2.0     ! Tfunc: the metabolic q10 for temperature response
c
      integer i,j,jrlv,k,l
      real    bm_n,bm_p,bm_z,bm_d,bn_n,bn_p,bn_z,bn_d,
     &        bu_n,bu_p,bu_z,bu_d,
     &        gamma1,gamma2,xnum,xkatheta,ynum,xkbtheta,
     &        tijk,tfn,vw,xw,yw,zw, 
     &        pgrw,zgrw,pref,prf2,qprf,ztgx,dofz,pofz,zofz,
     &        pijk,pijkp,par
c
      if     (ibio.lt.0) then !initialize only
c
c ---   read from tracer.input:
c ---   'biotyp' = type (90X=std.bio,X=3,4,9) must be 904
c
c ---   'pp    ' = zoopl: preference term for phytoplankton
c ---   'pz    ' = zoopl: preference term for zooplankton
c ---   'pd    ' = zoopl: preference term for detritus
c ---   'aa    ' = zoopl: assimilation efficiency
c ---   'am    ' = zoopl: metabolic    efficiency
c ---   'fkz   ' = zoopl: half-saturation coefficient (mg/m^3)
c ---   'gmax  ' = zoopl: maximum growth rate (1/day)
c ---   'zmor  ' = zoopl: mortality (1/day)
c
* ---   'ik    ' = phyto: light absorption efficiency scalar (einst/m^2/h)
c ---   'fkp   ' = phyto: half-saturation coefficient (mg/m^3)
c ---   'pmax  ' = phyto: maximum growth rate (1/day)
c ---   'psen  ' = phyto: senescence (1/day)
c
c ---   'remn  ' = detri: remineralization (1/day)
c
c ---   'spcflg' = tmpfn: species type (0=none,1=cold-water,2=warm-water)
c
        i = -ibio
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3,a,i3,a)')
     &    'Lima/Idrisi NPZD parameters for tracers',i,' to',i+3,':'
        endif !1st tile
c
        call blkini(k, 'biotyp')
        if     (k.ne.904) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')
     &        'error - biotyp must be 904'
          call flush(lp)
          endif !1st tile
          call xcstop('(trcini)')
                 stop '(trcini)'
        endif !biotyp.ne.904
c
        call blkinr(pp(    i), 'pp    ','(a6," =",f10.4," ")')
        call blkinr(pz(    i), 'pz    ','(a6," =",f10.4," ")')
        call blkinr(pd(    i), 'pd    ','(a6," =",f10.4," ")')
        call blkinr(aa(    i), 'aa    ','(a6," =",f10.4," ")')
        call blkinr(am(    i), 'am    ','(a6," =",f10.4," ")')
        call blkinr(fkz(   i), 'fkz   ','(a6," =",f10.4," mg/m^3")')
        call blkinr(gmax(  i), 'gmax  ','(a6," =",f10.4," 1/day")')
        call blkinr(zmor(  i), 'zmor  ','(a6," =",f10.4," 1/day")')
c
        call blkinr(fkp(   i), 'fkp   ','(a6," =",f10.4," mg/m^3")')
        call blkinr(pmax(  i), 'pmax  ','(a6," =",f10.4," 1/day")')
        call blkinr(psen(  i), 'psen  ','(a6," =",f10.4," 1/day")')
c
        call blkinr(remn(  i), 'remn  ','(a6," =",f10.4," 1/day")')
c
        call blkinl(spcflg(i),'spcflg')
c
        if (mnproc.eq.1) then
        write(lp,*)
        endif !1st tile
        return
      endif !ibio.lt.0
c
c --- leapfrog time step.
c
      margin = 0  ! no horizontal derivatives
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,jrlv,pijk,pijkp,par,
!$OMP&                    bm_n,bm_p,bm_z,bm_d,bn_n,bn_p,bn_z,bn_d,
!$OMP&                    bu_n,bu_p,bu_z,bu_d,
!$OMP&                    gamma1,gamma2,xnum,xkatheta,ynum,xkbtheta,
!$OMP&                    tijk,tfn,vw,xw,yw,zw,
!$OMP&                    pgrw,zgrw,pref,prf2,qprf,ztgx,dofz,pofz,zofz)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            jrlv=jerlov(i,j)
            pijkp=0.0
            do k=1,kk
              pijk  = pijkp
              pijkp = pijk+dp(i,j,k,n)
              par   = exp(-0.5*(pijk+pijkp)*betabl(jrlv))
c
              bm_n = tracer(i,j,k,m,ibio)
              bm_p = tracer(i,j,k,m,ibio+1)
              bm_z = tracer(i,j,k,m,ibio+2)
              bm_d = tracer(i,j,k,m,ibio+3)
              bn_n = tracer(i,j,k,n,ibio)
              bn_p = tracer(i,j,k,n,ibio+1)
              bn_z = tracer(i,j,k,n,ibio+2)
              bn_d = tracer(i,j,k,n,ibio+3)
c
              if (spcflg(ibio).eq.1) then
c ---           cold-water species temperature dependance
                tijk     = temp(i,j,k,n)
                gamma1   = 1.0/(theta2-q10l) *
     &                     log((xk2*(1.0-xk1))/(xk1*(1.0-xk2)))
                gamma2   = 1.0/(theta1-theta3) *
     &                     log((xk2*(1.0-xk3))/(xk4*(1.0-xk2)))
                xnum     = exp(gamma1*(tijk-q10l))
                xkatheta = (xk1*xnum)/(1.0+xk1*(xnum-1.0))
                ynum     = exp(gamma2*(theta1-tijk))
                xkbtheta = (xk4*ynum)/(1.0+xk3*(ynum-1.0))
                tfn      = xkatheta*xkbtheta
              elseif (spcflg(ibio).eq.2) then
c ---           warm-water species temperature dependance
                tijk     = temp(i,j,k,n)
                if (tijk.le.tmax) then
                  vw  = (tmax-tijk)/(tmax-topt)
                  yw  = log(q10w)*(tmax-topt+2.0)
                  zw  = log(q10w)*(tmax-topt)
                  xw  = (zw**2 * (1.0+sqrt(1.0+40.0/yw))**2)/400.0
                  tfn = vw**xw * exp(xw*(1.0-vw))
                else
                  tfn=0.0
                endif
              else
c ---           no temperature dependance
                tfn=1.0
              endif !spcflg
c
              pref = pp(ibio)*bm_p +
     &               pd(ibio)*bm_d +
     &               pz(ibio)*bm_z
              prf2 = pp(ibio)*bm_p**2 +
     &               pd(ibio)*bm_d**2 +
     &               pz(ibio)*bm_z**2
              qprf = 1.0/(fkz(ibio)*pref + prf2 + epsil)  !epsil prevents 1/0
              ztgx = bm_z*tfn*gmax(ibio)
c
              pgrw = bm_p*tfn*pmax(ibio)*bm_n*par/(fkp(ibio)+bm_n)
              zgrw = ztgx*(prf2            *qprf)*aa(ibio)*am(ibio)
              pofz = ztgx*(pp(ibio)*bm_p**2*qprf)
              zofz = ztgx*(pz(ibio)*bm_z**2*qprf)
              dofz = ztgx*(pd(ibio)*bm_d**2*qprf)
c
              ! limit negative terms to 10% of total per single time step
              pgrw = min(pgrw,bn_n*0.2*86400.0/delt1)
              zgrw = min(zgrw,bn_n*0.2*86400.0/delt1)
              pofz = min(pofz,bn_p*0.2*86400.0/delt1)
              zofz = min(zofz,bn_z*0.2*86400.0/delt1)
              dofz = min(dofz,bn_d*0.2*86400.0/delt1)
c
              bu_p =   pgrw
     &               - pofz
     &               - bm_p*psen(ibio)
              bu_z =   zgrw
     &               - zofz
     &               - bm_z*zmor(ibio)
              bu_d =   bm_p*psen(ibio)
     &               + bm_z*zmor(ibio)
     &               + (pofz+zofz+dofz)*(1.0-aa(ibio))
     &               - dofz
     &               - bm_d*remn(ibio)
              bu_n =   bm_d*remn(ibio)
     &               + (pofz+zofz+dofz)*     aa(ibio)
     &               - zgrw
     &               - pgrw
c
              tracer(i,j,k,n,ibio)   = bn_n + delt1/86400.0 * bu_n
              tracer(i,j,k,n,ibio+1) = bn_p + delt1/86400.0 * bu_p
              tracer(i,j,k,n,ibio+2) = bn_z + delt1/86400.0 * bu_z
              tracer(i,j,k,n,ibio+3) = bn_d + delt1/86400.0 * bu_d
            enddo
          enddo !i
        enddo !l
      enddo !j
      return
      end subroutine trcupd_904
c
c
c> Revision history:
c>
c> Aug  2002 - new routine to put all tracer interactions in one place
