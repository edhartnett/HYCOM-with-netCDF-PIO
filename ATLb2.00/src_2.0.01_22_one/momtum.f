      subroutine momtum(m,n)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
      logical    lpipe_momtum
      parameter (lpipe_momtum=.false.)
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &                 stress,stresx,stresy,dpmx,visu,visv,vort,oneta,
     &                 wgtia,wgtib,wgtja,wgtjb,
     &                 dl2u,dl2uja,dl2ujb,dl2v,dl2via,dl2vib
      common/momtumr4/ stress,stresx,stresy,dpmx,visu,visv,vort,oneta,
     &                 wgtia,wgtib,wgtja,wgtjb,
     &                 dl2u,dl2uja,dl2ujb,dl2v,dl2via,dl2vib
      save  /momtumr4/
c
      integer ifirst
      save    ifirst
c
      real dpia,dpib,dpja,dpjb,visca,viscb,scuya,scuyb,scvxa,scvxb,
     &     dpxy,ptopl,pbotl,cutoff,q,
     &     dt1inv,phi,plo,ubot,vbot,thkbop,dpup,dpdn,pstres
      integer i,ia,ib,j,ja,jb,k,ka,l,mbdy
c
*     real*8    wtime
*     external  wtime
*     real*8    wtime1(10),wtime2(20,kdm),wtimes
c
      character text*12
      integer, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & mask
c
      real hfharm,a,b
      include 'stmt_fns.h'
c
c --- harmonic mean divided by 2
      hfharm(a,b)=a*b/(a+b)
c
      data ifirst / 0 /
c
      mbdy = 6
c
      call xctilr(montg(  1-nbdy,1-nbdy,1  ),1,kk, 6,6, halo_ps)
      call xctilr(pu(     1-nbdy,1-nbdy,2  ),1,kk, 6,6, halo_us)
      call xctilr(pv(     1-nbdy,1-nbdy,2  ),1,kk, 6,6, halo_vs)
      call xctilr(dpmixl( 1-nbdy,1-nbdy,  m),1, 1, 6,6, halo_ps)
      call xctilr(dp(     1-nbdy,1-nbdy,1,m),1,kk, 6,6, halo_ps)
      call xctilr(dp(     1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_ps)
      call xctilr(dpu(    1-nbdy,1-nbdy,1,m),1,kk, 6,6, halo_us)
      call xctilr(dpu(    1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_us)
      call xctilr(dpv(    1-nbdy,1-nbdy,1,m),1,kk, 6,6, halo_vs)
      call xctilr(dpv(    1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_vs)
      call xctilr(u(      1-nbdy,1-nbdy,1,m),1,kk, 6,6, halo_uv)
      call xctilr(u(      1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_uv)
      call xctilr(v(      1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_vv)
      call xctilr(v(      1-nbdy,1-nbdy,1,m),1,kk, 6,6, halo_vv)
      call xctilr(ubavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_vv)
      call xctilr(pbavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_ps)
      call xctilr(dpoldm( 1-nbdy,1-nbdy,1  ),1,kk, 6,6, halo_ps)
      call xctilr(saln(   1-nbdy,1-nbdy,1,m),1,kk, 6,6, halo_ps)
      call xctilr(temp(   1-nbdy,1-nbdy,1,m),1,kk, 6,6, halo_ps)
      call xctilr(th3d(   1-nbdy,1-nbdy,1,m),1,kk, 6,6, halo_ps)
      call xctilr(uflx(   1-nbdy,1-nbdy,1  ),1,kk, 6,6, halo_uv)
      call xctilr(vflx(   1-nbdy,1-nbdy,1  ),1,kk, 6,6, halo_vv)
c
      if     (ifirst.eq.0) then
        ifirst=1
c ---   setup zero fill.
        margin = mbdy
c
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            visu(i,j)=0.0
            visv(i,j)=0.0
            dl2u(i,j)=0.0
            dl2v(i,j)=0.0
          enddo
        enddo
      endif
c
c --- --------------------
c --- hydrostatic equation
c --- --------------------
c
*        wtime1( 1) = wtime()
c
c --- rhs: th3d.m, temp.m, saln.m, p, pbavg.m
c --- lhs: thstar, p, oneta, montg
c
      margin = mbdy
c
!$OMP PARALLEL DO PRIVATE(j,l,k,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do k=1,kk
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              if (tbaric) then
c
c ---   sigma-star is virtual potential density, as defined in 
c ---   Sun et.al. (1999), 'Inclusion of thermobaricity in 
c ---   isopycnic-coordinate ocean models', JPO 29 pp 2719-2729.
c
c ---   use upper interface pressure in converting sigma to sigma-star.
c ---   to avoid density variations in layers intersected by bottom
c
                thstar(i,j,k)=th3d(i,j,k,m)
     &             +kappaf(temp(i,j,k,m),saln(i,j,k,m),p(i,j,k))
              else
                thstar(i,j,k)=th3d(i,j,k,m)
              endif
c
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,m)
            enddo
          enddo
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
c ---       store (1+eta) (= p_total/p_prime) in -oneta-
            oneta(i,j)=1.+pbavg(i,j,m)/p(i,j,kk+1)
c
c ---       m_prime in lowest layer:
            montg(i,j,kk)=psikk(i,j)+
     &        ( p(i,j,kk+1)*(thkk(i,j)-thstar(i,j,kk))
     &          -pbavg(i,j,m)*(thstar(i,j,kk)+thbase) )*thref**2
          enddo
c
c ---     m_prime in remaining layers:
          do k=kk-1,1,-1
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              montg(i,j,k)=montg(i,j,k+1)+p(i,j,k+1)*oneta(i,j)
     &                  *(thstar(i,j,k+1)-thstar(i,j,k))*thref**2
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
*        wtime1( 2) = wtime()
      call dpudpv(p,depthu,depthv,dpu(1-nbdy,1-nbdy,1,m),
     &                            dpv(1-nbdy,1-nbdy,1,m))
*        wtime1( 3) = wtime()
c
c --- account for temporal smoothing of mid-time dpmixl. calculate the vertical
c --- excursions of the coordinates immediately above and below the mixed
c --- layer base, then vertically interpolate this motion to dpmixl(i,j,m)
c
      if(hybrid .and. mxlkta .and. thermo) then
c
c ---   rhs: dpoldm, dpmixl.m
c ---   lhs: util1, util2
c
        margin = mbdy
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,k,dpup,dpdn,q)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              util1(i,j)=0.
              util2(i,j)=0.
            enddo
            do k=1,kk
              do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                util1(i,j)=util2(i,j)
                util2(i,j)=util2(i,j)+dpoldm(i,j,k)
                if     (util2(i,j).ge.dpmixl(i,j,m).and.
     &                  util1(i,j).lt.dpmixl(i,j,m)     ) then
                  dpup=p(i,j,k  )-util1(i,j)
                  dpdn=p(i,j,k+1)-util2(i,j)
                  q=(util2(i,j)-dpmixl(i,j,m))/max(onemm,dpoldm(i,j,k))
                  dpmixl(i,j,m)=dpmixl(i,j,m)+(dpdn+q*(dpup-dpdn))
                endif
              enddo
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
      endif
c
c +++ ++++++++++++++++++
c +++ momentum equations
c +++ ++++++++++++++++++
c
c --- bottom drag (standard bulk formula)
c
      thkbop=thkbot*g/thref
c
*        wtime1( 4) = wtime()
c
c --- rhs: p, u.n+, v.n+, ubavg.n+, vbavg.n+, depthv+, pvtrop+
c --- rhs: dpmixl.m+, taux+, dpu, depthu+, dpv, tauy+
c --- lhs: util1, util2, drag, ubrhs, stresx, vbrhs, stresy
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,phi,plo,ubot,vbot,pstres)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            util1(i,j)=0.
            util2(i,j)=0.
          enddo
          do k=1,kk
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              phi=max(p(i,j,k+1),p(i,j,kk+1)-thkbop)
              plo=max(p(i,j,k  ),p(i,j,kk+1)-thkbop)
              util1(i,j)=util1(i,j)+
     &                     (u(i,j,k,n)+u(i+1,j,k,n))*(phi-plo)
              util2(i,j)=util2(i,j)+
     &                     (v(i,j,k,n)+v(i,j+1,k,n))*(phi-plo)
            enddo
          enddo
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            ubot=ubavg(i,j,n)+ubavg(i+1,j,n)+util1(i,j)/thkbop
            vbot=vbavg(i,j,n)+vbavg(i,j+1,n)+util2(i,j)/thkbop
            drag(i,j)=cb*(.25*sqrt(ubot*ubot+vbot*vbot)+cbar)/thkbot
          enddo
        enddo
c
c ---   store r.h.s. of barotropic u/v eqn. in -ubrhs,vbrhs-
c ---   time-interpolate wind stress
c
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            ubrhs(i,j)=
     &        (vbavg(i  ,j,  m)*depthv(i  ,j)
     &        +vbavg(i  ,j+1,m)*depthv(i  ,j+1)
     &        +vbavg(i-1,j,  m)*depthv(i-1,j)
     &        +vbavg(i-1,j+1,m)*depthv(i-1,j+1))
     &        *(pvtrop(i,j)+pvtrop(i,j+1))*.125
c
            if     (wndflg.eq.1) then  ! taux on u grid
              if(hybrid .and. mxlkrt) then
                pstres=.5*(dpmixl(i,j,m)+dpmixl(i-1,j,m))
              else
                pstres=dpu(i,j,1,m)
              endif
              stresx(i,j)=(taux(i,j,l0)*w0
     &                    +taux(i,j,l1)*w1
     &                    +taux(i,j,l2)*w2
     &                    +taux(i,j,l3)*w3)*g/(pstres*thref)
            elseif (wndflg.eq.2) then  ! taux on p grid
              if(hybrid .and. mxlkrt) then
                pstres=.5*(dpmixl(i,j,m)+dpmixl(i-1,j,m))
              else
                pstres=dpu(i,j,1,m)
              endif
              stresx(i,j)=((taux(i,j,l0)+taux(i-1,j,l0))*w0
     &                    +(taux(i,j,l1)+taux(i-1,j,l1))*w1
     &                    +(taux(i,j,l2)+taux(i-1,j,l2))*w2
     &                    +(taux(i,j,l3)+taux(i-1,j,l3))*w3)*.5*g
     &                    /(pstres*thref)
            else  ! no taux
              stresx(i,j)=0.
            endif
          enddo
        enddo
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vbrhs(i,j)=
     &      -(ubavg(i,  j  ,m)*depthu(i,j  )
     &       +ubavg(i+1,j  ,m)*depthu(i+1,j  )
     &       +ubavg(i,  j-1,m)*depthu(i,j-1)
     &       +ubavg(i+1,j-1,m)*depthu(i+1,j-1))
     &       *(pvtrop(i,j)+pvtrop(i+1,j))*.125
c
            if     (wndflg.eq.1) then  ! tauy on v grid
              if(hybrid .and. mxlkrt) then
                pstres=.5*(dpmixl(i,j,m)+dpmixl(i,j-1,m))
              else
                pstres=dpv(i,j,1,m)
              endif
              stresy(i,j)=(tauy(i,j,l0)*w0
     &                    +tauy(i,j,l1)*w1
     &                    +tauy(i,j,l2)*w2
     &                    +tauy(i,j,l3)*w3)*g/(pstres*thref)
            elseif (wndflg.eq.2) then  ! tauy on p grid
              if(hybrid .and. mxlkrt) then
                pstres=.5*(dpmixl(i,j,m)+dpmixl(i,j-1,m))
              else
                pstres=dpv(i,j,1,m)
              endif
              stresy(i,j)=((tauy(i,j,l0)+tauy(i,j-1,l0))*w0
     &                    +(tauy(i,j,l1)+tauy(i,j-1,l1))*w1
     &                    +(tauy(i,j,l2)+tauy(i,j-1,l2))*w2
     &                    +(tauy(i,j,l3)+tauy(i,j-1,l3))*w3)*.5*g
     &                    /(pstres*thref)
            else  ! no tauy
              stresy(i,j)=0.
            endif
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
c --- the old  momeq2.f  starts here
c
      cutoff=.5*onem
c
*        wtime1( 5) = wtime()
c
c --- rhs: 0.0
c --- lhs: util1, util2
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
c ---   spatial weighting function for pressure gradient calculation:
          util1(i,j)=0.
          util2(i,j)=0.
        enddo
      enddo
c
      do 9 k=1,kk
c
c --- store total (barotropic plus baroclinic) flow at old and mid time in
c --- -utotn,vtotn- and -utotm,vtotm- respectively. store minimum thickness
c --- values for use in pot.vort. calculation in -dpmx-.
c
*         wtime2( 1,k) = wtime()
c
c --- rhs: dpmx, dp.m+
c --- lhs: dpmx
c
*     margin = mbdy - 2
c
      do i=1-margin,ii+margin
        dpmx(i,1)=2.*cutoff
      enddo
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          dpmx(i,j+1)=2.*cutoff
        enddo
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dpmx(i,j+1)=max(dpmx(i,j+1),dp(i,j,k,m)+dp(i-1,j,k,m))
          enddo
        enddo
      enddo
c
*         wtime2( 2,k) = wtime()
c
c --- rhs: ubavg.m, ubavg.n, dp.m+, dpu
c --- lhs: utotm, utotn, uflux, dpmx, pu
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
c
          i=ifu(j,l)-1
          if (iuopn(i,j).ne.0) then
*           uflux(i,j)=uflx(i,j,k)*scuxi(i,j)
*           utotm(i,j)=u(i,j,k,m)+ubavg(i,j,m)
*           utotn(i,j)=u(i,j,k,n)+ubavg(i,j,n)
            utotm(i,j)=ubavg(i,j,m)
            utotn(i,j)=ubavg(i,j,n)
            uflux(i,j)=utotm(i,j)*max(dp(i,j,k,m),cutoff)
          endif
          i=ilu(j,l)+1
          if (iuopn(i,j).ne.0) then
*           uflux(i,j)=uflx(i,j,k)*scuxi(i,j)
*           utotm(i,j)=u(i,j,k,m)+ubavg(i,j,m)
*           utotn(i,j)=u(i,j,k,n)+ubavg(i,j,n)
            utotm(i,j)=ubavg(i,j,m)
            utotn(i,j)=ubavg(i,j,n)
            uflux(i,j)=utotm(i,j)*max(dp(i-1,j,k,m),cutoff)
          endif
c
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dpmx(i,j)=max(dpmx(i,j),dp(i,j,k,m)+dp(i-1,j,k,m))
            utotm(i,j)=u(i,j,k,m)+ubavg(i,j,m)
            utotn(i,j)=u(i,j,k,n)+ubavg(i,j,n)
            uflux(i,j)=utotm(i,j)*max(dpu(i,j,k,m),cutoff)
            pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,k,m)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
*         wtime2( 3,k) = wtime()
c
c --- rhs: vbavg.m, vbavg.n, dp.m+, dpv
c --- lhs: vtotm, vtotn, vflux, dpmx, pv
c
*     margin = mbdy - 2
c
      do i=1-margin,ii+margin
        do l=1,jsv(i)
        j=jfv(i,l)-1
        if (ivopn(i,j).ne.0) then
*         vflux(i,j)=vflx(i,j,k)*scvyi(i,j)
*         vtotm(i,j)=v(i,j,k,m)+vbavg(i,j,m)
*         vtotn(i,j)=v(i,j,k,n)+vbavg(i,j,n)
          vtotm(i,j)=vbavg(i,j,m)
          vtotn(i,j)=vbavg(i,j,n)
          vflux(i,j)=vtotm(i,j)*max(dp(i,j,k,m),cutoff)
        endif
        j=jlv(i,l)+1
        if (ivopn(i,j).ne.0) then
*         vflux(i,j)=vflx(i,j,k)*scvyi(i,j)
*         vtotm(i,j)=v(i,j,k,m)+vbavg(i,j,m)
*         vtotn(i,j)=v(i,j,k,n)+vbavg(i,j,n)
          vtotm(i,j)=vbavg(i,j,m)
          vtotn(i,j)=vbavg(i,j,n)
          vflux(i,j)=vtotm(i,j)*max(dp(i,j-1,k,m),cutoff)
        endif
        enddo
      enddo
c
*         wtime2( 4,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            dpmx(i  ,j)=max(dpmx(i  ,j),dp(i,j,k,m)+dp(i,j-1,k,m))
            dpmx(i+1,j)=max(dpmx(i+1,j),dp(i,j,k,m)+dp(i,j-1,k,m))
            vtotm(i,j)=v(i,j,k,m)+vbavg(i,j,m)
            vtotn(i,j)=v(i,j,k,n)+vbavg(i,j,n)
            vflux(i,j)=vtotm(i,j)*max(dpv(i,j,k,m),cutoff)
            pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k,m)
          enddo
        enddo
      enddo
c
c --- define auxiliary velocity fields (via,vib,uja,ujb) to implement
c --- sidewall friction along near-vertical bottom slopes. wgtja,wgtjb,wgtia,
c --- wgtib indicate the extent to which a sidewall is present.
c
*         wtime2( 5,k) = wtime()
c
c --- rhs: pu, depthu+, utotn+, wgtja
c --- lhs: wgtja, wgtjb, uja, ujb, dl2u
c --- rhs: pv, depthv+, vtotn+, wgtia
c --- lhs: wgtia, wgtib, via, vib, dl2v
c --- rhs: vtotm, vort+, corio+, dp.m+, dpmx+, vtotn
c --- lhs: vort, potvor, defor2
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,jb,l,i,ia,ib)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c ---   assume margin<nblk
        ja=j-1
        jb=j+1
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            wgtja(i,j)=max(0.,min(1.,(pu(i,j,k+1)-depthu(i,ja))
     &                /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
            wgtjb(i,j)=max(0.,min(1.,(pu(i,j,k+1)-depthu(i,jb))
     &                /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
            uja(i,j)=(1.-wgtja(i,j))*utotn(i,ja)+
     &               wgtja(i,j)*slip*utotn(i,j)
            ujb(i,j)=(1.-wgtjb(i,j))*utotn(i,jb)+
     &               wgtjb(i,j)*slip*utotn(i,j)
            dl2u(i,j)=utotn(i,j)
     &          -biharm*.25*(utotn(i+1,j)+utotn(i-1,j)+
     &                         uja(i,  j)+  ujb(i,  j))
          enddo
        enddo
c
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
c ---       assume margin<nblk
            ia=i-1
            ib=i+1
            wgtia(i,j)=max(0.,min(1.,(pv(i,j,k+1)-depthv(ia,j))
     &                /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
            wgtib(i,j)=max(0.,min(1.,(pv(i,j,k+1)-depthv(ib,j))
     &                /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
            via(i,j)=(1.-wgtia(i,j))*vtotn(ia,j)+
     &                wgtia(i,j)*slip*vtotn(i,j)
            vib(i,j)=(1.-wgtib(i,j))*vtotn(ib,j)+
     &                wgtib(i,j)*slip*vtotn(i,j)
            dl2v(i,j)=vtotn(i,j)
     &          -biharm*.25*(vtotn(i,j+1)+vtotn(i,j-1)+
     &                         via(i,j  )+  vib(i,j  ))
          enddo
        enddo
c
c --- vorticity, pot.vort., defor. at lateral boundary points
        do l=1,isv(j)
          i=ifv(j,l)
          vort(i  ,j)= vtotm(i,j)*(1.-slip)*scvy(i,j)*scq2i(i  ,j)
          potvor(i  ,j)=(vort(i  ,j)+corio(i  ,j)) * 8.
     &    /max(8.*cutoff,
     &         4.*(dp(i,j,k,m)+dp(i,ja ,k,m)),
     &         dpmx(i,j),
     &         dpmx(i+1,j))
          defor2(i  ,j)=(vtotn(i,j)*(1.-slip)*scvy(i,j))**2*scq2i(i  ,j)
          i=ilv(j,l)
          vort(i+1,j)=-vtotm(i,j)*(1.-slip)*scvy(i,j)*scq2i(i+1,j)
          potvor(i+1,j)=(vort(i+1,j)+corio(i+1,j)) * 8.
     &    /max(8.*cutoff,
     &         4.*(dp(i,j,k,m)+dp(i,ja ,k,m)),
     &         dpmx(i,j),
     &         dpmx(i+1,j))
          defor2(i+1,j)=(vtotn(i,j)*(1.-slip)*scvy(i,j))**2*scq2i(i+1,j)
        enddo
      enddo
!$OMP END PARALLEL DO
c
*         wtime2( 6,k) = wtime()
c
c --- vorticity, pot.vort., defor. at lateral boundary points
c
c --- rhs: utotm, vort+, corio+, dp.m+, dpmx+, utotn
c --- lhs: vort, potvor, defor2
c
*     margin = mbdy - 2
c
      do i=1-margin,ii+margin
c ---   assume margin<nblk
        ia=i-1
        do l=1,jsu(i)
          j=jfu(i,l)
          vort(i,j  )=-utotm(i,j)*(1.-slip)*scux(i,j)*scq2i(i,j  )
          potvor(i,j  )=(vort(i,j  )+corio(i,j  )) * 8.
     &    /max(8.*cutoff,
     &         4.*(dp(i,j,k,m)+dp(ia ,j,k,m)),
     &         dpmx(i,j),
     &         dpmx(i,j+1))
          defor2(i,j  )=(utotn(i,j)*(1.-slip)*scux(i,j))**2*scq2i(i,j  )
          j=jlu(i,l)
          vort(i,j+1)= utotm(i,j)*(1.-slip)*scux(i,j)*scq2i(i,j+1)
          potvor(i,j+1)=(vort(i,j+1)+corio(i,j+1)) * 8.
     &    /max(8.*cutoff,
     &         4.*(dp(i,j,k,m)+dp(ia ,j,k,m)),
     &         dpmx(i,j),
     &         dpmx(i,j+1))
          defor2(i,j+1)=(utotn(i,j)*(1.-slip)*scux(i,j))**2*scq2i(i,j+1)
        enddo
      enddo
c
*         wtime2( 7,k) = wtime()
c
c --- rhs: p, utotn+, vtotn+
c --- lhs: util3, defor1
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            util3(i,j)=.5*(p(i,j,k+1)+p(i,j,k))*oneta(i,j)
            defor1(i,j)=((utotn(i+1,j)*scuy(i+1,j)
     &                   -utotn(i,  j)*scuy(i,  j))
     &                  -(vtotn(i,j+1)*scvx(i,j+1)
     &                   -vtotn(i,j  )*scvx(i,j  )))**2
     &                  *scp2i(i,j)
          enddo
        enddo
      enddo
c
c --- vorticity, pot.vort., defor. at interior points (incl. promontories)
*         wtime2( 8,k) = wtime()
c
c --- rhs: vtotm+, utotm+, vort, dp.m+, dpmx+, vib+, via, ujb+, uja
c --- lhs: vort, potvor, defor2
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isq(j)
          do i=max(1-margin,ifq(j,l)),min(ii+margin,ilq(j,l))
            vort(i,j)=(vtotm(i,j)*scvy(i,j)-vtotm(i-1,j)*scvy(i-1,j)
     &                -utotm(i,j)*scux(i,j)+utotm(i,j-1)*scux(i,j-1))
     &                *scq2i(i,j)
            potvor(i,j)=(vort(i,j)+corio(i,j)) * 8.
     &         /max(8.*cutoff,
     &              2.*(dp(i,j  ,k,m)+dp(i-1,j  ,k,m)+
     &                  dp(i,j-1,k,m)+dp(i-1,j-1,k,m)),
     &              dpmx(i,j),dpmx(i-1,j),dpmx(i+1,j),
     &                        dpmx(i,j-1),dpmx(i,j+1))
            defor2(i,j)=(vib(i-1,j)*scvy(i,j)-via(i,j)*scvy(i-1,j)
     &                  +ujb(i,j-1)*scux(i,j)-uja(i,j)*scux(i,j-1))**2
     &                  *scq2i(i,j)
          enddo
        enddo
      enddo
c
c --- define auxiliary del2 fields (dl2via,dl2vib,dl2uja,dl2ujb) to imple-
c --- ment biharmonic sidewall friction along near-vertical bottom slopes.
c
*         wtime2( 9,k) = wtime()
c
c --- rhs: wgtja, wgtjb, dlu2+, wgtia, wgtib, dlv2+
c --- lhs: dl2uja, dl2ujb, dl2via, dl2vib
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,jb,l,i,ia,ib)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        ja=j-1
        jb=j+1
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dl2uja(i,j)=(1.-wgtja(i,j))*dl2u(i,ja)+
     &                  wgtja(i,j)*slip*dl2u(i,j)
            dl2ujb(i,j)=(1.-wgtjb(i,j))*dl2u(i,jb)+
     &                  wgtjb(i,j)*slip*dl2u(i,j)
          enddo
        enddo
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ia=i-1
            ib=i+1
            dl2via(i,j)=(1.-wgtia(i,j))*dl2v(ia,j)+
     &                  wgtia(i,j)*slip*dl2v(i, j)
            dl2vib(i,j)=(1.-wgtib(i,j))*dl2v(ib,j)+
     &                  wgtib(i,j)*slip*dl2v(i, j)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
*         wtime2(10,k) = wtime()
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        do j=1,jj
          do i=1,ii
            mask(i,j)=min(1,iq(i,j)+iu(i,j  )+iv(i,j  )
     &                             +iu(i,j-1)+iv(i-1,j))
          enddo
        enddo
        write (text,'(a9,i3)') 'potvor k=',k
        call pipe_compare(potvor,mask,text)
      endif
c
c --- ----------
c --- u equation
c --- ----------
c
c --- deformation-dependent eddy viscosity coefficient
c
*         wtime2(11,k) = wtime()
c
c --- rhs: defor1+, defor2+
c --- lhs: visu
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            visu(i,j)=max(veldff,
     &                    viscos*sqrt(.5*(defor1(i,j)+defor1(i-1,j)+
     &                                    defor2(i,j)+defor2(i,j+1))))
          enddo
        enddo
      enddo
c
*         wtime2(12,k) = wtime()
c
c --- rhs: visu+, dl2u+, dpu.m+, dl2uja, wgtja,, dl2ujb, wgtjb
c --- lhs: visu, uflux1, uflux2, uflux3
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,jb,l,i,
!$OMP&                    dpxy,dpja,dpjb,visca,viscb,scuya,scuyb)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        ja=j-1
        jb=j+1
        do l=1,isu(j)
          visu(ifu(j,l)-1,j)=visu(ifu(j,l),j)
          visu(ilu(j,l)+1,j)=visu(ilu(j,l),j)
c
c ---     longitudinal turb. momentum flux (at mass points)
c
          do i=max(1-margin,ifu(j,l)-1),min(ii+margin,ilu(j,l))
            uflux1(i,j)=(visu(i,j)+visu(i+1,j))*(dl2u(i,j)-dl2u(i+1,j))
     &                   *hfharm(max(dpu(i  ,j,k,m),onemm),
     &                           max(dpu(i+1,j,k,m),onemm))
     &                   *scp2(i,j)*2./(scux(i,j)+scux(i+1,j))
          enddo
c
c --- lateral turb. momentum flux (at vorticity points)
c --- (left and right fluxes are evaluated separately because of sidewalls)
c
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dpxy=max(dpu(i,j ,k,m),onemm)
            dpja=max(dpu(i,ja,k,m),onemm)
            dpjb=max(dpu(i,jb,k,m),onemm)
c
c --- check whether variables along coast have been initialized correctly
cdiag if (k.eq.kk) then
cdiag   if (iu(i,ja).eq.0 .and. dpu(i,ja,k,m).ne.0.) then
cdiag!$OMP CRITICAL
cdiag     write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.      '  error - nonzero dpu(ja):',dpu(i,ja,k,m)
cdiag!$OMP END CRITICAL
cdiag   endif
cdiag   if (iu(i,jb).eq.0 .and. dpu(i,jb,k,m).ne.0.) then
cdiag!$OMP CRITICAL
cdiag     write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.      '  error - nonzero dpu(jb):',dpu(i,jb,k,m)
cdiag!$OMP END CRITICAL
cdiag   endif
cdiag endif
c
            if (iu(i,ja).eq.0) then
              visca=visu(i,j )
              scuya=scuy(i,j )
            else
              visca=visu(i,ja)
              scuya=scuy(i,ja)
            endif
            if (iu(i,jb).eq.0) then
              viscb=visu(i,j )
              scuyb=scuy(i,j )
            else
              viscb=visu(i,jb)
              scuyb=scuy(i,jb)
            endif
            uflux2(i,j)=(visu(i,j)+visca)*(dl2uja(i,j)-dl2u(i,j))
     &                  *hfharm(dpja+wgtja(i,j)*(dpxy-dpja),dpxy)
     &                  *scq2(i,j )*2./(scuy(i,j)+scuya)
            uflux3(i,j)=(visu(i,j)+viscb)*(dl2u(i,j)-dl2ujb(i,j))
     &                  *hfharm(dpjb+wgtjb(i,j)*(dpxy-dpjb),dpxy)
     &                  *scq2(i,jb)*2./(scuy(i,j)+scuyb)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
c --- pressure force in x direction
c --- ('scheme 2' from appendix -a- in bleck-smith paper)
c
*         wtime2(13,k) = wtime()
c
c --- rhs: depthu, pu, montg+, thstar+, p+, dp.m+
c --- lhs: util1, pgfx
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            util1(i,j)=max(0.,min(depthu(i,j)-pu(i,j,k),h1))
            pgfx(i,j)=util1(i,j)*
     &          (montg(i,j,k)-montg(i-1,j,k)+
     &           (thstar(i,j,k)-thstar(i-1,j,k))
     &          *(p(i,j,k+1)*p(i-1,j,k+1)-
     &            p(i,j,k  )*p(i-1,j,k  ))*thref**2
     &            /(dp(i,j,k,m)+dp(i-1,j,k,m)+epsil))
          enddo
        enddo
      enddo
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'pu     k=',k
        call pipe_compare(pu(1-nbdy,1-nbdy,k),iu,text)
        write (text,'(a9,i3)') 'depthu k=',k
        call pipe_compare(depthu,iu,text)
        write (text,'(a9,i3)') 'util1  k=',k
        call pipe_compare(util1, iu,text)
        util4(1:ii,1:jj) = montg( 1:ii,1:jj,k)-montg( 0:ii-1,1:jj,k)
        write (text,'(a9,i3)') 'montgX k=',k
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = thstar(1:ii,1:jj,k)-thstar(0:ii-1,1:jj,k)
        write (text,'(a9,i3)') 'thstaX k=',k
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k)*p(0:ii-1,1:jj,k)
        write (text,'(a9,i3)') 'pSQ    k=',k
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k+1)*p(0:ii-1,1:jj,k+1)
        write (text,'(a9,i3)') 'pSQ    k=',k+1
        call pipe_compare(util4, iu,text)
        util4(1:ii,1:jj) = dp(1:ii,1:jj,k,m)+dp(0:ii-1,1:jj,k,m)+epsil
        write (text,'(a9,i3)') 'dp.m+  k=',k
        call pipe_compare(util4, iu,text)
        write (text,'(a9,i3)') 'pgfx   k=',k
        call pipe_compare(pgfx,  iu,text)
      endif
c
cdiag if (k.eq.1) then
cdiag write (lp,100) nstep
cdiag do j=jtest-1,jtest+1
c
*         wtime2(14,k) = wtime()
c
c --- rhs: pgfx+, util1+, p+, dpmixl.m+, dpu.m, utotn, drag+
c --- rhs: stresx, u.m, u.n, dpu.n, gradx, utotm+, vtotm+,
c --- rhs: vflux+, potvor+, uflux1+, uflux2, uflux3
c --- lhs: gradx, stress, u.m, u.n
c
*     margin = mbdy - 2
c
!$OMP PARALLEL DO PRIVATE(j,ja,jb,l,i,ptopl,pbotl,pstres)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        ja=j-1
        jb=j+1
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
c
c --- check whether variables along coast have been initialized correctly
cdiag       if (k.eq.kk) then
cdiag         if (iu(i,ja).eq.0 .and. pgfx(i,ja).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.            '  error - nonzero pgfx(ja):',pgfx(i,ja)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag         if (iu(i,jb).eq.0 .and. pgfx(i,jb).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.            '  error - nonzero pgfx(jb):',pgfx(i,jb)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag       endif
c
            gradx(i,j)=(pgfx(i,j)+(h1-util1(i,j))*
     &        (pgfx (i-1,j)+pgfx (i+1,j)+pgfx (i,ja)+pgfx (i,jb))/
     &        (util1(i-1,j)+util1(i+1,j)+util1(i,ja)+util1(i,jb)+
     &          epsil))/h1
c
            ptopl=min(depthu(i,j),.5*(p(i,j,k  )+p(i-1,j,k  )))
            pbotl=min(depthu(i,j),.5*(p(i,j,k+1)+p(i-1,j,k+1)))
            if(hybrid .and. mxlkrt) then
              pstres=.5*(dpmixl(i,j,m)+dpmixl(i-1,j,m))
            else
              pstres=dpu(i,j,1,m)
            endif
c
c ---       top and bottom boundary layer stress
            stress(i,j)=(-utotn(i,j)*.5*(drag(i,j)+drag(i-1,j))*
     &           (max(depthu(i,j)-thkbop,          pbotl       )
     &           -max(depthu(i,j)-thkbop,min(ptopl,pbotl-onemm)))
     &           +stresx(i,j)*thref*(min(pstres,pbotl+onemm)
     &                              -min(pstres,ptopl      )))
     &           /max(dpu(i,j,k,m),onemm)
c
c ---       time smoothing of -u- field  (part 1)
            u(i,j,k,m)=u(i,j,k,m)*(wuv1*dpu(i,j,k,m)+onemm)
     &                +u(i,j,k,n)* wuv2*dpu(i,j,k,n)
c
cdiag       util4(i,j)=u(i,j,k,n)
            u(i,j,k,n)=u(i,j,k,n)+delt1*(-scuxi(i,j)*(gradx(i,j)
     &      +.25*(utotm(i+1,j)**2+vtotm(i  ,j)**2+vtotm(i  ,j+1)**2
     &           -utotm(i-1,j)**2-vtotm(i-1,j)**2-vtotm(i-1,j+1)**2))
     &      +.125*(vflux(i  ,j)+vflux(i  ,j+1)+
     &             vflux(i-1,j)+vflux(i-1,j+1))
     &           *(potvor(i,j)+potvor(i,j+1))
     &      - ubrhs(i,j) + stress(i,j)
     &      -(uflux1(i,j)-uflux1(i-1,j)
     &       +uflux3(i,j)-uflux2(i  ,j))
     &        /(scu2(i,j)*max(dpu(i,j,k,m),onemm)))
          enddo
        enddo
      enddo
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'vtotm  k=',k
        call pipe_compare(vtotm, iv,text)
        write (text,'(a9,i3)') 'utotm  k=',k
        call pipe_compare(utotm, iu,text)
        write (text,'(a9,i3)') 'gradx  k=',k
        call pipe_compare(gradx, iu,text)
        write (text,'(a9,i3)') 'ubrhs  k=',k
        call pipe_compare(ubrhs, iu,text)
        write (text,'(a9,i3)') 'stress k=',k
        call pipe_compare(stress,iu,text)
        write (text,'(a9,i3)') 'u.n    k=',k
        call pipe_compare(u(1-nbdy,1-nbdy,k,n),iu,text)
      endif
c
cdiag if (k.eq.1) then
cdiag write (lp,100) nstep
cdiag do j=jtest-1,jtest+1
cdiag do i=itest-1,itest+1
cdiag write (lp,'(2i5,i3,8f8.3)') i,j,k,
cdiag.  util4(i,j),u(i,j,k,n),-delt1*gradx(i,j)*scuxi(i,j),
cdiag.  -delt1*
cdiag. .25*(utotm(i+1,j)**2+vtotm(i  ,j)**2+vtotm(i  ,j+1)**2
cdiag.     -utotm(i-1,j)**2-vtotm(i-1,j)**2-vtotm(i-1,j+1)**2)*
cdiag.      scuxi(i,j),
cdiag.   delt1*(vflux(i  ,j)+vflux(i  ,j+1)+vflux(i-1,j)+vflux(i-1,j+1))
cdiag.        *(potvor(i,j)+potvor(i,j+1))*.125,
cdiag.  -delt1*ubrhs(i,j),delt1*stress(i,j),
cdiag.  -delt1*(uflux1(i,j)-uflux1(i-1,j)
cdiag.         +uflux3(i,j)-uflux2(i,j))/
cdiag.         (scu2(i,j)*max(dpu(i,j,k,m),onemm))
cdiag enddo
cdiag enddo
cdiag endif
 100    format(i9,8x,'uold    unew   gradp  nonlin   corio',
     &3x,'ubrhs  stress    fric')
c
c --- ----------
c --- v equation
c --- ----------
c
c --- deformation-dependent eddy viscosity coefficient
c
*         wtime2(15,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            visv(i,j)=max(veldff,viscos*
     &      sqrt(.5*(defor1(i,j)+defor1(i,j-1)+
     &               defor2(i,j)+defor2(i+1,j))))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
*         wtime2(16,k) = wtime()
      do i=1-margin,ii+margin
        do l=1,jsv(i)
          visv(i,jfv(i,l)-1)=visv(i,jfv(i,l))
          visv(i,jlv(i,l)+1)=visv(i,jlv(i,l))
        enddo
      enddo
c
!$OMP PARALLEL DO PRIVATE(j,l,i,ia,ib,
!$OMP&                    dpxy,dpia,dpib,visca,viscb,scvxa,scvxb)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c
c ---   longitudinal turb. momentum flux (at mass points)
c
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            if     (iv(i,j)+iv(i,j+1).gt.0) then
            vflux1(i,j)=(visv(i,j)+visv(i,j+1))*(dl2v(i,j)-dl2v(i,j+1))
     &                   *hfharm(max(dpv(i,j  ,k,m),onemm),
     &                           max(dpv(i,j+1,k,m),onemm))
     &                   *scp2(i,j)*2./(scvy(i,j)+scvy(i,j+1))
            endif
          enddo
        enddo
c
c ---   lateral turb. momentum flux (at vorticity points)
c ---   (left and right fluxes are evaluated separately because of sidewalls)
c
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ia=i-1
            ib=i+1
            dpxy=max(dpv(i ,j,k,m),onemm)
            dpia=max(dpv(ia,j,k,m),onemm)
            dpib=max(dpv(ib,j,k,m),onemm)
c
c --- check whether variables along coast have been initialized correctly
cdiag       if (k.eq.kk) then
cdiag         if (iv(ia,j).eq.0 .and. dpv(ia,j,k,m).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.            '  error - nonzero dpv(ia):',dpv(ia,j,k,m)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag         if (iv(ib,j).eq.0 .and. dpv(ib,j,k,m).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.            '  error - nonzero dpv(ib):',dpv(ib,j,k,m)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag       endif
            if (iv(ia,j).eq.0) then
              visca=visv(i ,j)
              scvxa=scvx(i, j)
            else
              visca=visv(ia,j)
              scvxa=scvx(ia,j)
            endif
            if (iv(ib,j).eq.0) then
              viscb=visv(i ,j)
              scvxb=scvx(i, j)
            else
              viscb=visv(ib,j)
              scvxb=scvx(ib,j)
            endif
            vflux2(i,j)=(visv(i,j)+visca)*(dl2via(i,j)-dl2v(i,j))
     &                  *hfharm(dpia+wgtia(i,j)*(dpxy-dpia),dpxy)
     &                  *scq2(i ,j)*2./(scvx(i,j)+scvxa)
            vflux3(i,j)=(visv(i,j)+viscb)*(dl2v(i,j)-dl2vib(i,j))
     &                  *hfharm(dpib+wgtib(i,j)*(dpxy-dpib),dpxy)
     &                  *scq2(ib,j)*2./(scvx(i,j)+scvxb)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
c --- pressure force in y direction
c --- ('scheme 2' from appendix -a- in bleck-smith paper)
c
*         wtime2(17,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            util2(i,j)=max(0.,min(depthv(i,j)-pv(i,j,k),h1))
            pgfy(i,j)=util2(i,j)*
     &          (montg(i,j,k)-montg(i,j-1,k)+
     &            (thstar(i,j,k)-thstar(i,j-1,k))
     &            *(p(i,j,k+1)*p(i,j-1,k+1)
     &             -p(i,j,k  )*p(i,j-1,k  ))*thref**2
     &            /(dp(i,j,k,m)+dp(i,j-1,k,m)+epsil))
          enddo
        enddo
      enddo
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        write (text,'(a9,i3)') 'pv     k=',k
        call pipe_compare(pv(1-nbdy,1-nbdy,k),iv,text)
        write (text,'(a9,i3)') 'depthv k=',k
        call pipe_compare(depthv,iv,text)
        write (text,'(a9,i3)') 'util2  k=',k
        call pipe_compare(util2, iv,text)
        util4(1:ii,1:jj) = montg( 1:ii,1:jj,k)-montg( 1:ii,0:jj-1,k)
        write (text,'(a9,i3)') 'montgY k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = thstar(1:ii,1:jj,k)-thstar(1:ii,0:jj-1,k)
        write (text,'(a9,i3)') 'thstaY k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k)  *p(1:ii,0:jj-1,k)
        write (text,'(a9,i3)') 'pSQ    k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = p(1:ii,1:jj,k+1)*p(1:ii,0:jj-1,k+1)
        write (text,'(a9,i3)') 'pSQ    k=',k+1
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) = dp(1:ii,1:jj,k,m)+dp(1:ii,0:jj-1,k,m)+epsil
        write (text,'(a9,i3)') 'dp.m+  k=',k
        call pipe_compare(util4, iv,text)
        write (text,'(a9,i3)') 'pgfy   k=',k
        call pipe_compare(pgfy,  iv,text)
      endif
c
*         wtime2(18,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i,ia,ib)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ia=i-1
            ib=i+1
c
c --- check whether variables along coast have been initialized correctly
cdiag       if (k.eq.kk) then
cdiag         if (iv(ia,j).eq.0 .and. pgfy(ia,j).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.            '  error - nonzero pgfy(ia):',pgfy(ia,j)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag         if (iv(ib,j).eq.0 .and. pgfy(ib,j).ne.0.) then
cdiag!$OMP CRITICAL
cdiag           write(lp,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.            '  error - nonzero pgfy(ib):',pgfy(ib,j)
cdiag!$OMP END CRITICAL
cdiag         endif
cdiag       endif
c
            grady(i,j)=(pgfy(i,j)+(h1-util2(i,j))*
     &        (pgfy (ia ,j)+pgfy (ib ,j)
     &        +pgfy (i,j-1)+pgfy (i,j+1))/
     &        (util2(ia ,j)+util2(ib ,j)
     &        +util2(i,j-1)+util2(i,j+1)+epsil))/h1
          enddo
        enddo
      enddo
c
*         wtime2(19,k) = wtime()
!$OMP PARALLEL DO PRIVATE(j,l,i,ptopl,pbotl,pstres)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
c
            ptopl=min(depthv(i,j),.5*(p(i,j,k  )+p(i,j-1,k  )))
            pbotl=min(depthv(i,j),.5*(p(i,j,k+1)+p(i,j-1,k+1)))
            if(hybrid .and. mxlkrt) then
              pstres=.5*(dpmixl(i,j,m)+dpmixl(i,j-1,m))
            else
              pstres=dpv(i,j,1,m)
            endif
c
c ---       top and bottom boundary layer stress
            stress(i,j)=(-vtotn(i,j)*.5*(drag(i,j)+drag(i,j-1))*
     &           (max(depthv(i,j)-thkbop,          pbotl       )
     &           -max(depthv(i,j)-thkbop,min(ptopl,pbotl-onemm)))
     &           +stresy(i,j)*thref*(min(pstres,pbotl+onemm)
     &                              -min(pstres,ptopl      )))
     &           /max(dpv(i,j,k,m),onemm)
c
c ---       time smoothing of -v- field  (part 1)
            v(i,j,k,m)=v(i,j,k,m)*(wuv1*dpv(i,j,k,m)+onemm)
     &                +v(i,j,k,n)* wuv2*dpv(i,j,k,n)
c
cdiag       util4(i,j)=v(i,j,k,n)
            v(i,j,k,n)=v(i,j,k,n)+delt1*(-scvyi(i,j)*(grady(i,j)
     &      +.25*(vtotm(i,  j+1)**2+utotm(i,  j  )**2
     &           +utotm(i+1,j  )**2-vtotm(i,  j-1)**2
     &           -utotm(i,  j-1)**2-utotm(i+1,j-1)**2))
     &      -.125*(uflux(i,  j  )+uflux(i+1,j  )
     &            +uflux(i,  j-1)+uflux(i+1,j-1))
     &           *(potvor(i,j)+potvor(i+1,j))
     &      - vbrhs(i,j) + stress(i,j)
     &      -(vflux1(i,j)-vflux1(i,j-1)
     &       +vflux3(i,j)-vflux2(i,j  ))
     &      /(scv2(i,j)*max(dpv(i,j,k,m),onemm)))
          enddo
        enddo
      enddo
c
      if     (lpipe .and. lpipe_momtum .and. k.eq.1) then
c ---   compare two model runs.
        util4(1:ii,1:jj) =  vtotm(1:ii,  2:jj+1)**2
     &                     +utotm(1:ii,  1:jj  )**2
     &                     +utotm(2:ii+1,1:jj  )**2
     &                     -vtotm(1:ii,  0:jj-1)**2
     &                     -utotm(1:ii,  0:jj-1)**2
     &                     -utotm(2:ii+1,0:jj-1)**2
        write (text,'(a9,i3)') 'totm   k=',k
        call pipe_compare(util4, iv,text)
        write (text,'(a9,i3)') 'grady  k=',k
        call pipe_compare(grady, iv,text)
        write (text,'(a9,i3)') 'vbrhs  k=',k
        call pipe_compare(vbrhs, iv,text)
        write (text,'(a9,i3)') 'stress k=',k
        call pipe_compare(stress,iv,text)
        util4(1:ii,1:jj) =  uflux(1:ii,  1:jj  )
     &                     +uflux(2:ii+1,1:jj  )
     &                     +uflux(1:ii,  0:jj-1)
     &                     +uflux(2:ii+1,0:jj-1)
        write (text,'(a9,i3)') 'uflux  k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) =  potvor(1:ii,  1:jj  )
     &                     +potvor(2:ii+1,1:jj  )
        write (text,'(a9,i3)') 'potvor k=',k
        call pipe_compare(util4, iv,text)
        util4(1:ii,1:jj) =  vflux1(1:ii,  1:jj  )
     &                     -vflux1(1:ii,  0:jj-1)
     &                     +vflux3(1:ii,  1:jj  )
     &                     -vflux2(1:ii,  1:jj  )
        write (text,'(a9,i3)') 'vflux  k=',k
        call pipe_compare(util4, iv,text)
        write (text,'(a9,i3)') 'v.n    k=',k
        call pipe_compare(v(1-nbdy,1-nbdy,k,n),iv,text)
      endif
c
cdiag if (k.eq.1) then
cdiag write (lp,101) nstep
cdiag do j=jtest-1,jtest+1
cdiag do i=itest-1,itest+1
cdiag write (lp,'(2i5,i3,8f8.3)') i,j,k,
cdiag.  util4(i,j),v(i,j,k,n),-delt1*grady(i,j)*scvyi(i,j),
cdiag.  -delt1*
cdiag. .25*(vtotm(i,j+1)**2+utotm(i,j  )**2+utotm(i+1,j  )**2
cdiag.     -vtotm(i,j-1)**2-utotm(i,j-1)**2-utotm(i+1,j-1)**2)*
cdiag.      scvyi(i,j),
cdiag.  -delt1*(uflux(i,j  )+uflux(i+1,j  )+uflux(i,j-1)+uflux(i+1,j-1))
cdiag.        *(potvor(i,j)+potvor(i+1,j))*.125,
cdiag.  -delt1*vbrhs(i,j),delt1*stress(i,j),
cdiag.  -delt1*(vflux1(i,j)-vflux1(i,j-1)
cdiag.         +vflux3(i,j)-vflux2(i,j))/
cdiag.         (scv2(i,j)*max(dpv(i,j,k,m),onemm))
cdiag enddo
cdiag enddo
cdiag endif
 101    format(i9,8x,'vold    vnew   gradp  nonlin   corio',
     &3x,'vbrhs  stress    fric')
c
*         wtime2(20,k) = wtime()
 9    continue  ! k=1,kk
c
      dt1inv = 1./delt1
c
*        wtime1( 6) = wtime()
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 14 j=1-margin,jj+margin
      do 14 k=1,kk
c
      do 12 l=1,isp(j)
      do 12 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
 12   continue
c
c --- compute new -dpu,dpv- field. save old -dpu,dpv- values in -pu,pv-.
c
      do 13 l=1,isu(j)
      do 13 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
      pu(i,j,k+1)=dpu(i,j,k,n)
 13   continue
c
      do 14 l=1,isv(j)
      do 14 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
      pv(i,j,k+1)=dpv(i,j,k,n)
 14   continue
!$OMP END PARALLEL DO
c
*        wtime1( 7) = wtime()
      call dpudpv(p,depthu,depthv,dpu(1-nbdy,1-nbdy,1,n),
     &                            dpv(1-nbdy,1-nbdy,1,n))
*        wtime1( 8) = wtime()
c
c --- substitute depth-weighted averages for (u,v) at massless grid points.
c --- (scan layers in top-down direction to save time.)
c --- extract barotropic velocities generated during most recent baroclinic
c --- time step and use them to force barotropic flow field.
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k,ka,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 30 j=1-margin,jj+margin
c
      do 31 l=1,isu(j)
c
      do 32 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
      utotn(i,j)=0.
 32   continue
      do 33 k=1,kk
      ka=max(1,k-1)
      do 33 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
      q=min(dpu(i,j,k,m),dpu(i,j,k,n),cutoff)
      u(i,j,k,n)=(u(i,j,k,n)*q+u(i,j,ka,n)*(cutoff-q))/cutoff
      utotn(i,j)=utotn(i,j)+u(i,j,k,n)*dpu(i,j,k,n)
 33   continue
      do 31 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
      utotn(i,j)=utotn(i,j)/depthu(i,j)
 31   continue
c
      do 30 l=1,isv(j)
c
      do 34 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
      vtotn(i,j)=0.
 34   continue
      do 35 k=1,kk
      ka=max(1,k-1)
      do 35 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
      q=min(dpv(i,j,k,m),dpv(i,j,k,n),cutoff)
      v(i,j,k,n)=(v(i,j,k,n)*q+v(i,j,ka,n)*(cutoff-q))/cutoff
      vtotn(i,j)=vtotn(i,j)+v(i,j,k,n)*dpv(i,j,k,n)
 35   continue
      do 30 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
      vtotn(i,j)=vtotn(i,j)/depthv(i,j)
 30   continue
!$OMP END PARALLEL DO
c
c --- time smoothing of -u,v- fields  (part 2)
c
*        wtime1( 9) = wtime()
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 22 j=1-margin,jj+margin
      do 22 k=1,kk
c
      do 24 l=1,isu(j)
      do 24 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
      u(i,j,k,n)=u(i,j,k,n)-utotn(i,j)
      u(i,j,k,m)=(u(i,j,k,m)+u(i,j,k,n)*wuv2*dpu(i,j,k,n))/
     &   (wuv1*dpu(i,j,k,m)+onemm+wuv2*(pu(i,j,k+1)+dpu(i,j,k,n)))
 24   continue
c
      do 22 l=1,isv(j)
      do 22 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
      v(i,j,k,n)=v(i,j,k,n)-vtotn(i,j)
      v(i,j,k,m)=(v(i,j,k,m)+v(i,j,k,n)*wuv2*dpv(i,j,k,n))/
     &   (wuv1*dpv(i,j,k,m)+onemm+wuv2*(pv(i,j,k+1)+dpv(i,j,k,n)))
 22   continue
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 867 j=1-margin,jj+margin
c
      do 865 l=1,isu(j)
      do 865 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
      utotn(i,j)=utotn(i,j)*dt1inv
      ubavg(i,j,n)=ubavg(i,j,m)
 865  continue
c
      do 866 l=1,isv(j)
      do 866 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
      vtotn(i,j)=vtotn(i,j)*dt1inv
      vbavg(i,j,n)=vbavg(i,j,m)
 866  continue
c
      do 867 l=1,isp(j)
      do 867 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      pbavg(i,j,n)=pbavg(i,j,m)
 867  continue
!$OMP END PARALLEL DO
c
*        wtime1(10) = wtime()
*        if     (mod(nstep,100).eq.0) then
* ---      timer printout.
*          do i= 2,10
*            write(lp,'(4x,a,i3,a,f9.5,a)')
*    .        'momtum   section',i,' time =',
*    .        wtime1(i) - wtime1(i-1),' seconds'
*          enddo
*          do i= 2,20
*            wtimes=0.0
*            do k=1,kk
*              wtimes = wtimes + (wtime2(i,k)-wtime2(i-1,k))
*            enddo
*            write(lp,'(4x,a,i3,a,f9.5,a)')
*    .        'momtum k section',i,' time =',
*    .        wtimes,' seconds'
*          enddo
*          call flush(lp)
*        endif
c
      return
      end
c
c
c> Revision history
c>
c> Mar. 1995 - changed min.depth in pot.vort. calculation from 1 mm to 1 m
c>             (loops 812,802,803)
c> July 1997 - transferred -visc- and -vort- from common to local
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c> Aug. 1997 - added some loops that used to be in barotp.f
c> Aug. 1997 - transferred -wgtia,wgtib,wgtja,wgtjb- from common to local
c> Oct. 1999 - ekman depth calculation added in loops 66-70 - based on ustar
c>             which is calculated from the stress
c> Oct. 1999 - surface stress terms handled differently depending on m.l.
c>             model
c> Dec. 1999 - switched to biharmonic friction
c> Jan. 2000 - added biharm to select between laplacian and biharmonic
c> Jan. 2000 - ekman depth calculation removed - surface momentum flux
c>             distributed within layer one for all mixed layer models
c> Jan. 2000 - added r. bleck's re-formulation of -pgfx- and -pgfy-
c>             to properly commute pressure torque from layer to layer
c> Jul. 2000 - loop reordering and logic changes for OpenMP
c> Aug. 2000 - loop 117 executed only when hybrid vertical coordinate is used
