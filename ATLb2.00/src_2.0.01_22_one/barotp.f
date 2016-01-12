      subroutine barotp(m,n)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
c
c --- micom version 2.8
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- ------------------------------------------------------------------------
c --- advance barotropic equations from baroclinic time level -m- to level -n-
c --- ------------------------------------------------------------------------
c
      logical    lpipe_barotp
      parameter (lpipe_barotp=.false.)
c
      real    q,utndcy,vtndcy
      real*8  sump
      integer i,j,l,lll,ml,nl,mn,ll,mbdy
      logical vthenu
c
      mbdy = 3
c
      call xctilr(utotn(  1-nbdy,1-nbdy    ),1, 1, 3,3, halo_uv)
      call xctilr(vtotn(  1-nbdy,1-nbdy    ),1, 1, 3,3, halo_vv)
c
      if     (lpipe .and. lpipe_barotp) then
c ---   compare two model runs.
        call pipe_compare(utotn,  iu,'barotp:utotn')
        call pipe_compare(vtotn,  iv,'barotp:vtotn')
        call pipe_compare(pvtrop, iq,'barotp:pvtrp')
      endif
c
      ml=n
      nl=3
c
c --- explicit time integration of barotropic flow (forward-backward scheme)
c --- in order to combine forward-backward scheme with leapfrog treatment of
c --- coriolis term, v-eqn must be solved before u-eqn every other time step
      vthenu=.false.
c
      do 840 lll=1,lstep
c
      call xctilr(pbavg(  1-nbdy,1-nbdy,1  ),1, 3, 3,3, halo_ps)
      call xctilr(ubavg(  1-nbdy,1-nbdy,1  ),1, 3, 3,3, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,1  ),1, 3, 3,3, halo_vv)
c
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare(pbavg(1-nbdy,1-nbdy,1), ip,'barot+:pbav1')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,2), ip,'barot+:pbav2')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,3), ip,'barot+:pbav3')
        call pipe_compare(ubavg(1-nbdy,1-nbdy,1), iu,'barot+:ubav1')
        call pipe_compare(ubavg(1-nbdy,1-nbdy,2), iu,'barot+:ubav2')
        call pipe_compare(ubavg(1-nbdy,1-nbdy,3), iu,'barot+:ubav3')
        call pipe_compare(vbavg(1-nbdy,1-nbdy,1), iv,'barot+:vbav1')
        call pipe_compare(vbavg(1-nbdy,1-nbdy,2), iv,'barot+:vbav2')
        call pipe_compare(vbavg(1-nbdy,1-nbdy,3), iv,'barot+:vbav3')
      endif
c
c --- continuity equation
c
c --- rhs: pbavg, ubavg+, vbavg+
c --- lhs: pbavg
c
      margin = mbdy - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isp(j)
          do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            pbavg(i,j,nl)=
     &        (1.-wbaro)*pbavg(i,j,ml)+
     &            wbaro *pbavg(i,j,nl)-
     &        (1.+wbaro)*dlt*
     &          (ubavg(i+1,j,ml)*depthu(i+1,j)*scuy(i+1,j)
     &          -ubavg(i  ,j,ml)*depthu(i  ,j)*scuy(i  ,j)
     &          +vbavg(i,j+1,ml)*depthv(i,j+1)*scvx(i,j+1)
     &          -vbavg(i,j  ,ml)*depthv(i,j  )*scvx(i,j  ))
     &        *scp2i(i,j)
        enddo
        enddo
      enddo
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare(pbavg(1-nbdy,1-nbdy,nl), ip,'barotp:pbavl')
      endif
c
      mn=ml
      if (vthenu) go to 901
c
c --- u momentum equation
c
 900  continue
c
c --- rhs: pbavg+, vbavg+, pvtrop+
c --- lhs: ubavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,utndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            utndcy=-thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)
     &      +(vbavg(i  ,j,  mn)*depthv(i  ,j)
     &       +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
     &       +vbavg(i-1,j,  mn)*depthv(i-1,j)
     &       +vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
     &       *(pvtrop(i,j)+pvtrop(i,j+1))*.125
c
            ubavg(i,j,nl)=
     &        (1.-wbaro)*ubavg(i,j,ml)+
     &            wbaro *ubavg(i,j,nl)+
     &        (1.+wbaro)*dlt*(utndcy+utotn(i,j))
c
cdiag       if (i.eq.itest.and.j.eq.jtest) then
cdiag         write (lp,'(i9,2i5,i3,3x,a,5f7.3)')
cdiag.          nstep,i+i0,j+j0,lll,
cdiag.          'u_old,u_new,p_grad,corio,u_star =',
cdiag.          ubavg(i,j,ml),ubavg(i,j,nl),
cdiag.          -thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)*dlt,
cdiag.          (vbavg(i  ,j,  mn)*depthv(i  ,j)
cdiag.          +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
cdiag.          +vbavg(i-1,j,  mn)*depthv(i-1,j)
cdiag.          +vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
cdiag.          *(pvtrop(i,j)+pvtrop(i,j+1))
cdiag.          *.125 * dlt,utotn(i,j) * dlt
cdiag       endif
*           if     (lpipe .and. lpipe_barotp) then
*             util1(i,j) = (pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)
*             util2(i,j) =  vbavg(i  ,j,  mn)*depthv(i  ,j)
*    &                     +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
*             util3(i,j) =  vbavg(i-1,j,  mn)*depthv(i-1,j)
*    &                     +vbavg(i-1,j+1,mn)*depthv(i-1,j+1)
*             util4(i,j) = (pvtrop(i,j)+pvtrop(i,j+1))
*           endif
          enddo
        enddo
      enddo
      if     (lpipe .and. lpipe_barotp) then
*       call pipe_compare(util1(1-nbdy,1-nbdy),    iu,'barotp:ubav1')
*       call pipe_compare(util2(1-nbdy,1-nbdy),    iu,'barotp:ubav2')
*       call pipe_compare(util3(1-nbdy,1-nbdy),    iu,'barotp:ubav3')
*       call pipe_compare(util4(1-nbdy,1-nbdy),    iu,'barotp:ubav4')
        call pipe_compare(ubavg(1-nbdy,1-nbdy,nl), iu,'barotp:ubavl')
      endif
c
      mn=nl
      if (vthenu) go to 902
c
c --- v momentum equation
c
 901  continue
c
c --- rhs: pbavg+, ubavg+, pvtrop+
c --- lhs: vbavg
c
      margin = margin - 1
c
!$OMP PARALLEL DO PRIVATE(j,l,i,vtndcy)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vtndcy=-thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)
     &      -(ubavg(i,  j  ,mn)*depthu(i,  j  )
     &       +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
     &       +ubavg(i,  j-1,mn)*depthu(i,  j-1)
     &       +ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
     &       *(pvtrop(i,j)+pvtrop(i+1,j))*.125
c
            vbavg(i,j,nl)=
     &        (1.-wbaro)*vbavg(i,j,ml)+
     &            wbaro *vbavg(i,j,nl)+
     &        (1.+wbaro)*dlt*(vtndcy+vtotn(i,j))
c
cdiag       if (i.eq.itest.and.j.eq.jtest) then
cdiag         write (lp,'(i9,2i5,i3,3x,a,5f7.3)')
cdiag.          nstep,i+i0,j+j0,lll,
cdiag.          'v_old,v_new,p_grad,corio,v_star =',
cdiag.          vbavg(i,j,ml),vbavg(i,j,nl),
cdiag.          -thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)*dlt,
cdiag.          -(ubavg(i,  j  ,mn)*depthu(i,j  )
cdiag.           +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
cdiag.           +ubavg(i,  j-1,mn)*depthu(i,j-1)
cdiag.           +ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
cdiag.          *(pvtrop(i,j)+pvtrop(i+1,j))
cdiag.          *.125 * dlt, vtotn(i,j) * dlt
cdiag       endif
*           if     (lpipe .and. lpipe_barotp) then
*             util1(i,j) = (pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)
*             util2(i,j) =  ubavg(i,  j  ,mn)*depthu(i,  j  )
*    &                     +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
*             util3(i,j) =  ubavg(i,  j-1,mn)*depthu(i,  j-1)
*    &                     +ubavg(i+1,j-1,mn)*depthu(i+1,j-1)
*             util4(i,j) = (pvtrop(i,j)+pvtrop(i+1,j))
*           endif
          enddo
        enddo
      enddo
      if     (lpipe .and. lpipe_barotp) then
*       call pipe_compare(util1(1-nbdy,1-nbdy),    iv,'barotp:vbav1')
*       call pipe_compare(util2(1-nbdy,1-nbdy),    iv,'barotp:vbav2')
*       call pipe_compare(util3(1-nbdy,1-nbdy),    iv,'barotp:vbav3')
*       call pipe_compare(util4(1-nbdy,1-nbdy),    iv,'barotp:vbav4')
        call pipe_compare(vbavg(1-nbdy,1-nbdy,nl), iv,'barotp:vbavl')
      endif
c
      mn=nl
      if (vthenu) go to 900
c
c --- switch order in which -u,v- equations are solved
 902  continue
      vthenu=.not.vthenu
c
      if     (lbflag.eq.1) then
        call latbdp(nl)
*     elseif (lbflag.eq.2) then
*       call latbdt(nl)
      endif
c
      ll=ml
      ml=nl
      nl=ll
c
 840  continue  ! lll=1,lstep
c
      if     (lbflag.eq.1) then
c
c ---   correct mean height.
c ---   this should not be required - so there may be a bug in the bc.
c
!$OMP   PARALLEL DO PRIVATE(j,l,i,sump)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1,jj
          sump = 0.d0
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              util1(i,j) = util1(i,j) + pbavg(i,j,nl)*scp2(i,j)
            enddo
          enddo
        enddo
        call xcsum(sump, util1,ip)
        q = sump/area
c
c ---   rhs: pbavg
c ---   lhs: pbavg
c
        margin = 0
c
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do l=1,isp(j)
            do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              pbavg(i,j,1) = pbavg(i,j,1) - q
              pbavg(i,j,2) = pbavg(i,j,2) - q
              pbavg(i,j,3) = pbavg(i,j,3) - q
            enddo
          enddo
        enddo
      endif
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare(pbavg(1-nbdy,1-nbdy,1), ip,'barotp:pbav1')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,2), ip,'barotp:pbav2')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,3), ip,'barotp:pbav3')
      endif
c
      return
      end subroutine barotp
c
c
c> Revision history:
c>
c> Mar. 1995 - changed vertical velocity averaging interval from 10 cm to 1 m
c>             (loops 33,35)
c> Mar. 1995 - changed order of loop nesting in loop 842
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c> Aug. 1997 - transferred loops preceding loop 840 to momeq2.f
c> Jan. 2000 - added latbdp for lateral boundary ports
