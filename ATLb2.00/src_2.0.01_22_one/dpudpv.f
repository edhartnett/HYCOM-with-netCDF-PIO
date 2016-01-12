      subroutine dpudpv(p,depthu,depthv,dpu,dpv)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
     & p
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     & dpu,dpv
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & depthu,depthv
c
c --- ----------------------------------
c --- define layer depth at  u,v  points
c --- ----------------------------------
c
      integer j
c
c --- using single row routine fixes SGI OpenMP bug.
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call dpudpvj(p,depthu,depthv,dpu,dpv, j)
      enddo
      return
      end
      subroutine dpudpvj(p,depthu,depthv,dpu,dpv, j)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c
      integer j
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
     & p
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     & dpu,dpv
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & depthu,depthv
c
c --- -----------------------------------------------
c --- define layer depth at  u,v  points,  single row
c --- -----------------------------------------------
c
      integer i,k,l
c
      do k=1,kk
c
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            dpu(i,j,k)=max(0.,
     .           min(depthu(i,j),.5*(p(i,j,k+1)+p(i-1,j,k+1)))-
     .           min(depthu(i,j),.5*(p(i,j,k  )+p(i-1,j,k  ))))
          enddo
        enddo
c
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            dpv(i,j,k)=max(0.,
     .           min(depthv(i,j),.5*(p(i,j,k+1)+p(i,j-1,k+1)))-
     .           min(depthv(i,j),.5*(p(i,j,k  )+p(i,j-1,k  ))))
          enddo
        enddo
c
      enddo
      return
      end
