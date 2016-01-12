      subroutine stencl(k1,n)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer k1,n
c
c --- write 5 x 5 point cluster of grid point values centered on (itest,jtest)
c
      integer i,j,k
c
      call xcsync(flush_lp)
c
      if     (itest.gt.0 .and. itest.le.ii .and.
     &        jtest.gt.0 .and. jtest.le.jj      ) then
        do k=max(1,k1-1),min(kk,k1)		!  print out 2 adjacent layers
          write (lp,100) 'u at k=',k,'v at k=',k,
     &                   i0+itest-2,i0+itest+2,j0+jtest+2,j0+jtest-2
          write (lp,110)
     &      ((u(i,j,k,n),i=itest-2,itest+2),
     &       (v(i,j,k,n),i=itest-2,itest+2),
     &       j=jtest+2,jtest-2,-1)
          write (lp,100) 'uflx k=',k,'vflx k=',k,
     &                   i0+itest-2,i0+itest+2,j0+jtest+2,j0+jtest-2
          write (lp,120)
     &      ((uflux(i,j),i=itest-2,itest+2),
     &       (vflux(i,j),i=itest-2,itest+2),
     &       j=jtest+2,jtest-2,-1)
          write (lp,100) 'temp k=',k,'saln k=',k,
     &                   i0+itest-2,i0+itest+2,j0+jtest+2,j0+jtest-2
          write (lp,130)
     &      ((temp(i,j,k,n),i=itest-2,itest+2),
     &       (saln(i,j,k,n),i=itest-2,itest+2),
     &       j=jtest+2,jtest-2,-1)
          write (lp,100) 'dp_o k=',k,'dp_n k=',k,
     &                   i0+itest-2,i0+itest+2,j0+jtest+2,j0+jtest-2
          write (lp,140)
     &      ((dpold(i,j,k)*qonem,i=itest-2,itest+2),
     &       (dp(i,j,k,n)*qonem,i=itest-2,itest+2),
     &       j=jtest+2,jtest-2,-1)
          write (lp,100) 'pres k=',k+1,'th3d k=',k,
     &                   i0+itest-2,i0+itest+2,j0+jtest+2,j0+jtest-2
          write (lp,150)
     &      ((p(i,j,k+1)*qonem, i=itest-2,itest+2),
     &       (th3d(i,j,k,n)+thbase,i=itest-2,itest+2),
     &       j=jtest+2,jtest-2,-1)
        end do
      endif
c
      call xcsync(flush_lp)
      return
c
 100  format( 1x,a7,i3,
     &       30x,a7,i3,
     &        4x,'(',i4,':',i4,',',i4,':',i4,':-1)')
 110  format(   5f7.3,5x,5f7.3)    ! vel in m/s
 120  format(1p,5e7.0,5x,5e7.0)    ! flx
 130  format(   5f7.2,5x,5f7.2)    ! T&S
 140  format(   5f7.1,5x,5f7.1)    ! dp  in m
 150  format(   5f7.1,5x,5f7.2)    ! p   in m, th in sigmaT
      end
c
c
c> Revision history:
c>
c> Nov. 2000 - rotated cluster so that i is across and j up the page
