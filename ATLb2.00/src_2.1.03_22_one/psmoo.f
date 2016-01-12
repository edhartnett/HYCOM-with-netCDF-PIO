      subroutine psmooth(a,margin_smooth)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      integer margin_smooth
      real    a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c --- ragged boundary version of basic 9-point smoothing routine.
c --- this routine is set up to smooth data carried at -p- points.
c --- see also psmooth_max (below).
c
c --- util1 used as workspace, so array a can't be util1
c
      integer i,ismth,j,jsmth,msmth
      real    qc,sh
c
      real    c(-1:1,-1:1)
      save    c
      data    c / 1.0, 2.0, 1.0,
     &            2.0, 4.0, 2.0,
     &            1.0, 2.0, 1.0 /
c
      qc = 1.0/sum(c(:,:))
c
      msmth = min(margin_smooth,nbdy-1)
c
      if     (margin.lt.msmth+1) then
c ---   update the halo
        call xctilr(a,1,1, msmth+1,msmth+1, halo_ps)
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC)
      do j=0-msmth,jj+msmth+1
        do i=0-msmth,ii+msmth+1
          util1(i,j) = a(i,j)
        enddo
      enddo
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,i,sh,jsmth,ismth)
!$OMP&         SCHEDULE(STATIC)
      do j=1-msmth,jj+msmth
        do i=1-msmth,ii+msmth
          if     (ip(i,j).eq.1) then
            sh = 0.0
            do jsmth= -1,1
              do ismth= -1,1
                if     (ip(i+ismth,j+jsmth).eq.1) then
                  sh = sh + c(ismth,jsmth)*util1(i+ismth,j+jsmth)
                else
                  sh = sh + c(ismth,jsmth)*util1(i,      j)
                endif
              enddo
            enddo
            a(i,j) = sh*qc
          endif
        enddo
      enddo
!$OMP END PARALLEL DO
      return
      end subroutine psmooth

      subroutine psmooth_max(a,margin_smooth)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      integer margin_smooth
      real    a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c --- ragged boundary version of basic 9-point smoothing routine.
c --- this routine is set up to smooth data carried at -p- points
c --- and to return the maximum of the original and smoothed value.
c --- see also psmooth (above).
c
c --- util1 used as workspace, so array a can't be util1
c
      integer i,ismth,j,jsmth,msmth
      real    qc,sh
c
      real    c(-1:1,-1:1)
      save    c
      data    c / 1.0, 2.0, 1.0,
     &            2.0, 4.0, 2.0,
     &            1.0, 2.0, 1.0 /
c
      qc = 1.0/sum(c(:,:))
c
      msmth = min(margin_smooth,nbdy-1)
c
      if     (margin.lt.msmth+1) then
c ---   update the halo
        call xctilr(a,1,1, msmth+1,msmth+1, halo_ps)
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC)
      do j=0-msmth,jj+msmth+1
        do i=0-msmth,ii+msmth+1
          util1(i,j) = a(i,j)
        enddo
      enddo
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,i,sh,jsmth,ismth)
!$OMP&         SCHEDULE(STATIC)
      do j=1-msmth,jj+msmth
        do i=1-msmth,ii+msmth
          if     (ip(i,j).eq.1) then
            sh = 0.0
            do jsmth= -1,1
              do ismth= -1,1
                if     (ip(i+ismth,j+jsmth).eq.1) then
                  sh = sh + c(ismth,jsmth)*util1(i+ismth,j+jsmth)
                else
                  sh = sh + c(ismth,jsmth)*util1(i,      j)
                endif
              enddo
            enddo
            a(i,j) = max( a(i,j), sh*qc )
          endif
        enddo
      enddo
!$OMP END PARALLEL DO
      return
      end subroutine psmooth_max
