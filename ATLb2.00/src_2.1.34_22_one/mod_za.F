      module mod_za
      use mod_xc  ! HYCOM communication API
c
      implicit none
c
c --- HYCOM I/O interface.
c
c --- See README.src.mod_za for more details.
c
      private zaiordd,zaiowrd
c
      contains
#if ! defined(MPI) && ! defined(SHMEM)
# include "mod_za_sm.F"
#elif defined(SERIAL_IO)
# include "mod_za_mp1.F"
#else
# include "mod_za_mp.F"
#endif
      end module mod_za

#if defined(ENDIAN_IO)
      subroutine zaio_endian(a,n)
      implicit none
c
      integer,         intent(in)    :: n
      integer(kind=4), intent(inout) :: a(n)  ! 4-bytes
c
c**********
c*
c 1)  swap the endian-ness of the array.
c
c 2)  assumes integer(kind=1) and integer(kind=4) ocupy one and four
c     bytes respectively.
c*
c**********
c
      integer k
c
      integer(kind=4) ii4,   io4     ! 4-bytes
      common/czioxe/  ii4,   io4     ! helps prevent unwanted optimization
      save  /czioxe/
c
      integer(kind=1) ii1(4),io1(4)  ! 1-byte
      equivalence    (ii4,ii1(1)), (io4,io1(1))  ! non-standard f90
c
      do k= 1,n
        ii4 = a(k)
        io1(1) = ii1(4)
        io1(2) = ii1(3)
        io1(3) = ii1(2)
        io1(4) = ii1(1)
        a(k) = io4
      enddo
      return
      end subroutine zaio_endian
#endif /* ENDIAN_IO */
