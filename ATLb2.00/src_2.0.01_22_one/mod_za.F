      module mod_za
      use mod_xc  ! HYCOM communication API
c
      implicit none
c
c --- HYCOM I/O interface.
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
