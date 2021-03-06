c-----------------------------------------------------------------------------
c --- START OF REGION AND TILING SPECIFIC PARAMETERS
c --- See: README.src.newregion for more details.
c
c --- itdm  = total grid dimension in i direction
c --- jtdm  = total grid dimension in j direction
c --- kdm   =       grid dimension in k direction
      integer    itdm,jtdm,kdm
      parameter (itdm=  93,jtdm=  65,kdm=22)  ! IASb0.50
c
c --- iqr   = maximum number of tiles in i direction
c --- jqr   = maximum number of tiles in j direction
      integer    iqr,jqr
      parameter (iqr= 1,jqr= 1)  ! single tile (TYPE=one or omp)
c
c --- idm   = maximum single tile grid dimension in i direction
c --- jdm   = maximum single tile grid dimension in j direction
      integer    idm,jdm
      parameter (idm=itdm,jdm=jtdm)  ! single tile (TYPE=one or omp)
c
c --- mxthrd= maximum number of OpenMP threads
      integer    mxthrd
      parameter (mxthrd=8)  ! NOMP=1,2,4,8
c
c --- kkwall= grid dimension in k direction for wall relax arrays
c --- kknest= grid dimension in k direction for nest relax arrays
      integer    kkwall,kknest
      parameter (kkwall=kdm)  ! must be 1 or kdm
      parameter (kknest=kdm)  ! must be 1 or kdm
c
c --- kkmy25= grid dimension in k direction for M-Y 2.5 arrays
      integer    kkmy25
      parameter (kkmy25=kdm)  ! must be -1 or kdm
c
c --- mxtrcr= maximum number of tracers
      integer    mxtrcr
      parameter (mxtrcr=1)
c
c ---   END OF REGION AND TILING SPECIFIC PARAMETERS
c-----------------------------------------------------------------------------
c
c --- halo size
      integer    nbdy
      parameter (nbdy=6)
c
c --- OpenMP will allocate jblk rows to each thread in turn
      integer    jblk
      parameter (jblk=(jdm+2*nbdy+mxthrd-1)/mxthrd)
c
c --- how far out the halo is valid (margin<=nbdy)
      integer      margin
      common/edge/ margin
      save  /edge/
c
c --- actual extent of this tile is (i0+1:i0+ii,j0+1:j0+jj,1:kk)
      integer      i0,j0,ii,jj
      common/dimi/ i0,j0,ii,jj
      save  /dimi/
      integer      kk
      parameter   (kk=kdm)
c
c --- ijqr  = maximum total number of active tiles (= ipr*jpr)
      integer    ijqr
      parameter (ijqr=iqr*jqr)
c
c --- ms-1  = max. number of interruptions of any tile row or column by land
      integer    ms
      parameter (ms=31)  ! should be enough for any region
c
c --- information in /gindex/ keeps do loops from running into land
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &         ip,iu,iv,iq, iuopn,ivopn
      integer, dimension (1-nbdy:jdm+nbdy,ms) :: 
     &         ifp,ilp,ifq,ilq,ifu,ilu,ifv,ilv
      integer, dimension (1-nbdy:idm+nbdy,ms) :: 
     &         jfp,jlp,jfq,jlq,jfu,jlu,jfv,jlv
      integer, dimension (1-nbdy:jdm+nbdy) :: 
     &         isp,isq,isu,isv
      integer, dimension (1-nbdy:idm+nbdy) :: 
     &         jsp,jsq,jsu,jsv
      common/gindex/ ip,iu,iv,iq, iuopn,ivopn,
     &               ifp,ilp,isp,jfp,jlp,jsp,ifq,ilq,isq,jfq,jlq,jsq,
     &               ifu,ilu,isu,jfu,jlu,jsu,ifv,ilv,isv,jfv,jlv,jsv
      save  /gindex/
c
c --- line printer unit (stdout)
      integer        lp
      common/linepr/ lp
      save  /linepr/
c-----------------------------------------------------------------------------
