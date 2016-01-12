      real function poflat(theta,xlat)
      implicit none
c
      real theta,xlat
c
c --- this routine returns pressure as function of density (sigma) and latitude
c
      integer        lp
      common/linepr/ lp
      save  /linepr/
c
      integer ix,kz
      real    p1,p2,pinthi,pintlo,thet,x,xla,z
c
      integer    kdpth,klat
      parameter (kdpth=15,klat=20)
      real pdat(kdpth,klat)
c
      real onem,thet1,thet2,dthet,xlat1,xlat2,dlat
c
      data onem/9806./  ! SI units
      data thet1,thet2,dthet/21.0,28.0,0.5/
      data xlat1,xlat2,dlat/-30.,65.,5./
c
c---  depth (m) of isopycnals of potential density 21.0, 21.5, ... , 28.0
c---  at latitudes  30s ... 65n  for ATLa (source: levitus atlas)
c
      data  pdat/
     + 0.,0.,0.,0., 0., 0., 0., 0., 7., 54.,152.,217.,443.,938.,6000.  ! 30s
     +,0.,0.,0.,0., 0., 0., 0., 0., 7., 54.,152.,217.,443.,938.,6000.  ! 25s
     +,0.,0.,0.,0., 0., 0., 0., 0., 7., 54.,152.,217.,443.,938.,6000.  ! 20s
     +,0.,0.,0.,0., 0., 0., 0., 2.,25., 76.,145.,193.,393.,947.,6000.  ! 15s
     +,0.,0.,0.,0., 0., 0., 6.,33.,74., 98.,122.,160.,339.,990.,6000.  ! 10s
     +,0.,0.,0.,0., 0., 0.,20.,43.,68., 77., 96.,139.,380.,933.,6000.  !  5s
     +,1.,1.,1.,1., 2.,10.,37.,51.,66., 74., 92.,130.,413.,901.,6000.  !  0 
     +,1.,2.,2.,4., 9.,39.,57.,65.,78., 87.,103.,145.,348.,895.,6000.  !  5n
     +,0.,0.,3.,5.,10.,25.,43.,61.,76., 89.,111.,147.,296.,909.,6000.  ! 10n
     +,0.,0.,1.,3., 7.,20.,39.,59.,90.,113.,144.,191.,366.,881.,6000.  ! 15n
     +,0.,0.,0.,0., 2.,19.,36.,48.,76.,111.,162.,228.,472.,820.,6000.  ! 20n
     +,0.,0.,0.,0., 2., 7.,19.,28.,45., 71.,132.,223.,536.,794.,6000.  ! 25n
     +,0.,0.,0.,1., 1., 2., 7.,16.,27., 46., 85.,203.,602.,828.,6000.  ! 30n
     +,0.,0.,0.,0., 0., 1., 2.,11.,21., 37., 66.,147.,529.,868.,6000.  ! 35n
     +,0.,0.,0.,0., 0., 0., 1., 7.,12., 24., 47., 89.,286.,812.,6000.  ! 40n
     +,0.,0.,0.,0., 0., 1., 2., 3., 6., 10., 23., 59.,111.,676.,6000.  ! 45n
     +,0.,0.,0.,0., 0., 0., 0., 0., 0.,  1.,  3., 22., 74.,434.,6000.  ! 50n
     +,0.,0.,0.,0., 0., 0., 0., 0., 0.,  0.,  1.,  4., 37.,203.,6000.  ! 55n
     +,0.,0.,0.,0., 0., 0., 0., 0., 0.,  0.,  0.,  2., 10.,133.,6000.  ! 60n
     +,0.,0.,0.,0., 0., 0., 0., 0., 0.,  0.,  0.,  2., 12.,108.,6000.  ! 65n
     +/
c
c---  quasi-hermite interpolation function (0 < xx < 1)
c
      real parabl,xx,a,b,c
ccc   real qsihmt,d
      parabl(xx,a,b,c)=b+.5*xx*(c-a+xx*(a+c-b-b))
ccc   qsihmt(xx,a,b,c,d)=(1.-xx)*parabl(xx,a,b,c)+xx*parabl(1.-xx,d,c,b)
c
      xla=(xlat*57.29578-xlat1)/dlat+1.
      ix=max(2,min(klat-2,int(xla)))
      x=max(0.,min(1.,xla-float(ix)))
      thet=(theta-thet1)/dthet+1.
      kz=max(1,min(kdpth-1,int(thet)))
      z=max(0.,min(1.,thet-float(kz)))
c
c --- horizontal interpolation: quasi-hermite. vertical interpol.: linear
c
      p1=parabl(   x,pdat(kz  ,ix-1),pdat(kz  ,ix  ),pdat(kz  ,ix+1))
      p2=parabl(1.-x,pdat(kz  ,ix+2),pdat(kz  ,ix+1),pdat(kz  ,ix  ))
      pintlo=p1*(1.-x)+p2*x
      p1=parabl(   x,pdat(kz+1,ix-1),pdat(kz+1,ix  ),pdat(kz+1,ix+1))
      p2=parabl(1.-x,pdat(kz+1,ix+2),pdat(kz+1,ix+1),pdat(kz+1,ix  ))
      pinthi=p1*(1.-x)+p2*x
      poflat=(pintlo*(1.-z)+pinthi*z)*onem
ccc   write (lp,'('' poflat'',2f7.2,2i5,2f7.2,f7.1)')
ccc  .  theta,xlat*57.29578,ix,kz,x,z,poflat/onem
      return
      end
c
c> Revision history
c>
c> May  2000 - conversion to SI units
