      subroutine profile_lat(theta,press,xlat)
      implicit none
c
      real theta,press,xlat
c
      integer        lp
      common/linepr/ lp
      save  /linepr/
c
c --- this routine returns either:
c
c ---    pressure as function of density  and latitude
c ---    density  as function of pressure and latitude
c
c --- set press < 0.0 on input to return density
c
c --- typically invoked via either poflat or roflat.
c
      integer ix,kz
      real    p1,p2,pinthi,pintlo,pz,thet,thetlo,thethi,x,xla,z
c
      integer    kdpth,klat
      parameter (kdpth=14,klat=27)  ! kdpth>1, klat>3
c
      real onem,thet1,thet2,dthet,xlat1,xlat2,dlat
      real pdat(kdpth,klat)
c
      data onem/9806./  ! SI units
      data thet1,thet2,dthet/21.5,28.0,0.5/
      data xlat1,xlat2,dlat/5.,31.,27./
c
c---  depth (m) of isopycnals of potential density 22.0, ... , 28.0
c---  at latitudes  5n ... 31n  for IASd (source: levitus atlas)
c
      data pdat/
     + 0.,9.,19.,27.,49.,61., 82., 92.,114.,130.,164.,234.,756.,4000.,  ! 5n
     + 0.,9.,19.,27.,49.,61., 82., 92.,114.,130.,164.,234.,756.,4000.,  ! 6n
     + 0.,9.,19.,27.,49.,61., 82., 92.,114.,130.,164.,234.,756.,4000.,  ! 7n
     + 0.,9.,19.,27.,49.,61., 82., 92.,114.,130.,164.,234.,756.,4000.,  ! 8n
     + 0.,8.,15.,30.,47.,64., 80., 92.,111.,130.,171.,261.,811.,4000.,  ! 9n
     + 0.,6., 9.,23.,43.,58., 77., 89.,111.,132.,184.,290.,868.,4000.,  !10n
     + 0.,3., 7.,12.,35.,45., 65., 75., 98.,119.,184.,300.,822.,4000.,  !11n
     + 0.,1., 7.,13.,34.,49., 67., 81.,103.,127.,185.,317.,872.,4000.,  !12n
     + 0.,0., 9.,18.,40.,65., 79.,100.,119.,150.,201.,341.,909.,4000.,  !13n
     + 0.,0., 8.,21.,50.,79., 94.,117.,137.,171.,220.,369.,896.,4000.,  !14n
     + 0.,0., 4.,21.,52.,88.,103.,128.,152.,187.,239.,396.,882.,4000.,  !15n
     + 0.,0., 1.,18.,50.,92.,107.,131.,156.,190.,245.,407.,835.,4000.,  !16n
     + 0.,0., 0.,14.,50.,96.,111.,137.,162.,200.,259.,445.,856.,4000.,  !17n
     + 0.,0., 0., 9.,54.,86.,107.,136.,165.,205.,275.,486.,875.,4000.,  !18n
     + 0.,0., 0., 3.,53.,68.,100.,121.,160.,193.,278.,500.,849.,4000.,  !19n
     + 0.,0., 0., 0.,40.,62., 89.,113.,151.,186.,273.,511.,840.,4000.,  !20n
     + 0.,0., 0., 0.,25.,60., 77.,108.,141.,185.,266.,536.,830.,4000.,  !21n
     + 0.,0., 0., 0.,15.,50., 63., 98.,127.,177.,252.,531.,776.,4000.,  !22n
     + 0.,0., 0., 0.,14.,44., 62., 91.,128.,176.,260.,526.,777.,4000.,  !23n
     + 0.,0., 0., 1.,13.,34., 60., 78.,126.,169.,269.,536.,826.,4000.,  !24n
     + 0.,0., 0., 3.,10.,24., 53., 67.,113.,156.,264.,531.,798.,4000.,  !25n
     + 0.,0., 1., 4., 9.,22., 46., 63., 99.,148.,258.,543.,783.,4000.,  !26n
     + 0.,0., 3., 5., 9.,18., 39., 57., 87.,137.,255.,558.,785.,4000.,  !27n
     + 0.,0., 3., 5., 8.,16., 33., 51., 75.,122.,247.,559.,719.,4000.,  !28n
     + 0.,0., 3., 4., 6.,14., 28., 46., 65.,107.,243.,571.,680.,4000.,  !29n
     + 0.,0., 1., 1., 3.,10., 22., 46., 59.,106.,258.,684.,742.,4000.,  !30n
     + 0.,0., 0., 0., 1., 8., 17., 47., 57.,106.,264.,770.,804.,4000./  !31n
c
c---  quasi-hermite interpolation function (0 < xx < 1)
c
      real parabl,xx,a,b,c
      parabl(xx,a,b,c)=b+.5*xx*(c-a+xx*(a+c-b-b))
c
      xla=(xlat*57.29578-xlat1)/dlat+1.
      ix=max(2,min(klat-2,int(xla)))
      x=max(0.,min(1.,xla-float(ix)))
c
      if     (press.lt.0.0) then
c
c ----  pressure from density.
c
        thet=(theta-thet1)/dthet+1.
        if     (thet.lt.1.0) then
        press=0.0
        else  ! normal case
        kz=max(1,min(kdpth-1,int(thet)))
        z=max(0.,min(1.,thet-float(kz)))
c
c ---   horizontal/vertical interpolation: quasi-hermite/linear
c
        p1=parabl(   x,pdat(kz  ,ix-1),pdat(kz  ,ix  ),pdat(kz  ,ix+1))
        p2=parabl(1.-x,pdat(kz  ,ix+2),pdat(kz  ,ix+1),pdat(kz  ,ix  ))
        pintlo=p1*(1.-x)+p2*x
        p1=parabl(   x,pdat(kz+1,ix-1),pdat(kz+1,ix  ),pdat(kz+1,ix+1))
        p2=parabl(1.-x,pdat(kz+1,ix+2),pdat(kz+1,ix+1),pdat(kz+1,ix  ))
        pinthi=p1*(1.-x)+p2*x
        press =(pintlo*(1.-z)+pinthi*z)*onem
        endif
cdiag   write (lp,'('' poflat'',2f7.2,2i5,2f7.2,f7.1)')
cdiag&    theta,xlat*57.29578,ix,kz,x,z,press/onem
      else
c
c ----  density from pressure.
c
        pz=press/onem
        kz=1
        p1=parabl(   x,pdat(kz,ix-1),pdat(kz,ix  ),pdat(kz,ix+1))
        p2=parabl(1.-x,pdat(kz,ix+2),pdat(kz,ix+1),pdat(kz,ix  ))
        pinthi=p1*(1.-x)+p2*x
        if     (pinthi.ge.pz) then
        theta=thet1
        else  ! normal range
        do kz= 2,kdpth
          pintlo=pinthi
          p1=parabl(   x,pdat(kz,ix-1),pdat(kz,ix  ),pdat(kz,ix+1))
          p2=parabl(1.-x,pdat(kz,ix+2),pdat(kz,ix+1),pdat(kz,ix  ))
          pinthi=p1*(1.-x)+p2*x
          if     (pinthi.ge.pz) then
            exit
          elseif (kz.eq.kdpth) then
            exit
          endif
        enddo
        z=max((pinthi-pz)/(pinthi-pintlo),0.0)
        theta=thet1+(kz-z-1.0)*dthet
        endif
cdiag   write (lp,'('' roflat'',2f7.2,2i5,2f7.2,f7.1)')
cdiag&    theta,xlat*57.29578,ix,kz,x,z,pz
      endif
      return
      end

      real function poflat(theta,xlat)
      implicit none
c
      real theta,xlat
c
c --- returns pressure as function of density and latitude
c
      real press
      press = -1.0
      call profile_lat(theta,press,xlat)
      poflat = press
      return
      end

      real function roflat(press,xlat)
      implicit none
c
      real press,xlat
c
c --- returns density as function of pressure and latitude
c
      real theta
c
      call profile_lat(theta,press,xlat)
      roflat = theta
      return
      end

c
c> Revision history
c>
c> May  2000 - conversion to SI units
c> Aug  2001 - added roflat and profile_lat to poflat.
