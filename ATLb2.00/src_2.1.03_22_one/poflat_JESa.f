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
      parameter (kdpth=14,klat=22)  ! kdpth>1, klat>3
c
      real onem,thet1,thet2,dthet,xlat1,xlat2,dlat
      real pdat(kdpth,klat)
c
      data onem/9806./  ! SI units
      data thet1,thet2,dthet/21.0,27.5,0.5/
      data xlat1,xlat2,dlat/32.,53.,1./
c
c---  depth (m) of isopycnals of potential density 21.0, 21.5, ... , 27.5
c---  at latitudes  32N ... 53N  for JESa (source: NRL MODAS)
c
      data  pdat/
     + 0.,0.,0.,0.,0.,14.,24.,37.,54.,81.,99.,103.,105.,4000. ! 32N
     +,0.,0.,0.,0.,0.,14.,24.,37.,54.,81.,99.,103.,105.,4000. ! 33N
     +,0.,0.,0.,0.,0.,14.,24.,37.,54.,81.,99.,103.,105.,4000. ! 34N
     +,0.,0.,0.,0.,0.,12.,24.,35.,53.,78.,99.,106.,111.,4000. ! 35N
     +,0.,0.,0.,0.,0., 1.,19.,23.,42.,57.,92.,117.,159.,4000. ! 36N
     +,0.,0.,0.,0.,0., 0., 8.,21.,27.,44.,65.,108.,158.,4000. ! 37N
     +,0.,0.,0.,0.,0., 0., 1.,16.,20.,35.,48., 87.,143.,4000. ! 38N
     +,0.,0.,0.,0.,0., 0., 0., 8.,19.,24.,38., 59.,116.,4000. ! 39N
     +,0.,0.,0.,0.,0., 0., 0., 2.,13.,18.,27., 39., 75.,4000. ! 40N
     +,0.,0.,0.,0.,0., 0., 0., 0., 6.,15.,21., 33., 53.,4000. ! 41N
     +,0.,0.,0.,0.,0., 0., 0., 0., 3.,13.,19., 30., 47.,4000. ! 42N
     +,0.,0.,0.,0.,0., 0., 0., 0., 1.,13.,21., 32., 58.,4000. ! 43N
     +,0.,0.,0.,0.,0., 0., 0., 0., 0.,14.,23., 40., 73.,4000. ! 44N
     +,0.,0.,0.,0.,0., 0., 0., 0., 0.,11.,22., 36., 66.,4000. ! 45N
     +,0.,0.,0.,0.,0., 0., 0., 0., 0.,11.,21., 35., 67.,4000. ! 46N
     +,0.,0.,0.,0.,0., 0., 0., 0., 0., 9.,20., 30., 65.,4000. ! 47N
     +,0.,0.,0.,0.,0., 0., 0., 0., 0., 5.,20., 24., 64.,4000. ! 48N
     +,0.,0.,0.,0.,0., 0., 0., 0., 0., 5.,21., 24., 62.,4000. ! 49N
     +,0.,0.,0.,0.,0., 0., 0., 0., 0.,11.,19., 31., 58.,4000. ! 50N
     +,0.,0.,0.,0.,0., 0., 0., 0., 4.,16.,19., 35., 49.,4000. ! 51N
     +,0.,0.,0.,0.,0., 0., 0., 5.,15.,16.,24., 33., 33.,4000. ! 52N
     +,0.,0.,0.,0.,0., 0., 0., 5.,15.,16.,24., 33., 33.,4000. ! 53N
     +/
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
