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
      parameter (kdpth=13,klat=5)  ! kdpth>1, klat>3
c
      real onem,thet1,thet2,dthet,xlat1,xlat2,dlat
      real pdat(kdpth,klat)
c
      data onem/9806./  ! SI units
      data thet1,thet2,dthet/22.0,28.0,0.5/
      data xlat1,xlat2,dlat/-90.0,90.0,45.0/
c
c---  depth (m) of isopycnals of potential density 22.0, 22.5, ... , 28.0
c---  at latitudes  90s ... 90n  SYMa (no latitudinal variation).
c
      data pdat /
     + 0., 0., 0., 0., 0., 0.,0., 30.,150.,250.,500.,1250.,3000.  !90s
     +,0., 0., 0., 0., 0., 0.,0., 30.,150.,250.,500.,1250.,3000.  !45s
     +,0., 0., 0., 0., 0., 0.,0., 30.,150.,250.,500.,1250.,3000.  ! 0
     +,0., 0., 0., 0., 0., 0.,0., 30.,150.,250.,500.,1250.,3000.  !45n
     +,0., 0., 0., 0., 0., 0.,0., 30.,150.,250.,500.,1250.,3000.  !90n
     &/
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
