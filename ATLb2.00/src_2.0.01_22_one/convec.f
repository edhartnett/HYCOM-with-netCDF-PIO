      subroutine convch(m,n)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
cdiag real sigup,uup,vup,siglo,ulo,vlo,coluin(idm),colout(idm)
      real q,tem,sal,thet,trc
      integer i,iter,j,k,ks,kp,l
c
      integer itmax
      data itmax/4/
c
      include 'stmt_fns.h'
c
c --- ---------------------
c --- convective adjustment
c --- ---------------------
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f8.2,f8.1))
cdiag if (itest.gt.0 .and. jtest.gt.0) write (lp,103) nstep,itest,jtest,
cdiag.  '  entering convec:  temp    saln    dens    thkns    dpth',
cdiag.  (k,temp(itest,jtest,k,n),saln(itest,jtest,k,n),
cdiag.  th3d(itest,jtest,k,n)+thbase,dp(itest,jtest,k,n)/onem,
cdiag.  p(itest,jtest,k+1)/onem,k=1,kk)
c
!$OMP PARALLEL DO PRIVATE(j,k,l,i,q)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 26 j=1-margin,jj+margin
c
c --- convection of u
c
      do 16 k=2,kk
c
      do 16 l=1,isu(j)
      do 16 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
      if (pu(i,j,k+1).lt.depthu(i,j)-onemm .and.
     .    th3d(i,j,k,  n)+th3d(i-1,j,k,  n).lt.
     .    th3d(i,j,k-1,n)+th3d(i-1,j,k-1,n)) then
cdiag   if (i.eq.itest .and. j.eq.jtest) then
cdiag   uup=u(i,j,k-1,n)
cdiag   ulo=u(i,j,k,  n)
cdiag   endif
        q=1.0-max(dpu(i,j,k,n),0.)/
     .   (max(dpu(i,j,k,n),0.)+max(dpu(i,j,k-1,n),onemm))
        u(i,j,k,  n)=u(i,j,k,n)+q*(u(i,j,k-1,n)-u(i,j,k,n))
        u(i,j,k-1,n)=u(i,j,k,n)
cdiag   if (i.eq.itest .and. j.eq.jtest) write (lp,100) nstep,i,j,k,
cdiag.  1,'  upr,lwr,final u:',uup,ulo,u(i,j,k,n),q
      end if
 16   continue
c
c --- convection of v
c
      do 26 k=2,kk
c
      do 26 l=1,isv(j)
      do 26 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
      pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k,n)
      if (pv(i,j,k+1).lt.depthv(i,j)-onemm .and.
     .    th3d(i,j,k,  n)+th3d(i,j-1,k,  n).lt.
     .    th3d(i,j,k-1,n)+th3d(i,j-1,k-1,n)) then
cdiag   vup=v(i,j,k-1,n)
cdiag   vlo=v(i,j,k,  n)
        q=1.0-max(dpv(i,j,k,n),0.)/
     .   (max(dpv(i,j,k,n),0.)+max(dpv(i,j,k-1,n),onemm))
        v(i,j,k,  n)=v(i,j,k,n)+q*(v(i,j,k-1,n)-v(i,j,k,n))
        v(i,j,k-1,n)=v(i,j,k,n)
cdiag   if (i.eq.itest .and. j.eq.jtest) write (lp,100) nstep,i,j,k,
cdiag.  1,'  upr,lwr,final v:',vup,vlo,v(i,j,k,n),q
      end if
 26   continue
!$OMP END PARALLEL DO
c
c --- convection of thermodynamical variables and tracer
c
!$OMP PARALLEL DO PRIVATE(j,l,k,i,iter,ks,kp,
!$OMP&                    q,tem,sal,thet,trc)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 1 j=1-margin,jj+margin
      do 1 l=1,isp(j)
c
ccc      do 9 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
ccc      coluin(i)=0.
ccc 9    colout(i)=0.
c
      do 12 k=1,kk
      do 12 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
ccc      coluin(i)=coluin(i)+temp(i,j,k,n)*dp(i,j,k,n)
 12   p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
c
      do 6 iter=1,itmax
c
      do 11 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      klist(i,j)=1
 11   util1(i,j)=dp(i,j,1,n)
c
      do 6 k=2,kk
c
      do 6 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      ks=klist(i,j)
      if (th3d(i,j,k,n).lt.th3d(i,j,ks,n)) then
        if (dp(i,j,k,n).lt.onemm) then
          th3d(i,j,k,n)=max(th3d(i,j,k,n),th3d(i,j,k-1,n))
          if(th3d(i,j,k-1,n).eq.th3d(i,j,k,n)) then
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=temp(i,j,k-1,n)
          end if
        else				!  dp > onemm
          if (iter.eq.itmax) then
            if (th3d(i,j,k,n).lt.th3d(i,j,ks,n)-sigjmp .and.
     .          dp(i,j,k,n).gt.onem) then
!$OMP         CRITICAL
              write (lp,'(i9,2i5,i3,a,i3,a,2f10.4)')
     .         nstep,i,j,k,' colmn unstbl after',iter-1,' iterations',
     .         th3d(i,j,ks,n)+thbase,th3d(i,j,k,n)+thbase
!$OMP         END CRITICAL
            endif
          else				!  it < itmax
cdiag       sigup=th3d(i,j,ks,n)
cdiag       siglo=th3d(i,j,k, n)
            util1(i,j)=util1(i,j)+dp(i,j,k,n)
            q=1.0-max(dp(i,j,k,n),.5*onemm)/max(util1(i,j),onemm)
            tem=temp(i,j,k,n)+q*(temp(i,j,ks,n)-temp(i,j,k,n))
            sal=saln(i,j,k,n)+q*(saln(i,j,ks,n)-saln(i,j,k,n))
            if(trcout)
     .       trc=tracer(i,j,k)+q*(tracer(i,j,ks)-tracer(i,j,k))
            thet=sig(tem,sal)-thbase
            do 10 kp=ks,k
            temp(i,j,kp,n)=tem
            saln(i,j,kp,n)=sal
            th3d(i,j,kp,n)=thet
            if(trcout) tracer(i,j,kp)=trc
 10         continue
c
cdiag if (i.eq.itest .and. j.eq.jtest) write (lp,100) nstep,i,j,k,iter,
cdiag.  '  upr,lwr,final dens:',(sigup+thbase),
cdiag.  (siglo+thbase),(th3d(i,j,k,n)+thbase),q
 100    format (i9,2i5,i3,'  it',i2,a,3f8.3,f5.2)
c
          end if
        end if
      else				!  th3d(kn) > th3d(ksn)
        klist(i,j)=k
        util1(i,j)=dp(i,j,k,n)
      end if
 6    continue
c
ccc      do 8 k=1,kk
ccc      do 8 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
ccc 8    colout(i)=colout(i)+temp(i,j,k,n)*dp(i,j,k,n)
c
ccc      do 1 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
ccc      if (abs((colout(i)-coluin(i))/coluin(i)).gt.1.e-6)
ccc     .  write (lp,'(i9,2i5,a/1p,3e14.6)') nstep,i,j,
ccc     .  '  column integral not conserved in convec:',
ccc     .  coluin(i),colout(i),(colout(i)-coluin(i))/coluin(i)
 1    continue
!$OMP END PARALLEL DO
c
cdiag if (itest.gt.0 .and. jtest.gt.0) write (lp,103) nstep,itest,jtest,
cdiag.  '  exiting  convec:  temp    saln    dens    thkns    dpth',
cdiag.  (k,temp(itest,jtest,k,n),saln(itest,jtest,k,n),
cdiag.  th3d(itest,jtest,k,n)+thbase,dp(itest,jtest,k,n)/onem,
cdiag.  p(itest,jtest,k+1)/onem,k=1,kk)
c
      return
      end
c
c
c> Revision history:
c>
c> Oct. 1999 - convection of u and v added
      subroutine convcm(m,n)
      use mod_xc  ! HYCOM communication interface
c
c --- hycom version 1.0 (adapted from micom version 2.8)
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
      integer i,j,k,l
      real    dthet,delp
c
      include 'stmt_fns.h'
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,3p,f8.3,0p,f8.2,f8.1))
cdiag if (itest.gt.0 .and. jtest.gt.0) write (lp,103) nstep,itest,jtest,
cdiag.  '  entering convec:  temp    saln    dens    thkns    dpth',
cdiag.  (k,temp(itest,jtest,k,n),saln(itest,jtest,k,n),
cdiag.  th3d(itest,jtest,k,n)+thbase,dp(itest,jtest,k,n)/onem,
cdiag.  p(itest,jtest,k+1)/onem,k=1,kk)
c
c --- -------------------------------------------------------------
c --- entrain water lighter than mixed-layer water into mixed layer
c --- -------------------------------------------------------------
c
!$OMP PARALLEL DO PRIVATE(j,k,l,i,dthet,delp)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 4 j=1-margin,jj+margin
c
      do 2 l=1,isp(j)
c
      do 9 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      klist(i,j)=0
      p(i,j,2)=dp(i,j,1,n)
 9    dpold(i,j,1)=dp(i,j,1,n)
c
      do 2 k=2,kk
      do 2 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
 2    dpold(i,j,k)=0.
c
      do 4 k=2,kk
c
      do 4 l=1,isp(j)
      do 4 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
      if (dp(i,j,k,n).le.0. .or.
     .   (klist(i,j).gt.0 .and. k.gt.klist(i,j)+1)) go to 4
      dthet=th3d(i,j,k,n)-th3d(i,j,1,n)
      if (dthet.lt.0. .and. p(i,j,k+1).le.p(i,j,kk+1)) then
c
cdiag   if (i.eq.itest .and. j.eq.jtest)
cdiag.  write (lp,100) nstep,i,j,' convec',1,th3d(i,j,1,n)
cdiag.  +thbase,dp(i,j,1,n)/onem,temp(i,j,1,n),saln(i,j,1,n),k,
cdiag.  th3d(i,j,k,n)+thbase,dp(i,j,k,n)/onem,temp(i,j,k,n),saln(i,j,k,n)
 100    format (i9,2i5,a,i3,'  th3d,dp,t,s =',3pf7.3,0pf7.1,2f8.3
     .    /26x,i3,15x,3pf7.3,0pf7.1,2f8.3)
c
c --- layer -k- contains mass less dense than mixed layer. entrain it.
        delp=dp(i,j,1,n)+dp(i,j,k,n)
        saln(i,j,1,n)=(saln(i,j,1,n)*dp(i,j,1,n)
     .                +saln(i,j,k, n)*dp(i,j,k, n))/delp
        temp(i,j,1,n)=(temp(i,j,1,n)*dp(i,j,1,n)
     .                +temp(i,j,k, n)*dp(i,j,k, n))/delp
        th3d(i,j,1,n)=sig(temp(i,j,1,n),saln(i,j,1,n))-thbase
c
        diaflx(i,j,1)=diaflx(i,j,1)+dp(i,j,k,n)			! diapyc.flux
        diaflx(i,j,k)=diaflx(i,j,k)-dp(i,j,k,n)			! diapyc.flux
c
c --- mass in layer -k- transferred to mixed layer is stored in -dpold-.
        dp(i,j,1,n)=delp
        dpold(i,j,k)=dp(i,j,k,n)
        dp(i,j,k,n)=0.
        klist(i,j)=k
      end if                 !  dthet < 0
c
 4    continue
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j,k,l,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do 7 j=1-margin,jj+margin
c
      do 20 l=1,isp(j)
      do 20 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      dpmixl(i,j,n)=dp(i,j,1,n)
 20   continue
c
c --- entrain -u- momentum
c
      do 5 l=1,isu(j)
c
      do 6 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
 6    util2(i,j)=min(.5*(dpold(i,j,1)+dpold(i-1,j,1)),depthu(i,j))
c
      do 5 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
      do 5 k=2,max(klist(i,j),klist(i-1,j))
      util1(i,j)=util2(i,j)
      util2(i,j)=min(util2(i,j)+
     .               .5*(dpold(i,j,k)+dpold(i-1,j,k)),depthu(i,j))
 5    u(i,j,1,n)=(u(i,j,1,n)*util1(i,j)
     .           +u(i,j,k,n)*(util2(i,j)-util1(i,j)))/util2(i,j)
c
c --- entrain -v- momentum
c
      do 7 l=1,isv(j)
c
      do 8 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
 8    util2(i,j)=min(.5*(dpold(i,j,1)+dpold(i,j-1,1)),depthv(i,j))
c
      do 7 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
      do 7 k=2,max(klist(i,j),klist(i,j-1))
      util1(i,j)=util2(i,j)
      util2(i,j)=min(util2(i,j)+
     .               .5*(dpold(i,j,k)+dpold(i,j-1,k)),depthv(i,j))
      v(i,j,1,n)=(v(i,j,1,n)*util1(i,j)
     .           +v(i,j,k,n)*(util2(i,j)-util1(i,j)))/util2(i,j)
 7    continue
!$OMP END PARALLEL DO
c
cdiag if (itest.gt.0 .and. jtest.gt.0) write (lp,103) nstep,itest,jtest,
cdiag.  '  exiting  convec:  temp    saln    dens    thkns    dpth',
cdiag.  (k,temp(itest,jtest,k,n),saln(itest,jtest,k,n),
cdiag.  th3d(itest,jtest,k,n)+thbase,dp(itest,jtest,k,n)/onem,
cdiag.  p(itest,jtest,k+1)/onem,k=1,kk)
c
      return
      end
c
c
c> Revision history:
c>
c> July 1997 - deleted final pressure calculation loop (moved to to thermf.f)
c> Apr. 1999 - added calculation of -th3d-
c> Aug. 2000 - adapted from micom 2.8 to run within hycom 1.0
