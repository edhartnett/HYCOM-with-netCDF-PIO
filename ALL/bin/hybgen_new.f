      subroutine hybgen(temp,saln,th3d,dp,theta,kdm,
     &                  nhybrd,hybflg,
     &                  dp00i,dp0k,dp0kp,ds0k,dssk,depths, ldebug)
      implicit none
c
      logical ldebug
      integer kdm, nhybrd,hybflg
      real     temp(1,1,kdm,1),
     &         saln(1,1,kdm,1),
     &         th3d(1,1,kdm,1),
     &           dp(1,1,kdm,1),
     &        theta(1,1,kdm)
      real    dp00i,dp0k(kdm),dp0kp(kdm),ds0k(kdm),dssk(kdm),depths(1,1)
c
      integer     mnproc,lp
      common/xxx/ mnproc,lp
c
c --- ---------------------
c --- hybrid grid generator
c --- ---------------------
c
      logical, parameter :: lunmix=.true.     !unmix a too light deepest layer
      logical, parameter :: lpcm=.false.      !PCM, instead of PLM, remapping
      logical, parameter :: lconserve=.false. !explicitly conserve each column
c
      integer i,j,kk,n, lp,nstep, i0,j0,itest,jtest
      real    epsil,onem,tencm,onemm,qonem,thbase
      real    p(1,1,kdm+1)
c
      double precision tsum,ssum,thsum,psum,rpsum
      double precision asum(  3,3)
      real             offset(3)
      real ttem(kdm+1,2),tsal(kdm+1,2),tthe(kdm+1,2),
     &     dprs(kdm+1),pres(kdm+2),dp0ij(kdm),dp0cum(kdm+1)
      real pwidth,p_hat,p_hat0,p_hat2,p_hat3,
     &     delt,deltm,dels,delsm,thnew,q,qtr,qts,thkbop,
     &     qbot,qcen,qtop,zbot,zbox,zcen,ztop,zthk,dpthin
      integer k,k1,ka,kbot,kbox,ktop,ksubl,kp,kq,ktr,l,iter
ccc   real colint,colins,colout,colous
      character*12 cinfo
c
      double precision, parameter :: dsmll=1.0d-8
      double precision, parameter ::   zp5=0.5    !for sign function
c
c --- c u s h i o n   function (from Bleck & Benjamin, 1992):
c --- if delp >= qqmx*dp0 >>  dp0, -cushn- returns -delp-
c --- if delp <= qqmn*dp0 << -dp0, -cushn- returns  -dp0-
c
      real       qqmn,qqmx,cusha,cushb
      parameter (qqmn=-4.0, qqmx=2.0)  ! shifted range
*     parameter (qqmn=-2.0, qqmx=4.0)  ! traditional range
*     parameter (qqmn=-4.0, qqmx=6.0)  ! somewhat wider range
      parameter (cusha=qqmn**2*(qqmx-1.0)/(qqmx-qqmn)**2)
      parameter (cushb=1.0/qqmn)
c
      external sig,dsigdt,dsigds,tofsig,sofsig
      real     sig,dsigdt,dsigds,tofsig,sofsig
      real     r,s,t
c
      real qq,cushn,delp,dp0
*     include 'stmt_fns.h'
      qq(   delp,dp0)=max(qqmn, min(qqmx, delp/dp0))
      cushn(delp,dp0)=dp0*
     &                (1.0+cusha*(1.0-cushb*qq(delp,dp0))**2)*
     &                max(1.0, delp/(dp0*qqmx))
c
      i      = 1
      j      = 1
      kk     = kdm
      n      = 1
      nstep  = 1
      thbase = 0.0
      i0     = 0
      j0     = 0
      if     (ldebug) then
        itest = 1
        jtest = 1
      else
        itest = 0
        jtest = 0
      endif
c
      epsil = 1.0e-11
      onem  = 9806.0 !g/thref
      qonem = 1.0/onem
      tencm = onem/10.0
      onemm = onem/1000.0
c
      dpthin = 0.001*onemm
*     thkbop = thkbot*onem
      thkbop =   10.0*onem
c
*     do 1 l=1,isp(j)
c
*     do 2 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
      dp0cum(1)=0.0
      dp0ij( 1)=min( dp0k(1), max( ds0k(1), dssk(1)*depths(i,j) ) )
      dp0cum(2)=dp0cum(1)+dp0ij(1)
      p(i,j, 2)=p(i,j,1)+dp(i,j,1,n)
      do k=2,kk
        if     (dp0k(k).le.dp00i) then
          q =      dp0k(k)
        else
          q = max( dp00i,
     &             dp0k(k)*
     &             (dp0k(k)/max(dp0k(k),p(i,j,k)-dp0cum(k)))**2 )
        endif
        dp0ij( k)  =min( q,max( ds0k(k), dssk(k)*depths(i,j) ) )
        dp0cum(k+1)=dp0cum(k)+dp0ij(k)
        p(i,j, k+1)=p(i,j,k)+dp(i,j,k,n)
      enddo
c
c --- mark the layer containing the mixed layer base by -ksubl- if the
c --- kraus-turner mixed layer model is selected. otherwise set ksubl to -1.
*     if     (mxlkta .and. thermo) then
*       ksubl=-1
*       do k=1,kk
*         if (p(i,j,k  ).lt.dpmixl(i,j,n) .and. 
*    &        p(i,j,k+1).ge.dpmixl(i,j,n)+onemm) then
*           ksubl=k
*         endif
*       enddo
*       klist(i,j)=ksubl
*     else
        ksubl=-1
*     endif
c
c --- prevent too many extremely thin layers from existing at the bottom.
c --- if one extremely thin layer exists with a density that is not more than
c --- 0.0001 sigma units less than the target density, it is left alone. if
c --- this is not true, or if two or more extremely thin layers exist, the
c --- thin layers are absorbed into the first layer above that has significant
c --- thickness (> 10 cm).
c
      ka=0
      do k=3,kk
        if (k.le.nhybrd) then
          if (p(i,j,k+1).eq.p(i,j,kk+1) .and.
     &        dp(i,j,k,n).ge.tencm) then
            ka=k   ! deepest thick layer, no thin layers at bottom
            exit
          elseif (p(i,j,k).gt.p(i,j,kk+1)-tencm) then
            if (abs(th3d(i,j,k,n)-theta(i,j,k)).gt.0.0001 .or.
     &          p(i,j,k+1).lt.p(i,j,kk+1)) then
              ka=k-1   ! deepest thick layer
c ---         absorb thin layer(s) into deepest thick layer
              psum=dp(i,j,ka,n)
              if     (hybflg.eq.0) then
                 tsum=temp(i,j,ka,n)*psum
                 ssum=saln(i,j,ka,n)*psum
              elseif (hybflg.eq.1) then
                thsum=th3d(i,j,ka,n)*psum
                 ssum=saln(i,j,ka,n)*psum
              elseif (hybflg.eq.2) then
                thsum=th3d(i,j,ka,n)*psum
                 tsum=temp(i,j,ka,n)*psum
              endif
              do k1=k,kk
                if (p(i,j,k1).lt.p(i,j,kk+1)) then
                  psum=psum+dp(i,j,k1,n)
                  if     (hybflg.eq.0) then
                     tsum= tsum+temp(i,j,k1,n)*dp(i,j,k1,n)
                     ssum= ssum+saln(i,j,k1,n)*dp(i,j,k1,n)
                  elseif (hybflg.eq.1) then
                    thsum=thsum+th3d(i,j,k1,n)*dp(i,j,k1,n)
                     ssum= ssum+saln(i,j,k1,n)*dp(i,j,k1,n)
                  elseif (hybflg.eq.2) then
                    thsum=thsum+th3d(i,j,k1,n)*dp(i,j,k1,n)
                     tsum= tsum+temp(i,j,k1,n)*dp(i,j,k1,n)
                  endif
                endif
                if (th3d(i,j,k1,n).ne.theta(i,j,k1)) then
                  th3d(i,j,k1,n)=theta(i,j,k1)
                  if     (hybflg.eq.0) then
                    saln(i,j,k1,n)=saln(i,j,k1-1,n)
                    temp(i,j,k1,n)=tofsig(th3d(i,j,k1,n)+thbase,
     &                                    saln(i,j,k1,n))
                  elseif (hybflg.eq.1) then
                    saln(i,j,k1,n)=saln(i,j,k1-1,n)
                    temp(i,j,k1,n)=tofsig(th3d(i,j,k1,n)+thbase,
     &                                    saln(i,j,k1,n))
                  elseif (hybflg.eq.2) then
                    temp(i,j,k1,n)=temp(i,j,k1-1,n)
                    saln(i,j,k1,n)=sofsig(th3d(i,j,k1,n)+thbase,
     &                                    temp(i,j,k1,n))
                  endif
                endif
                p(i,j,k1)=p(i,j,kk+1)
              enddo
              if     (hybflg.eq.0) then
                temp(i,j,ka,n)= tsum/psum
                saln(i,j,ka,n)= ssum/psum
                th3d(i,j,ka,n)=sig(temp(i,j,ka,n),saln(i,j,ka,n))-thbase
              elseif (hybflg.eq.1) then
                th3d(i,j,ka,n)=thsum/psum
                saln(i,j,ka,n)= ssum/psum
                temp(i,j,ka,n)=tofsig(th3d(i,j,ka,n)+thbase,
     &                                saln(i,j,ka,n))
              elseif (hybflg.eq.2) then
                th3d(i,j,ka,n)=thsum/psum
                temp(i,j,ka,n)= tsum/psum
                saln(i,j,ka,n)=sofsig(th3d(i,j,ka,n)+thbase,
     &                                temp(i,j,ka,n))
              endif
            endif
            exit  ! deepest thick layer
          endif
        endif !k<=nhybrd
      enddo !k
c
c --- if ka>0, the deepest inflated layer has significant thickness (>10 cm)
      if     (ka.gt.0) then
      k=ka !at least 2
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,i3)')
     &        'hybgen, deepest inflated layer:',k
        call flush(lp)
      endif !debug
      if     ((theta(i,j,k) -th3d(i,j,k,  n)).gt. 0.0001 .and.
     &        (th3d(i,j,k,n)-th3d(i,j,k-1,n)).lt.-0.0001      ) then
c
c ---   water in the deepest inflated layer with significant thickness
c ---   is too light, and is lighter than the layer above (unstable).
c
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3)')
     &      'hybgen, deepest inflated layer too light (unstable):',k
          call flush(lp)
        endif !debug
        if     (th3d(i,j,k,n).le.theta(i,j,k-1)) then
c ---     adsorb the (very light) bottom layer into the layer above
          if     (hybflg.eq.0) then
             psum=                dp(i,j,k-1,n) +
     &                            dp(i,j,k,  n)
             tsum=temp(i,j,k-1,n)*dp(i,j,k-1,n) +
     &            temp(i,j,k,  n)*dp(i,j,k,  n)
             ssum=saln(i,j,k-1,n)*dp(i,j,k-1,n) +
     &            saln(i,j,k,  n)*dp(i,j,k,  n)
            temp(i,j,k-1,n)= tsum/psum
            saln(i,j,k-1,n)= ssum/psum
            th3d(i,j,k-1,n)=sig(temp(i,j,k-1,n),saln(i,j,k-1,n))-thbase
          elseif (hybflg.eq.1) then
             psum=                dp(i,j,k-1,n) +
     &                            dp(i,j,k,  n)
            thsum=th3d(i,j,k-1,n)*dp(i,j,k-1,n) +
     &            th3d(i,j,k,  n)*dp(i,j,k,  n)
             ssum=saln(i,j,k-1,n)*dp(i,j,k-1,n) +
     &            saln(i,j,k,  n)*dp(i,j,k,  n)
            th3d(i,j,k-1,n)=thsum/psum
            saln(i,j,k-1,n)= ssum/psum
            temp(i,j,k-1,n)=tofsig(th3d(i,j,k-1,n)+thbase,
     &                             saln(i,j,k-1,n))
          elseif (hybflg.eq.2) then
             psum=                dp(i,j,k-1,n) +
     &                            dp(i,j,k,  n)
            thsum=th3d(i,j,k-1,n)*dp(i,j,k-1,n) +
     &            th3d(i,j,k,  n)*dp(i,j,k,  n)
             tsum=temp(i,j,k-1,n)*dp(i,j,k-1,n) +
     &            temp(i,j,k,  n)*dp(i,j,k,  n)
            th3d(i,j,k-1,n)=thsum/psum
            temp(i,j,k-1,n)= tsum/psum
            saln(i,j,k-1,n)=sofsig(th3d(i,j,k-1,n)+thbase,
     &                             temp(i,j,k-1,n))
          endif !hybflg
             p(i,j,k)    =p(i,j,k+1)
            dp(i,j,k-1,n)=psum
            dp(i,j,k,  n)=0.0
          th3d(i,j,k,  n)=th3d(i,j,k-1,n)
          temp(i,j,k,  n)=temp(i,j,k-1,n)
          saln(i,j,k,  n)=saln(i,j,k-1,n)
          if (i.eq.itest .and. j.eq.jtest) then
            write(lp,'(a,i3,f6.3,5f8.3)')
     &        'hybgen, adsor:',
     &        k,0.0,temp(i,j,k,n),saln(i,j,k,n),
     &            th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
            call flush(lp)
          endif !debug
        else
c ---     entrain up to 50% of near-bottom layer into bottom layer
          if     ((th3d(i,j,k-1,n)-theta(i,j,k)).le.epsil) then
            q=0.5*dp(i,j,k-1,n)
          else
            q=(th3d(i,j,k,  n)-theta(i,j,k))/
     &        (th3d(i,j,k-1,n)-theta(i,j,k))  ! -ve
            q=min( 0.5*dp(i,j,k-1,n),
     &              -q*dp(i,j,k,  n) )
          endif
          if     (hybflg.eq.0) then
             psum=                q +
     &                            dp(i,j,k,  n)
             tsum=temp(i,j,k-1,n)*q +
     &            temp(i,j,k,  n)*dp(i,j,k,  n)
             ssum=saln(i,j,k-1,n)*q +
     &            saln(i,j,k,  n)*dp(i,j,k,  n)
            temp(i,j,k,n)= tsum/psum
            saln(i,j,k,n)= ssum/psum
            th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
          elseif (hybflg.eq.1) then
             psum=                q +
     &                            dp(i,j,k,  n)
            thsum=th3d(i,j,k-1,n)*q +
     &            th3d(i,j,k,  n)*dp(i,j,k,  n)
             ssum=saln(i,j,k-1,n)*q +
     &            saln(i,j,k,  n)*dp(i,j,k,  n)
            th3d(i,j,k,n)=thsum/psum
            saln(i,j,k,n)= ssum/psum
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,
     &                           saln(i,j,k,n))
          elseif (hybflg.eq.2) then
             psum=                q +
     &                            dp(i,j,k,  n)
            thsum=th3d(i,j,k-1,n)*q +
     &            th3d(i,j,k,  n)*dp(i,j,k,  n)
             tsum=temp(i,j,k-1,n)*q +
     &            temp(i,j,k,  n)*dp(i,j,k,  n)
            th3d(i,j,k,n)=thsum/psum
            temp(i,j,k,n)= tsum/psum
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
     &                           temp(i,j,k,n))
          endif !hybflg
           p(i,j,k)    = p(i,j,k)    -q
          dp(i,j,k-1,n)=dp(i,j,k-1,n)-q
          dp(i,j,k,  n)=psum
          if (i.eq.itest .and. j.eq.jtest) then
            write(lp,'(a,i3,f6.3,5f8.3)')
     &        'hybgen, unstb:',
     &        k,psum*qonem,temp(i,j,k,n),saln(i,j,k,n),
     &            th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
            call flush(lp)
          endif !debug
        endif !adsorb bottom layer:entrain near-bottom layer
      endif !too light and unstable
      if     (lunmix .and. !usually .true.
     &        (theta(i,j,k) -th3d(i,j,k,  n)).gt.0.0001 .and.
     &        (th3d(i,j,k,n)-th3d(i,j,k-1,n)).gt.0.0001      ) then
c
c ---   water in the deepest inflated layer with significant thickness
c ---   is too light. split layer into 2 sublayers, one near the desired
c ---   density and one exactly matching the T&S properties of layer k-1.
c ---   To prevent "runaway" T or S, the result satisfies either
c ---     abs(T.k - T.k-1) <= abs(T.k-2 - T.k-1) or
c ---     abs(S.k - S.k-1) <= abs(S.k-2 - S.k-1) 
c ---   It is also limited to a 50% change in layer thickness.
c
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3)')
     &      'hybgen, deepest inflated layer too light   (stable):',k
          call flush(lp)
        endif !debug
c
        delsm=abs(saln(i,j,k-2,n)-saln(i,j,k-1,n))
        dels =abs(saln(i,j,k-1,n)-saln(i,j,k,  n))
        deltm=abs(temp(i,j,k-2,n)-temp(i,j,k-1,n))
        delt =abs(temp(i,j,k-1,n)-temp(i,j,k,  n))
c ---   sanity check on deltm and delsm
        q=min(temp(i,j,k-2,n),temp(i,j,k-1,n),temp(i,j,k,n))
        if     (q.gt. 6.0) then
          deltm=min( deltm,  6.0*(theta(i,j,k)-theta(i,j,k-1)) )
        elseif (q.gt. 0.0) then
          deltm=min( deltm, 10.0*(theta(i,j,k)-theta(i,j,k-1)) )
        else  !(q.le. 0.0)
          deltm=min( deltm, 25.0*(theta(i,j,k)-theta(i,j,k-1)) )
        endif
        delsm=min( delsm, 1.3*(theta(i,j,k)-theta(i,j,k-1)) )
        qts=0.0
        if     (delt.gt.epsil) then
          qts=max(qts, (min(deltm, 2.0*delt)-delt)/delt)  ! qts<=1.0
        endif
        if     (dels.gt.epsil) then
          qts=max(qts, (min(delsm, 2.0*dels)-dels)/dels)  ! qts<=1.0
        endif
        q=(theta(i,j,k)-th3d(i,j,k,  n))/
     &    (theta(i,j,k)-th3d(i,j,k-1,n))
        q=min(q,qts/(1.0+qts))  ! upper sublayer <= 50% of total
        p_hat=q*(p(i,j,k+1)-p(i,j,k))
        p(i,j,k)=p(i,j,k)+p_hat
c
        if     (hybflg.eq.0) then  !T&S
          temp(i,j,k,n)=temp(i,j,k,n)+(q/(1.0-q))*(temp(i,j,k,n)  -
     &                                             temp(i,j,k-1,n) )
          saln(i,j,k,n)=saln(i,j,k,n)+(q/(1.0-q))*(saln(i,j,k,n)  -
     &                                             saln(i,j,k-1,n) )
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
        elseif (hybflg.eq.1) then  !th&S
          th3d(i,j,k,n)=th3d(i,j,k,n)+(q/(1.0-q))*(th3d(i,j,k,n)  -
     &                                             th3d(i,j,k-1,n) )
          saln(i,j,k,n)=saln(i,j,k,n)+(q/(1.0-q))*(saln(i,j,k,n)  -
     &                                             saln(i,j,k-1,n) )
          temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
        elseif (hybflg.eq.2) then  !th&T
          th3d(i,j,k,n)=th3d(i,j,k,n)+(q/(1.0-q))*(th3d(i,j,k,n)  -
     &                                             th3d(i,j,k-1,n) )
          temp(i,j,k,n)=temp(i,j,k,n)+(q/(1.0-q))*(temp(i,j,k,n)  -
     &                                             temp(i,j,k-1,n) )
          saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
        endif
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3,f6.3,5f8.3)')
     &      'hybgen, 10(+):',
     &      k,q,temp(i,j,k,n),saln(i,j,k,n),
     &          th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
          call flush(lp)
        endif !debug
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3,f6.3,5f8.3)')
     &      'hybgen, 10(-):',
     &      k,0.0,temp(i,j,k,n),saln(i,j,k,n),
     &          th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
          call flush(lp)
        endif !debug
      endif !too light and stable
      endif !ka>0
c
c --- store one-dimensional arrays of -temp-, -saln-, and -p- for the 'old'
c --- vertical grid before attempting to restore isopycnic conditions
      pres(1)=p(i,j,1)
      kp=0
      do k=1,kk
        k1=k+kp
*       if (k.ne.ksubl) then
          tthe(k1,1)=th3d(i,j,k,n)
          ttem(k1,1)=temp(i,j,k,n)
          tsal(k1,1)=saln(i,j,k,n)
          pres(k1+1)=p(i,j,k+1)
          dprs(k1)  =pres(k1+1)-pres(k1)
*       else				!  k = ksubl
c ---     expand layer into two sublayers, above and below mixed layer base
*         tthe(k,1)=tthe(k-1,1)
*         ttem(k,1)=ttem(k-1,1)
*         tsal(k,1)=tsal(k-1,1)
*         pres(k+1)=dpmixl(i,j,n)
*         dprs(k)  =pres(k+1)-pres(k)
*         q=(p(i,j,k)-dpmixl(i,j,n))/(p(i,j,k+1)-dpmixl(i,j,n))
*         tthe(k+1,1)=th3d(i,j,k,n)+q*(tthe(k,1)-th3d(i,j,k,n))
*         ttem(k+1,1)=temp(i,j,k,n)+q*(ttem(k,1)-temp(i,j,k,n))
*         tsal(k+1,1)=saln(i,j,k,n)+q*(tsal(k,1)-saln(i,j,k,n))
*         pres(k+2)=p(i,j,k+1)
*         dprs(k+1)=pres(k+2)-pres(k+1)
*         kp=1
*       endif !k.ne.ksubl:else
      enddo !k
c
      if     (lpcm) then  !usually .false.
c
c ---   PCM (zero slope, recovers original hybgen behaviour).
c
        do k= 1,kk+kp
          tthe(k,2)=0.0
          ttem(k,2)=0.0
          tsal(k,2)=0.0
        enddo
      else
c
c ---   PLM (non-zero slope, but no new extrema)
c ---   layer value is (:,1)-0.5*(:,2) at top    interface,
c ---              and (:,1)+0.5*(:,2) at bottom interface.
c ---   still use PCM for isopycnal layers, because we don't
c ---   want detrainment to change the density of these layers.
c
c ---   monotonized central-difference limiter (van Leer, 1977,
c ---   JCP 23 pp 276-299).  For a discussion of PLM limiters, see
c ---   Finite Volume Methods for Hyperbolic Problems by R.J. Leveque.
c
        do k= 1,kk+kp
          if     (k.eq.1 .or. k.eq.kk+kp .or.
     &            dprs(k).le.dpthin      .or.
     &            dprs(k).gt.dp0kp(k)        ) then
c ---       top, bottom, thin and isopycnal layers have zero slope.
            tthe(k,2)=0.0
            ttem(k,2)=0.0
            tsal(k,2)=0.0
          else
c ---       interior non-isopycnal layer
c ---       use qcen in place of 0.5 to allow for non-uniform grid
            qcen = dprs(k)/(dprs(k)+0.5*(dprs(k-1)+dprs(k+1)))
c
            ztop = 2.0*(tthe(k,  1)-tthe(k-1,1))
            zbot = 2.0*(tthe(k+1,1)-tthe(k,  1))
            zcen =qcen*(tthe(k+1,1)-tthe(k-1,1))
            if     (ztop*zbot.gt.0.0) then
              tthe(k,2)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
            else
              tthe(k,2)=0.0  !local extrema, so no slope
            endif
c
            ztop = 2.0*(ttem(k,  1)-ttem(k-1,1))
            zbot = 2.0*(ttem(k+1,1)-ttem(k,  1))
            zcen =qcen*(ttem(k+1,1)-ttem(k-1,1))
            if     (ztop*zbot.gt.0.0) then
              ttem(k,2)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
            else
              ttem(k,2)=0.0  !local extrema, so no slope
            endif
c
            ztop = 2.0*(tsal(k,  1)-tsal(k-1,1))
            zbot = 2.0*(tsal(k+1,1)-tsal(k,  1))
            zcen =qcen*(tsal(k+1,1)-tsal(k-1,1))
            if     (ztop*zbot.gt.0.0) then
              tsal(k,2)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
            else
              tsal(k,2)=0.0  !local extrema, so no slope
            endif
c
          endif !top/bottom/thin/isopycnal or PLM layer
        enddo !k
      endif !PCM:PLM
c
ccc   colint=temp(i,j,1,n)*(p(i,j,2)-p(i,j,1))
ccc   colins=saln(i,j,1,n)*(p(i,j,2)-p(i,j,1))
c
c --- try to restore isopycnic conditions by moving layer interfaces
c
      do 88 k=2,nhybrd
            if (i.eq.itest .and. j.eq.jtest) then
              write(cinfo,'(a9,i2.2,1x)') '  do 88 k=',k
 109          format (i9,2i5,a,a/(i9,8x,a,a,i3,f9.2,f8.2,f9.2,f8.2))
              write(lp,109) nstep,itest+i0,jtest+j0,
     .         cinfo,':    othkns  odpth    nthkns  ndpth',
     .        (nstep,cinfo,':',k1,
     .         (pres(k1+1)-
     .          pres(k1)   )*qonem,
     .          pres(k1+1)  *qonem,
     .         (p(itest,jtest,k1+1)-
     .          p(itest,jtest,k1)   )*qonem,
     .          p(itest,jtest,k1+1)  *qonem,k1=1,kk)
              call flush(lp)
            endif !debug
c
ccc   colint=colint+temp(i,j,k,n)*(p(i,j,k+1)-p(i,j,k))
ccc   colins=colins+saln(i,j,k,n)*(p(i,j,k+1)-p(i,j,k))
c
c --- maintain constant thickness in layer 1
      if (k.eq.2) then
        p_hat=-999.0*dp0k(1)
        go to 9
      end if
c
c --- are we dealing with a massless layer on the sea floor?
      if (p(i,j,k).gt.p(i,j,kk+1)-dpthin) then
        if (k.le.nhybrd) then
c ---     fill massless layers on sea floor with fluid from above
          th3d(i,j,k,n)=th3d(i,j,k-1,n)
          saln(i,j,k,n)=saln(i,j,k-1,n)
          temp(i,j,k,n)=temp(i,j,k-1,n)
        elseif (th3d(i,j,k,n).ne.theta(i,j,k)) then
          if (hybflg.ne.2) then
c ---       fill with saln from above
            th3d(i,j,k,n)=max(theta(i,j,k), th3d(i,j,k-1,n))
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          else
c ---       fill with temp from above
            th3d(i,j,k,n)=max(theta(i,j,k), th3d(i,j,k-1,n))
            temp(i,j,k,n)=temp(i,j,k-1,n)
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
        endif
        go to 8
      endif
c
c --- is density noticeably different from isopycnic reference value?
      if (abs(th3d(i,j,k,n)-theta(i,j,k)).lt.0.0001) go to 8
c
      if (th3d(i,j,k,n).le.theta(i,j,k)) go to 7	!  layer too light
c
c --- water in layer k is too dense. try to dilute with water from layer k-1
c
c --- if layer k-1 is too light, perform this adjustment every other pair
c --- of time steps beneath the z-coordinate domain. otherwise, layers
c --- k-1 and k may not be able to reach their target densities
      if (mod(nstep-1,4).le.1 .and.
     &    p(i,j,k).gt.dp0cum(k)+onem .and.
     &    th3d(i,j,k-1,n).lt.theta(i,j,k-1)) go to 8
c
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,i3,f8.4)') 'hybgen, too dense:',
     &                          k,th3d(i,j,k,n)-theta(i,j,k)
        call flush(lp)
      endif !debug
c
      if     ((theta(i,j,k)-th3d(i,j,k-1,n)).le.epsil) then
        p_hat=min(p(i,j,k-1),
     &            p(i,j,k)-999.0*dp0k(k-1))  ! take entire layer k-1
      else
        q=(theta(i,j,k)-th3d(i,j,k,n))/(theta(i,j,k)-th3d(i,j,k-1,n))  ! -ve
        p_hat=p(i,j,k)+q*(p(i,j,k+1)-p(i,j,k))
c
c ---   correct for cabbeling by performing an iterative procedure to make
c ---   sure that the target density is achieved to a tolerance of 0.0001.
c ---   if the water becomes more dense, do not move the interface
c
        if     (hybflg.eq.0 .and.
     &          p_hat.gt.p(i,j,k-1) .and. 
     &          p(i,j,k).gt.dp0cum(k)+onem) then
          do 58 iter=1,5
            q=(p(i,j,k)-p_hat)/max(p(i,j,k+1)-p_hat,onemm)
            delt=q*(temp(i,j,k-1,n)-temp(i,j,k,n))
            dels=q*(saln(i,j,k-1,n)-saln(i,j,k,n))
            thnew=sig(temp(i,j,k,n)+delt,saln(i,j,k,n)+dels)-thbase
            if (thnew.gt.th3d(i,j,k,n)) then  ! more dense
              p_hat=p(i,j,k)
              go to 158
            endif
            if (abs(thnew-theta(i,j,k)).le.0.0001) then
              go to 158
            endif
            q=(theta(i,j,k)-thnew)/(theta(i,j,k)-th3d(i,j,k-1,n))
            p_hat=min(p(i,j,k),p_hat+q*(p(i,j,k+1)-p_hat))
  58      continue
 158      continue
        end if
      end if
c
c --- maintain minimum layer thickess of layer k-1.
c
 9    continue
      p_hat0=p_hat
      p_hat=max(p(i,j,k-1)+dp0ij(k-1),
     &      min(p(i,j,k+1)-dp0ij(k),
     &          p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1)) ))
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,3f8.2)') 'hybgen, 9: ',
     &   (p_hat0-p(i,j,k-1))*qonem,
     &   dp0ij(k-1)*qonem,
     &   cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))*qonem
        call flush(lp)
      endif !debug
c
c --- if isopycnic conditions cannot be achieved because of a blocking
c --- layer in the interior ocean, move interface k-1 (and k-2 if
c --- necessary) upward
c
      if     (k.le.2) then
c ---   do nothing.
      else if (p_hat.ge.p(i,j,k) .and. 
     &         p(i,j,k-1).gt.dp0cum(k-1)+onem .and.
     &         (p(i,j,kk+1)-p(i,j,k-1).lt.thkbop .or.
     &          p(i,j,k-1) -p(i,j,k-2).gt.qqmx*dp0ij(k-2))) then  ! k.gt.2
        p_hat2=p(i,j,k-2)+cushn(p(i,j,k-1)-p_hat+p_hat0-p(i,j,k-2),
     &                          dp0ij(k-2))
        if (p_hat2.lt.p(i,j,k-1)-onemm) then
          p(i,j,k-1)=max(p_hat2,2.0*p(i,j,k-1)-p_hat)
          if (i.eq.itest .and. j.eq.jtest) then
            write(lp,'(a,i3.2,f8.2)') 'hybgen,  1blocking :',
     .                              k-1,p(i,j,k-1)*qonem
            call flush(lp)
          endif !debug
          p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
        else if (k.le.3) then
c ---     do nothing.
        else if (p(i,j,k-2).gt.dp0cum(k-2)+onem .and.
     &           (p(i,j,kk+1)-p(i,j,k-2).lt.thkbop .or.
     &            p(i,j,k-2) -p(i,j,k-3).gt.qqmx*dp0ij(k-3))) then  ! k.gt.3
          p_hat3=p(i,j,k-3)+cushn(p(i,j,k-2)-p_hat+p_hat0-p(i,j,k-3),
     &                            dp0ij(k-3))
          if (p_hat3.lt.p(i,j,k-2)-onemm) then
            p(i,j,k-2)=max(p_hat3,2.0*p(i,j,k-2)-p(i,j,k-1))
            if (i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,i3.2,f8.2)') 'hybgen,  2blocking :',
     .                                k-2,p(i,j,k-2)*qonem
              call flush(lp)
            endif !debug
            p_hat2=p(i,j,k-2)+cushn(p(i,j,k-1)-p_hat+p_hat0-p(i,j,k-2),
     &                              dp0ij(k-2))
            if (p_hat2.lt.p(i,j,k-1)-onemm) then
              p(i,j,k-1)=max(p_hat2,2.0*p(i,j,k-1)-p_hat)
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i3.2,f8.2)') 'hybgen,  3blocking :',
     .                                  k-1,p(i,j,k-1)*qonem
                call flush(lp)
              endif !debug
              p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
            end if
          end if
        end if
      end if
c
      if (p_hat.le.p(i,j,k)) then
c
c --- upper intfc moves up. entrain layer k-1 water into layer k
c
        p(i,j,k)=p_hat
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k) :',
     .                            k,p(i,j,k)*qonem
          call flush(lp)
        endif !debug
c
      else				!  p_hat > p(i,j,k)
c
c --- move upper interface down and entrain layer k water into layer k-1
c
        p_hat=min(p_hat,p(i,j,k+1))
        p(i,j,k)=p_hat
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k-):',
     .                            k,p(i,j,k)*qonem
          call flush(lp)
        endif !debug
c
c --- do we need to inflate layer k by lowering  l o w e r  interface?
c
        if (k.lt.kk) then
          p_hat=p(i,j,k+1)
          go to 6
        end if
      end if
      go to 8
c
c --- water in layer k is too light. try to dilute with water from layer k+1
c
 7    continue
c
c --- is this layer touching the sea floor?
      if (p(i,j,k+1).gt.p(i,j,kk+1)-dpthin .and. k.gt.2) go to 8
c
c --- are we below any KT mixed layer?
*     if (mxlkta .and. k.le.klist(i,j)) go to 8
c
c --- if layer k+1 is too dense, perform this adjustment every other pair
c --- of time steps beneath the z-coordinate domain. otherwise, layers
c --- k and k+1 may not be able to reach their target densities
      if (mod(nstep-1,4).ge.2 .and.
     &    p(i,j,k+1).gt.dp0cum(k+1)+onem   .and.
     &    th3d(i,j,k+1,n).gt.theta(i,j,k+1)     ) go to 8
c
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,i3,f8.4)') 'hybgen, too light:',
     &                          k,theta(i,j,k)-th3d(i,j,k,n)
        call flush(lp)
      endif !debug
c
      if     ((th3d(i,j,k+1,n)-theta(i,j,k)).le.epsil) then
        p_hat=max(p(i,j,k+2),
     &            p(i,j,k+1)+999.0*dp0k(k))  ! take entire layer k+1
      else
        q=(th3d(i,j,k,n)-theta(i,j,k))/(th3d(i,j,k+1,n)-theta(i,j,k))  ! -ve
        p_hat=p(i,j,k+1)+q*(p(i,j,k)-p(i,j,k+1))
c
c ---   correct for cabbeling by performing an iterative procedure to make
c ---   sure that the target density is achieved to a tolerance of 0.0001.
c ---   if the water becomes less dense, do not move the interface
c
        if     (hybflg.eq.0 .and.
     &          p_hat.lt.p(i,j,k+2) .and. 
     &          p(i,j,k+1).gt.dp0cum(k+1)) then
          do 59 iter=1,5
            q=(p_hat-p(i,j,k+1))/max(p_hat-p(i,j,k),onemm)
            delt=q*(temp(i,j,k+1,n)-temp(i,j,k,n))
            dels=q*(saln(i,j,k+1,n)-saln(i,j,k,n))
            thnew=sig(temp(i,j,k,n)+delt,saln(i,j,k,n)+dels)-thbase
            if (thnew.lt.th3d(i,j,k,n)) then  ! less dense
              p_hat=p(i,j,k+1)
              go to 159
            endif
            if (abs(thnew-theta(i,j,k)).le.0.0001) then
              go to 159
            endif
            q=(thnew-theta(i,j,k))/(th3d(i,j,k+1,n)-theta(i,j,k))
            p_hat=max(p(i,j,k+1),p_hat+q*(p(i,j,k)-p_hat))
  59      continue
 159      continue
        endif
      endif
c
c --- maintain minimum thickess of layers k and k+1
 6    continue
*     if (mxlkta .and. k.eq.klist(i,j)) go to 8
      if     (k.lt.nhybrd) then
        p_hat=p(i,j,k+2)-cushn(p(i,j,k+2)-p_hat,dp0ij(k+1))
      endif
      p_hat=p(i,j,k)  +  max(p_hat    -p(i,j,k),dp0ij(k))
      p_hat=min(p_hat,
     &          max(0.5*(p(i,j,k+1)+p(i,j,k+2)),
     &              p(i,j,k+2)-dp0ij(k+1)))
c
      if (p_hat.gt.p(i,j,k+1)) then
c
c --- entrain layer k+1 water into layer k
c
        p(i,j,k+1)=p_hat
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k+):',
     .                            k+1,p(i,j,k+1)*qonem
          call flush(lp)
        endif !debug
      end if
 8    continue
c
c --- if layer above is too thin, move upper interface down
      p_hat0=p(i,j,k-1)+dp0ij(k-1)
      if (p_hat0.gt.p(i,j,k)   .and.
     &    p_hat0.lt.p(i,j,kk+1)     ) then
        p(i,j,k)=p_hat0
      endif
c
c --- inforce interface order (is this necessary?), usually inexpensive.
      do ka= k+1,kk
        if     (p(i,j,ka).ge.p(i,j,k)) then
          exit  ! usually get here quickly
        else
          p(i,j,ka) = p(i,j,k)
        endif
      enddo
 88   continue
c
ccc   write (lp,'(a,1p,2e16.8,e9.1)') 'temp. column integral:',
ccc  .  colint,colout,(colout-colint)/colint
ccc   write (lp,'(a,1p,2e16.8,e9.1)') 'saln. column integral:',
ccc  .  colins,colous,(colous-colins)/colins
c
c --- remap scalar field profiles from the 'old' vertical
c --- grid onto the 'new' vertical grid, using PLM
c
      if     (lconserve) then  !usually .false.
        do ktr=1,5
          asum(ktr,1) = 0.d0
          asum(ktr,2) = 0.d0
          asum(ktr,3) = 0.d0
        enddo !ktr
      endif
c
      zbot=0.0
      kbot=1
      do k=1,kk
        ztop=zbot  !top is bottom of previous layer
        ktop=kbot
        if     (ztop.ge.pres(ktop+1)) then
          ktop=ktop+1
        endif
       if (ktop.eq.kbot .and. i.eq.itest .and. j.eq.jtest) then
         write(lp,'(a,2i3,2f23.16)')
     &          'k,ktop,dp =',k,ktop,
     &          ztop*qonem,(pres(ktop+1)-ztop)*qonem
       endif !debug
        zbot=p(i,j,k+1)
        zthk=zbot-ztop
        dp(i,j,k,n)=zthk
            if (i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,2i3,3f23.16)')
     &          'k,z =',k,k,ztop*qonem,zbot*qonem,zthk*qonem
            endif !debug
        if     (zthk.gt.dpthin .and. ztop.lt.p(i,j,kk+1)) then
c ---     normal layer
          kbot=ktop
          do while (pres(kbot+1).lt.zbot.and.kbot.lt.kk+kp)
            kbot=kbot+1
          enddo
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,2i3,2f23.16)')
     &            'k,p =',ktop,kbot,
     &            pres(ktop)*qonem,pres(kbot+1)*qonem
              endif !debug
c
c ---     include thin adjacent layers in sum
          zbox=zbot
          do k1= k+1,kk
            if     (p(i,j,k1+1)-p(i,j,k1).gt.dpthin) then
              exit !thick layer
            else
              zbox=p(i,j,k1+1)  !include thin adjacent layers
              if     (zbox.eq.p(i,j,kk+1)) then
                exit !at bottom
              endif
            endif
          enddo
          zthk=zbox-ztop
c
          kbox=ktop
          do while (pres(kbox+1).lt.zbox.and.kbox.lt.kk+kp)
            kbox=kbox+1
          enddo
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,2i3,2f23.16)')
     &            'k,px=',ktop,kbox,
     &            pres(ktop)*qonem,pres(kbox+1)*qonem
              endif !debug
          if     (ktop.eq.kbox) then
c ---       single layer
            if     (p(i,j,k)  .ne.pres(kbox)   .or.
     &              p(i,j,k+1).ne.pres(kbox+1)     ) then
c ---         part of a single layer
              zcen = 0.5*(ztop+zbox)
              if     (dprs(kbox).gt.dpthin) then
                q = 0.5 - (pres(kbox+1)-zcen)/dprs(kbox)
              else
                q = 0.0
              endif
              if     (hybflg.eq.0) then  !T&S
                temp(i,j,k,n)=ttem(kbox,1)+q*ttem(kbox,2)
                saln(i,j,k,n)=tsal(kbox,1)+q*tsal(kbox,2)
                th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
              elseif (hybflg.eq.1) then  !th&S
                th3d(i,j,k,n)=tthe(kbox,1)+q*tthe(kbox,2)
                saln(i,j,k,n)=tsal(kbox,1)+q*tsal(kbox,2)
                temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,
     &                               saln(i,j,k,n))
              elseif (hybflg.eq.2) then  !th&T
                th3d(i,j,k,n)=tthe(kbox,1)+q*tthe(kbox,2)
                temp(i,j,k,n)=ttem(kbox,1)+q*ttem(kbox,2)
                saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
     &                               temp(i,j,k,n))
              endif
            else
c ---         all of a single layer
              temp(i,j,k,n)=ttem(kbox,1)
              saln(i,j,k,n)=tsal(kbox,1)
              th3d(i,j,k,n)=tthe(kbox,1)
              q = 0.0 !for debugging only
            endif !part:all of single layer
          else
c ---       multiple layers.
            if     (ktop.le.k .and. kbox.ge.k) then
              ka = k
            elseif (kbox-ktop.ge.3) then
              ka = (kbox+ktop)/2
            elseif (dprs(ktop).ge.dprs(kbox)) then
              ka = ktop
            else
              ka = kbox
            endif !choose ka
c ---       calculate as perturbation from layer ka (reduces roundoff)
            offset(1)=ttem(ka,1)
            offset(2)=tsal(ka,1)
            offset(3)=tthe(ka,1)
c
            qtop = pres(ktop+1)-ztop !partial layer thickness
            zcen = 0.5*(ztop+pres(ktop+1))
            if     (dprs(ktop).gt.dpthin) then
              q = 0.5 - (pres(ktop+1)-zcen)/dprs(ktop)
            else
              q = 0.0
            endif
            tsum =((ttem(ktop,1)+q*ttem(ktop,2))-offset(1))*qtop
            ssum =((tsal(ktop,1)+q*tsal(ktop,2))-offset(2))*qtop
            thsum=((tthe(ktop,1)+q*tthe(ktop,2))-offset(3))*qtop
c
            do k1= ktop+1,kbox-1
              tsum =tsum +(ttem(k1,1)-offset(1))*dprs(k1)
              ssum =ssum +(tsal(k1,1)-offset(2))*dprs(k1)
              thsum=thsum+(tthe(k1,1)-offset(3))*dprs(k1)
            enddo !k1
c
            qbot = zbox-pres(kbox) !partial layer thickness
            zcen = 0.5*(pres(kbox)+zbox)
            if     (dprs(kbox).gt.dpthin) then
              q = 0.5 - (pres(kbox+1)-zcen)/dprs(kbox)
            else
              q = 0.0
            endif
            tsum =tsum +((ttem(kbox,1)+q*ttem(kbox,2))-offset(1))*qbot
            ssum =ssum +((tsal(kbox,1)+q*tsal(kbox,2))-offset(2))*qbot
            thsum=thsum+((tthe(kbox,1)+q*tthe(kbox,2))-offset(3))*qbot
c
            rpsum=1.0d0/zthk
            if     (hybflg.eq.0) then  !T&S
              temp(i,j,k,n)=offset(1)+tsum*rpsum
              saln(i,j,k,n)=offset(2)+ssum*rpsum
              th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
            elseif (hybflg.eq.1) then  !th&S
              th3d(i,j,k,n)=offset(3)+thsum*rpsum
              saln(i,j,k,n)=offset(2)+ ssum*rpsum
              temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            elseif (hybflg.eq.2) then  !th&T
              th3d(i,j,k,n)=offset(3)+thsum*rpsum
              temp(i,j,k,n)=offset(1)+ tsum*rpsum
              saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
            endif
          endif !single or multiple layers
        else
c ---     thin or bottomed layer
            if (i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,i3,f7.3)')
     &          'thin k =',k,th3d(i,j,k-1,n)
            endif !debug
          if (k.le.nhybrd) then
c ---       fill with fluid from above
            th3d(i,j,k,n)=th3d(i,j,k-1,n)
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=temp(i,j,k-1,n)
          elseif (hybflg.ne.2) then
c ---       fill with saln from above
            th3d(i,j,k,n)=theta(i,j,k)
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          else
c ---       fill with temp from above
            th3d(i,j,k,n)=theta(i,j,k)
            temp(i,j,k,n)=temp(i,j,k-1,n)
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
        endif !normal:thin layer
c
        if     (lconserve) then  !usually .false.
          asum(1,1) = asum(1,1) + ttem(    k,1)*dprs(k)
          asum(1,2) = asum(1,2) + temp(i,j,k,n)*zthk
          asum(2,1) = asum(2,1) + tsal(    k,1)*dprs(k)
          asum(2,2) = asum(2,2) + saln(i,j,k,n)*zthk
          asum(3,1) = asum(3,1) + tthe(    k,1)*dprs(k)
          asum(3,2) = asum(3,2) + th3d(i,j,k,n)*zthk
        endif !lconserve
c
      enddo !k
c
      if (i.eq.itest .and. j.eq.jtest) then
        write (lp,'(i9,3a/(i9,3f23.17))')
     &  nstep,
     &  '                   dens',
     &  '                  thkns',
     &  '                   dpth',
     &  (k,tthe(k,1),    dprs(k)*qonem,    pres(k+1)*qonem,
     &   k,th3d(i,j,k,n),dp(i,j,k,n)*qonem,p(i,j,k+1)*qonem,
     &  k=1,kk)
        call flush(lp)
      endif !debug
c
      if     (lconserve) then  !usually .false.
c
c ---   inforce water column conservation
c
        do ktr=1,5
          q = asum(ktr,1)-asum(ktr,2)
          if     (q.eq.0.0) then
            offset(ktr) = 0.0
          elseif (abs(asum(ktr,2)).lt.2.0*abs(q)) then
            offset(ktr) = sign(zp5,q*asum(ktr,2))  !        -0.5 or  +0.5
          else
            offset(ktr) =          q/asum(ktr,2)   !between -0.5 and +0.5
          endif
        enddo !ktr
        do k=1,kk
          if     (hybflg.eq.0) then  !T&S
            temp(i,j,k,n)=temp(i,j,k,n)*(1.0+offset(1))
            saln(i,j,k,n)=saln(i,j,k,n)*(1.0+offset(2))
            th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
          elseif (hybflg.eq.1) then  !th&S
            saln(i,j,k,n)=saln(i,j,k,n)*(1.0+offset(2))
            th3d(i,j,k,n)=th3d(i,j,k,n)*(1.0+offset(3))
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,
     &                           saln(i,j,k,n))
          elseif (hybflg.eq.2) then  !th&T
            temp(i,j,k,n)=temp(i,j,k,n)*(1.0+offset(1))
            th3d(i,j,k,n)=th3d(i,j,k,n)*(1.0+offset(3))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
     &                           temp(i,j,k,n))
          endif
c
          if     (.false.) then !debugging
            zthk = dp(i,j,k,n)
            asum(1,3) = asum(1,3) + temp(i,j,k,n)*zthk
            asum(2,3) = asum(2,3) + saln(i,j,k,n)*zthk
            asum(3,3) = asum(3,3) + th3d(i,j,k,n)*zthk
          endif !debuging
        enddo !k
c
        if     (.false. .and. !debugging
     &          i.eq.itest .and. j.eq.jtest) then
          do ktr= 1,5
            write(lp,'(a,1p4e16.8,i3)')
     &        'hybgen,sum:',
     &        asum(ktr,1)/p(i,j,kk+1),
     &        asum(ktr,2)/p(i,j,kk+1),
     &        asum(ktr,3)/p(i,j,kk+1),
     &        offset(ktr),ktr
          enddo !ktr
        endif !debugging .and. i.eq.itest .and. j.eq.jtest
        if     (.false. .and. !debugging
     &          j.eq.jtest) then
          ktr=6
*         if     (abs(offset(ktr)).gt.1.e-08) then
          if     (abs(offset(ktr)).gt.1.e-12) then
            write(lp,'(a,1p4e16.8,i3)')
     &        'hybgen,sum:',
     &        asum(ktr,1)/p(i,j,kk+1),
     &        asum(ktr,2)/p(i,j,kk+1),
     &        asum(ktr,3)/p(i,j,kk+1),
     &        offset(ktr),i
          endif !large offset
        endif !debugging .and. j.eq.jtest
      endif !lconserve
c
      if (i.eq.itest .and. j.eq.jtest) then
        write (lp,'(i9,3a/(i9,3f23.17))')
     &  nstep,
     &  '                   dens',
     &  '                  thkns',
     &  '                   dpth',
     &  (k,tthe(k,1),    dprs(k)*qonem,    pres(k+1)*qonem,
     &   k,th3d(i,j,k,n),dp(i,j,k,n)*qonem,p(i,j,k+1)*qonem,
     &  k=1,kk)
        call flush(lp)
      endif !debug
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f9.3,f9.2))
      if (i.eq.itest .and. j.eq.jtest) then
        write (lp,103) nstep,itest+i0,jtest+j0,
     &  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
     &  (k,ttem(k,1),tsal(k,1),tthe(k,1)+thbase,
     &   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
     &   k,temp(i,j,k,n),saln(i,j,k,n),
     &   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
     &   p(i,j,k+1)*qonem,
     &  k=1,kk),
     &  (k,ttem(k,1),tsal(k,1),tthe(k,1)+thbase,
     &   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
     &   k=kk+1,kp)
        call flush(lp)
      endif !debug
c
 2    continue  !i
c
c --- to avoid roundoff errors in -dpudpv- after restart, make sure p=p(dp)
c
      do 1 k=1,kk
*     do 1 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
 1    continue
c
      return
      end
