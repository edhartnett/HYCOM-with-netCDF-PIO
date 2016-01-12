      subroutine latbdp(n)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      integer n
c
c --- apply lateral boundary conditions to   barotropic  flow field
c
c --- port flow version:
c --- similar to the standard 'Browning and Kreiss' MICOM/HYCOM open
c --- boundary condition, except that the exterior normal velocity
c --- is constant in time and exterior pressure = interior pressure.
c --- tangential velocity is not constrained.
c
c --- the code is as similar as possible to that for the standard case.
c --- so for example, 'speed' is in fact 1/SQRT(gH) which represents
c --- c1/g in the notation of (Bleck and Sun, Open boundary conditions
c --- for MICOM).  The 1/g allows for the use of pressure fields.
c
c --- the first call is made during initialization.
c
c --- Alan J. Wallcraft,  NRL,  November 1999.
c
      logical    ldebug_latbdp
      parameter (ldebug_latbdp=.false.)
c
      integer    nchar
      parameter (nchar=120)
c
      logical     lfatal,lfatalp
      integer     i,j,isec,ifrst,ilast,l
      real        aline(nchar),
     &            dline(itdm+jtdm),xline(itdm+jtdm),
     &            pline(itdm+jtdm),uline(itdm+jtdm,2)
      real        crs,fin,sum,svspin,uvscl,uvscl2
      real*8      tstep
      character*3 char3
c
      integer nports,kdport(9),
     &        ifport(9),ilport(9),jfport(9),jlport(9),lnport(9)
      real    pefold,svpnow(9),svport(9)
      real*8  refold
      save    nports,kdport,ifport,ilport,jfport,jlport,lnport
      save    pefold,svpnow,svport,refold
c
      real    uportw(jtdm),speedw(jtdm),rspedw(jtdm),
     &        uporte(jtdm),speede(jtdm),rspede(jtdm),
     &        vportn(itdm),speedn(itdm),rspedn(itdm),
     &        vports(itdm),speeds(itdm),rspeds(itdm)
      save    uportw,speedw,rspedw,uporte,speede,rspede,
     &        vportn,speedn,rspedn,vports,speeds,rspeds
c
      character*13 fmt
      save         fmt
      data         fmt / '(i4,1x,120i1)' /
c
      integer lcount
      save    lcount
      data    lcount / 0 /
c
      lcount = lcount + 1
c
c --- the first call just initializes data structures.
c
      if     (lcount.eq.1) then
c
        open(unit=99,file='ports.input')
c
c ---   'nports' = number of boundary port sections.
        call blkini(nports,'nports')
        if     (mnproc.eq.1) then
        write(lp,*)
        endif
        if     (nports.lt.0 .or. nports.gt.9) then
          if     (mnproc.eq.1) then
          write(lp,*) 
          write(lp,*) 'error in latbdp - illegal nports value'
          write(lp,*) 
          call flush(lp)
          endif
          call xcstop('(latbdp)')
                 stop '(latbdp)'
        endif
c
c ---   'pefold' = port transport e-folding time in days
        call blkinr(pefold,'pefold','(a6," =",f10.4," days")')
        if     (mnproc.eq.1) then
        write(lp,*)
        endif
c
c ---   switch units from days to baroclinic time steps
c ---   shift lcount to prevent underflow (lcount*refold.ge.0.001)
c
        tstep  = pefold*(86400.d0/batrop)
        refold = 1.d0/tstep
        lcount = lcount + int(tstep)/1000
c
c ---   read in the ports one at a time
c
        do l= 1,nports
c
c ---     port location is w.r.t. u (EW) or v (NS) grid
c ---     and identifies the sea at the port
c ---     the minimum index is 0
c
c ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)
c ---     'ifport' = first i-index
c ---     'ilport' = last  i-index (=ifport for N or S orientation)
c ---     'jfport' = first j-index
c ---     'jlport' = last  j-index (=jfport for E or W orientation)
c ---     'svpnow' = existing port transport in Sv (+ve towards E or S)
c ---     'svport' = target   port transport in Sv (+ve towards E or S)
c ---     'lnport' = port length (calculated, not input)
          call blkini(kdport(l),'kdport')
          call blkini(ifport(l),'ifport')
          call blkini(ilport(l),'ilport')
          call blkini(jfport(l),'jfport')
          call blkini(jlport(l),'jlport')
          call blkinr(svpnow(l),'svpnow','(a6," =",f10.4," Sv")')
          call blkinr(svport(l),'svport','(a6," =",f10.4," Sv")')
          if     (mnproc.eq.1) then
          write(lp,*)
          endif
c
          lnport(l) = ilport(l)-ifport(l)+jlport(l)-jfport(l)+1
c
c ---     sanity check.
c
          if     (kdport(l).gt.2) then
            if     (ifport(l).ne.ilport(l)) then
              if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error in latbdp - port direction',
     &                     ' and orientation are not consistent'
              write(lp,*) 
              call flush(lp)
              endif
              call xcstop('(latbdp)')
                     stop '(latbdp)'
            endif
          else
            if     (jfport(l).ne.jlport(l)) then
              if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error in latbdp - port direction',
     &                     ' and orientation are not consistent'
              write(lp,*) 
              call flush(lp)
              endif
              call xcstop('(latbdp)')
                     stop '(latbdp)'
            endif
          endif
          if     (ifport(l).gt.ilport(l) .or.
     &            jfport(l).gt.jlport(l)     ) then
            if     (mnproc.eq.1) then
            write(lp,*) 
            write(lp,*) 'error in latbdp - port',
     &                   ' location is not consistent'
            write(lp,*) 
            call flush(lp)
            endif
            call xcstop('(latbdp)')
                   stop '(latbdp)'
          endif
        enddo
c
        close(unit=99)
c
c ---   check ports against masks,
c ---   mark the port locations on masks and print them out.
c
        lfatal = .false.
        do l= 1,nports
          lfatalp = .false.
c
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.1 .or. i.gt.itdm-2 .or.
     &                j.lt.1 .or. j.gt.jtdm       ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0+1,j-j0).ne.1 .or.
     &                iu(i-i0+2,j-j0).ne.1     ) then
                lfatalp = .true.
              endif
            enddo
c
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.3 .or. i.gt.itdm .or.
     &                j.lt.1 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0-1,j-j0).ne.1 .or.
     &                iu(i-i0-2,j-j0).ne.1     ) then
                lfatalp = .true.
              endif
            enddo
c
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm .or.
     &                j.lt.3 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0-1).ne.1 .or.
     &                iv(i-i0,j-j0-2).ne.1     ) then
                lfatalp = .true.
              endif
            enddo
c
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm   .or.
     &                j.lt.1 .or. j.gt.jtdm-2     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0+1).ne.1 .or.
     &                iv(i-i0,j-j0+2).ne.1     ) then
                lfatalp = .true.
              endif
            enddo
c
          endif
c
          if     (lfatalp) then
            write(lp,*) 
            write(lp,*) 'error in latbdp - port ',l,' mislocated',
     &                  '  (mnproc = ',mnproc,')'
            write(lp,*) 
            call flush(lp)
          endif
          lfatal = lfatal .or. lfatalp
        enddo
c
c ---   write out  -iu-  and -iv- arrays, if they are not too big
c ---   data are written in strips nchar points wide
        if     (lfatal .or. max(itdm,jtdm).le.2*nchar) then
          util1(1:ii,1:jj) = iu(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'iu array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
c
          util1(1:ii,1:jj) = iv(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'iv array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
        endif  ! small region
c
        if     (lfatal) then
          write(lp,*) 
          write(lp,*) 'error in latbdp - bad port(s)'
          write(lp,*) 
          call flush(lp)
          call xchalt('(latbdp)')
                 stop '(latbdp)'
        endif
c
c ---   restore iu and iv, and zero iuopn and ivopn.
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j= 1,jj
          do i= 1,ii
            iu(i,j) = max( iu(i,j), 0 )
            iv(i,j) = max( iv(i,j), 0 )
          enddo
        enddo
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            iuopn(i,j) = 0
            ivopn(i,j) = 0
          enddo
        enddo
c
c ---   initialize the ports
c
        do l= 1,nports
          if     (kdport(l).eq.4) then
c
c           western port
c
            sum = 0.0
            i = ifport(l)
            j = jfport(l)
            call xclget(dline(j),lnport(l), depths,i+1,j,0,1, 0)
            call xclget(xline(j),lnport(l), scuy,  i,  j,0,1, 0)
            do j= jfport(l),jlport(l)
              sum = sum + dline(j)*xline(j)
            enddo
            sum = 1.e6/sum
            do j= jfport(l),jlport(l)
              uportw(j) = sum
              speedw(j) = sqrt(1.0/(onem*dline(j)))
              rspedw(j) = 1.0/speedw(j)
              if     (mnproc.eq.1) then
              write(lp,'(a,i2,2i4,1p2e12.5)') 
     &          'w port: ',l,i,j,uportw(j),speedw(j)
              endif
c
              if     (i.ge.i0+ 1-nbdy .and.
     &                i.le.i0+ii+nbdy .and.
     &                j.ge.j0+ 1-nbdy .and.
     &                j.le.j0+jj+nbdy      ) then
                iuopn(i-i0,j-j0) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            sum = 0.0
            i = ifport(l)-1
            j = jfport(l)
            call xclget(dline(j),lnport(l), depths,i,  j,0,1, 0)
            call xclget(xline(j),lnport(l), scuy,  i+1,j,0,1, 0)
            do j= jfport(l),jlport(l)
              sum = sum + dline(j)*xline(j)
            enddo
            sum = 1.e6/sum
            do j= jfport(l),jlport(l)
              uporte(j) = sum
              speede(j) = sqrt(1.0/(onem*dline(j)))
              rspede(j) = 1.0/speede(j)
              if     (mnproc.eq.1) then
              write(lp,'(a,i2,2i4,1p2e12.5)') 
     &          'e port: ',l,i,j,uporte(j),speede(j)
              endif
c
              if     (i+1.ge.i0+ 1-nbdy .and.
     &                i+1.le.i0+ii+nbdy .and.
     &                j  .ge.j0+ 1-nbdy .and.
     &                j  .le.j0+jj+nbdy      ) then
                iuopn(i-i0+1,j-j0) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            sum = 0.0
            j = jfport(l)-1
            i = ifport(l)
            call xclget(dline(i),lnport(l), depths,i,j,  1,0, 0)
            call xclget(xline(i),lnport(l), scuy,  i,j+1,1,0, 0)
            do i= ifport(l),ilport(l)
              sum = sum + dline(i)*xline(i)
            enddo
            sum = 1.e6/sum
            do i= ifport(l),ilport(l)
              vportn(i) = sum
              speedn(i) = sqrt(1.0/(onem*dline(i)))
              rspedn(i) = 1.0/speedn(i)
              if     (mnproc.eq.1) then
              write(lp,'(a,i2,2i4,1p2e12.5)') 
     &          'n port: ',l,i,j,vportn(i),speedn(i)
              endif
c
              if     (i  .ge.i0+ 1-nbdy .and.
     &                i  .le.i0+ii+nbdy .and.
     &                j+1.ge.j0+ 1-nbdy .and.
     &                j+1.le.j0+jj+nbdy      ) then
                ivopn(i-i0,j-j0+1) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            sum = 0.0
            j = jfport(l)
            i = ifport(l)
            call xclget(dline(i),lnport(l), depths,i,j+1,1,0, 0)
            call xclget(xline(i),lnport(l), scuy,  i,j,  1,0, 0)
            do i= ifport(l),ilport(l)
              sum = sum + dline(i)*xline(i)
            enddo
            sum = 1.e6/sum
            do i= ifport(l),ilport(l)
              vports(i) = sum
              speeds(i) = sqrt(1.0/(onem*dline(i)))
              rspeds(i) = 1.0/speeds(i)
              if     (mnproc.eq.1) then
              write(lp,'(a,i2,2i4,1p2e12.5)') 
     &          's port: ',l,i,j,vports(i),speeds(i)
              endif
c
              if     (i.ge.i0+ 1-nbdy .and.
     &                i.le.i0+ii+nbdy .and.
     &                j.ge.j0+ 1-nbdy .and.
     &                j.le.j0+jj+nbdy      ) then
                ivopn(i-i0,j-j0) = 1
              endif
            enddo
c
          endif
c
          if     (mnproc.eq.1) then
          write(lp,*) 'port, now/target velocity = ',
     &                l,svpnow(l)*sum,svport(l)*sum
          call flush(lp)
          endif
        enddo
        if     (mnproc.eq.1) then
        write(lp,*) 
        call flush(lp)
        endif
c
c       end of initialization
c
        call xcsync(flush_lp)
        return
      endif
c
c --- 'wellposed' treatment of pressure and normal velocity fields
c --- not in fact wellposed with this exterior data
c
      tstep  = lcount
      svspin = exp( -tstep*refold )
      do l= 1,nports
        uvscl = svport(l) + svspin*(svpnow(l)-svport(l))
c
        if     (kdport(l).eq.4) then
c
c         western port
c
          i = ifport(l)
          j = jfport(l)
          call xclget(dline(j),  lnport(l),
     &                depthu,                 i+1,j,0,1, 0)
          call xclget(xline(j),  lnport(l),
     &                scuy,                   i,  j,0,1, 0)
          call xclget(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,0,1, 0)
          call xclget(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,0,1, 0)
          call xclget(uline(j,2),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+2,j,0,1, 0)
          sum = 0.0
          do j= jfport(l),jlport(l)
            crs=uvscl*uportw(j)+speedw(j)*pline(j)
            fin=1.5*uline(j,1)-.5*uline(j,2)-speedw(j)*pline(j)
            sum=sum+((crs+fin)-uline(j,1))*dline(j)*xline(j)
          enddo
          uvscl2 = uvscl + (uvscl - sum/(onem*1.e6))
          sum = 0.0
          do j= jfport(l),jlport(l)
            crs=uvscl2*uportw(j)+speedw(j)*pline(j)
            fin=1.5*uline(j,1)-.5*uline(j,2)-speedw(j)*pline(j)
            pline(j)  =.5*(crs-fin)*rspedw(j)
            uline(j,1)=(crs+fin)-uline(j,1)
            sum=sum+uline(j,1)*dline(j)*xline(j)
          enddo
          j = jfport(l)
          call xclput(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,0,1)
          call xclput(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i,  j,0,1)
c
              if     (ldebug_latbdp .and. mnproc.eq.1) then
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(pb - ',
     &                                      l,lnport(l),i,  j,0,1
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(ub - ',
     &                                      l,lnport(l),i,  j,0,1
                call flush(lp)
              endif
c
        elseif (kdport(l).eq.3) then
c
c         eastern port
c
          i = ifport(l)-1
          j = jfport(l)
          call xclget(dline(j),  lnport(l),
     &                depthu,                 i+1,j,0,1, 0)
          call xclget(xline(j),  lnport(l),
     &                scuy,                   i+1,j,0,1, 0)
          call xclget(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,0,1, 0)
          call xclget(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i,  j,0,1, 0)
          call xclget(uline(j,2),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i-1,j,0,1, 0)
          sum = 0.0
          do j= jfport(l),jlport(l)
            crs=uvscl*uporte(j)-speede(j)*pline(j)
            fin=1.5*uline(j,1)-.5*uline(j,2)+speede(j)*pline(j)
            sum=sum+((crs+fin)-uline(j,1))*dline(j)*xline(j)
          enddo
          uvscl2 = uvscl + (uvscl - sum/(onem*1.e6))
          sum = 0.0
          do j= jfport(l),jlport(l)
            crs=uvscl2*uporte(j)-speede(j)*pline(j)
            fin=1.5*uline(j,1)-.5*uline(j,2)+speede(j)*pline(j)
            pline(j)  =.5*(fin-crs)*rspede(j)
            uline(j,1)=(fin+crs)-uline(j,1)
            sum=sum+uline(j,1)*dline(j)*xline(j)
*             if     (mnproc.eq.1) then
*             write(lp,'(a,i2,2i4,1p2e12.5)') 
*    &          'e port: ',l,i,j,pline(j),uline(j,1)
*             endif
          enddo
          j = jfport(l)
          call xclput(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,0,1)
          call xclput(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,0,1)
c
              if     (ldebug_latbdp .and. mnproc.eq.1) then
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(pb - ',
     &                                      l,lnport(l),i,  j,0,1
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(ub - ',
     &                                      l,lnport(l),i+1,j,0,1
                call flush(lp)
              endif
c
        elseif (kdport(l).eq.1) then
c
c         northern port
c
          j = jfport(l)-1
          i = ifport(l)
          call xclget(dline(i),  lnport(l),
     &                depthv,                 i,j+1,1,0, 0)
          call xclget(xline(i),  lnport(l),
     &                scux,                   i,j+1,1,0, 0)
          call xclget(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,j,  1,0, 0)
          call xclget(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j,  1,0, 0)
          call xclget(uline(i,2),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j-1,1,0, 0)
          sum = 0.0
          do i= ifport(l),ilport(l)
            crs=uvscl*vportn(i)-speedn(i)*pline(i)
            fin=1.5*uline(i,1)-.5*uline(i,2)+speedn(i)*pline(i)
            sum=sum+((fin+crs)-uline(i,1))*dline(i)*xline(i)
          enddo
          uvscl2 = uvscl + (uvscl - sum/(onem*1.e6))
          sum = 0.0
          do i= ifport(l),ilport(l)
            crs=uvscl2*vportn(i)-speedn(i)*pline(i)
            fin=1.5*uline(i,1)-.5*uline(i,2)+speedn(i)*pline(i)
            pline(i)  =.5*(fin-crs)*rspedn(i)
            uline(i,1)=(fin+crs)-uline(i,1)
            sum=sum+uline(i,1)*dline(i)*xline(i)
          enddo
          i = ifport(l)
          call xclput(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,j,  1,0)
          call xclput(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j+1,1,0)
c
              if     (ldebug_latbdp .and. mnproc.eq.1) then
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(pb - ',
     &                                      l,lnport(l),i,j,  1,0
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(vb - ',
     &                                      l,lnport(l),i,j+1,1,0
                call flush(lp)
              endif
c
        elseif (kdport(l).eq.2) then
c
c         southern port
c
          j = jfport(l)
          i = ifport(l)
          call xclget(dline(i),  lnport(l),
     &                depthv,                 i,j,  1,0, 0)
          call xclget(xline(i),  lnport(l),
     &                scux,                   i,j,  1,0, 0)
          call xclget(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,j,  1,0, 0)
          call xclget(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j+1,1,0, 0)
          call xclget(uline(i,2),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j+2,1,0, 0)
          sum = 0.0
          do i= ifport(l),ilport(l)
            crs=uvscl*vports(i)+speeds(i)*pline(i)
            fin=1.5*uline(i,1)-.5*uline(i,2)-speeds(i)*pline(i)
            sum=sum+((crs+fin)-uline(i,1))*dline(i)*xline(i)
          enddo
          uvscl2 = uvscl + (uvscl - sum/(onem*1.e6))
          sum = 0.0
          do i= ifport(l),ilport(l)
            crs=uvscl2*vports(i)+speeds(i)*pline(i)
            fin=1.5*uline(i,1)-.5*uline(i,2)-speeds(i)*pline(i)
            pline(i)  =.5*(crs-fin)*rspeds(i)
            uline(i,1)=(crs+fin)-uline(i,1)
            sum=sum+uline(i,1)*dline(i)*xline(i)
          enddo
          i = ifport(l)
          call xclput(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,j,  1,0)
          call xclput(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,j,  1,0)
c
              if     (ldebug_latbdp .and. mnproc.eq.1) then
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(pb - ',
     &                                      l,lnport(l),i,j,  1,0
                write(lp,'(a,i2,3i5,2i2)') 'l,xclput(vb - ',
     &                                      l,lnport(l),i,j,  1,0
                call flush(lp)
              endif
c
        endif
c
*       if     (mod(lcount,512).eq.0) then
*         if     (mnproc.eq.1) then
*         write(lp,*) 'latbdp - l,sv,sum = ',l,uvscl,sum/(onem*1.e6)
*         call flush(lp)
*         endif
*       endif
      enddo
c
      return
      end
      subroutine latbdt(n,lll)
      use mod_xc  ! HYCOM communication interface
      implicit none
      include 'common_blocks.h'
c
      integer n,lll
c
c --- apply lateral boundary conditions to   barotropic  flow field
c
c --- nested sub-region version:
c --- Uses the 'Browning and Kreiss' MICOM open boundary condition.
c
c --- Note that 'speed' is in fact thref*SQRT(gH) which represents
c --- c1*thref/g in the notation of (Bleck and Sun, Open boundary conditions
c --- for MICOM).  The thref/g allows for the use of pressure fields.
c
c --- the first call is made during initialization.
c
c --- Alan J. Wallcraft,  NRL,  July, 2001.
c
      logical    ldebug_latbdt
      parameter (ldebug_latbdt=.false.)
c
      integer    nchar
      parameter (nchar=120)
c
      logical     lfatal,lfatalp
      integer     i,j,isec,ifrst,ilast,l,npf,npi,npl
      real        aline(nchar),
     &            pline(itdm+jtdm),uline(itdm+jtdm,2)
      real        crs,fin
      character*3 char3
c
      integer nports,kdport(9),
     &        ifport(9),ilport(9),jfport(9),jlport(9),lnport(9)
      save    nports,kdport,ifport,ilport,jfport,jlport,lnport
c
      real    speedw(jtdm),rspedw(jtdm),
     &        speede(jtdm),rspede(jtdm),
     &        speedn(itdm),rspedn(itdm),
     &        speeds(itdm),rspeds(itdm)
      real    plnstw(jtdm),ulnstw(jtdm),vlnstw(jtdm),
     &        plnste(jtdm),ulnste(jtdm),vlnste(jtdm),
     &        plnstn(itdm),ulnstn(itdm),vlnstn(itdm),
     &        plnsts(itdm),ulnsts(itdm),vlnsts(itdm)
      save    speedw,rspedw,plnstw,ulnstw,vlnstw,
     &        speede,rspede,plnste,ulnste,vlnste,
     &        speedn,rspedn,plnstn,ulnstn,vlnstn,
     &        speeds,rspeds,plnsts,ulnsts,vlnsts
c
      character*13 fmt
      save         fmt
      data         fmt / '(i4,1x,120i1)' /
c
      integer lcount
      save    lcount
      data    lcount / 0 /
c
      lcount = lcount + 1
c
c --- the first call just initializes data structures.
c
      if     (lcount.eq.1) then
c
        open(unit=99,file='ports.input')
c
c ---   'nports' = number of boundary port sections.
        call blkini(nports,'nports')
        if     (mnproc.eq.1) then
        write(lp,*)
        endif
        if     (nports.lt.0 .or. nports.gt.9) then
          if     (mnproc.eq.1) then
          write(lp,*) 
          write(lp,*) 'error in latbdt - illegal nports value'
          write(lp,*) 
          call flush(lp)
          endif
          call xcstop('(latbdt)')
                 stop '(latbdt)'
        endif
c
c ---   read in the ports one at a time
c
        do l= 1,nports
c
c ---     port location is w.r.t. u (EW) or v (NS) grid
c ---     and identifies the sea at the port
c ---     the minimum index is 0
c
c ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)
c ---     'ifport' = first i-index
c ---     'ilport' = last  i-index (=ifport for N or S orientation)
c ---     'jfport' = first j-index
c ---     'jlport' = last  j-index (=jfport for E or W orientation)
c ---     'lnport' = port length (calculated, not input)
          call blkini(kdport(l),'kdport')
          call blkini(ifport(l),'ifport')
          call blkini(ilport(l),'ilport')
          call blkini(jfport(l),'jfport')
          call blkini(jlport(l),'jlport')
          if     (mnproc.eq.1) then
          write(lp,*)
          endif
c
          lnport(l) = ilport(l)-ifport(l)+jlport(l)-jfport(l)+1
c
c ---     sanity check.
c
          if     (kdport(l).gt.2) then
            if     (ifport(l).ne.ilport(l)) then
              if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error in latbdt - port direction',
     &                     ' and orientation are not consistent'
              write(lp,*) 
              call flush(lp)
              endif
              call xcstop('(latbdt)')
                     stop '(latbdt)'
            endif
          else
            if     (jfport(l).ne.jlport(l)) then
              if     (mnproc.eq.1) then
              write(lp,*) 
              write(lp,*) 'error in latbdt - port direction',
     &                     ' and orientation are not consistent'
              write(lp,*) 
              call flush(lp)
              endif
              call xcstop('(latbdt)')
                     stop '(latbdt)'
            endif
          endif
          if     (ifport(l).gt.ilport(l) .or.
     &            jfport(l).gt.jlport(l)     ) then
            if     (mnproc.eq.1) then
            write(lp,*) 
            write(lp,*) 'error in latbdt - port',
     &                   ' location is not consistent'
            write(lp,*) 
            call flush(lp)
            endif
            call xcstop('(latbdt)')
                   stop '(latbdt)'
          endif
        enddo
c
        close(unit=99)
c
c ---   check ports against masks,
c ---   mark the port locations on masks and print them out.
c
        lfatal = .false.
        do l= 1,nports
          lfatalp = .false.
c
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.1 .or. i.gt.itdm-2 .or.
     &                j.lt.1 .or. j.gt.jtdm       ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0+1,j-j0).ne.1 .or.
     &                iu(i-i0+2,j-j0).ne.1     ) then
                lfatalp = .true.
              endif
            enddo
c
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.3 .or. i.gt.itdm .or.
     &                j.lt.1 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iu(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iu(i-i0,j-j0) =  9
              else
                iu(i-i0,j-j0) = -1
              endif
              if     (iu(i-i0-1,j-j0).ne.1 .or.
     &                iu(i-i0-2,j-j0).ne.1     ) then
                lfatalp = .true.
              endif
            enddo
c
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm .or.
     &                j.lt.3 .or. j.gt.jtdm     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0-1).ne.1 .or.
     &                iv(i-i0,j-j0-2).ne.1     ) then
                lfatalp = .true.
              endif
            enddo
c
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.itdm   .or.
     &                j.lt.1 .or. j.gt.jtdm-2     ) then
                lfatalp = .true.
              elseif (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              elseif (iv(i-i0,j-j0).ne.0) then
                lfatalp = .true.
                iv(i-i0,j-j0) =  9
              else
                iv(i-i0,j-j0) = -1
              endif
              if     (iv(i-i0,j-j0+1).ne.1 .or.
     &                iv(i-i0,j-j0+2).ne.1     ) then
                lfatalp = .true.
              endif
            enddo
c
          endif
c
          if     (lfatalp) then
            write(lp,*) 
            write(lp,*) 'error in latbdt - port ',l,' mislocated',
     &                  '  (mnproc = ',mnproc,')'
            write(lp,*) 
            call flush(lp)
          endif
          lfatal = lfatal .or. lfatalp
        enddo  !l=1,nports
c
c ---   write out  -iu-  and -iv- arrays, if they are not too big
c ---   data are written in strips nchar points wide
c
        if     (lfatal .or. max(itdm,jtdm).le.2*nchar) then
          util1(1:ii,1:jj) = iu(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'iu array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
c
          util1(1:ii,1:jj) = iv(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'iv array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
        endif  ! small region
c
        if     (lfatal) then
          write(lp,*) 
          write(lp,*) 'error in latbdt - bad port(s)'
          write(lp,*) 
          call flush(lp)
          call xchalt('(latbdt)')
                 stop '(latbdt)'
        endif
c
c ---   restore iu and iv, and zero iuopn and ivopn.
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j= 1,jj
          do i= 1,ii
            iu(i,j) = max( iu(i,j), 0 )
            iv(i,j) = max( iv(i,j), 0 )
          enddo
        enddo
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            iuopn(i,j) = 0
            ivopn(i,j) = 0
          enddo
        enddo
c
c ---   define the nested boundary input mask.
c
        do j= 1,jj
          do i= 1,ii
            maskbc(i,j) = 0
          enddo
        enddo
c
        do l= 1,nports
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              else
                maskbc(i-i0,j-j0) = 1
              endif
            enddo
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)-1
            do j= jfport(l),jlport(l)
              if     (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              else
                maskbc(i-i0,j-j0) = 1
              endif
            enddo
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)-1
            do i= ifport(l),ilport(l)
              if     (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              else
                maskbc(i-i0,j-j0) = 1
              endif
            enddo
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.le.i0 .or. i.gt.i0+ii .or.
     &                j.le.j0 .or. j.gt.j0+jj     ) then
                cycle  ! not on this tile.
              else
                maskbc(i-i0,j-j0) = 1
              endif
            enddo
          endif
        enddo  !l=1,nports
c
        if     (ldebug_latbdt) then
          util1(1:ii,1:jj) = maskbc(1:ii,1:jj)  ! xclget is for real arrays
          isec=(itdm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(itdm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            if     (mnproc.eq.1) then
            write (lp,'(a,i5,a,i5)')
     &        'bc array, cols',ifrst+1,' --',ilast
            endif
            do j= jtdm,1,-1
              call xclget(aline,ilast-ifrst, util1,ifrst+1,j,1,0, 1)
              if     (mnproc.eq.1) then
              write (lp,fmt) j,(nint(aline(i)),i=1,ilast-ifrst)
              endif
            enddo
          enddo
          if     (mnproc.eq.1) then
          write (lp,*)
          endif
          call xcsync(flush_lp)
        endif !ldebug_latbdt
c
c ---   initialize the ports
c
        do l= 1,nports
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            j = jfport(l)
            call xclget(pline(j),lnport(l), depths,i,j,0,1, 0)
            do j= jfport(l),jlport(l)
              speedw(j) = sqrt(thref/(onem*pline(j)))
              rspedw(j) = 1.0/speedw(j)
              if     (ldebug_latbdt .and. mnproc.eq.1) then
              write(lp,'(a,i2,2i4,1pe12.5)') 
     &          'w port: ',l,i,j,speedw(j)
              endif
c
              if     (i.ge.i0+ 1-nbdy .and.
     &                i.le.i0+ii+nbdy .and.
     &                j.ge.j0+ 1-nbdy .and.
     &                j.le.j0+jj+nbdy      ) then
                iuopn(i-i0,j-j0) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)-1
            j = jfport(l)
            call xclget(pline(j),lnport(l), depths,i,j,0,1, 0)
            do j= jfport(l),jlport(l)
              speede(j) = sqrt(thref/(onem*pline(j)))
              rspede(j) = 1.0/speede(j)
              if     (ldebug_latbdt .and. mnproc.eq.1) then
              write(lp,'(a,i2,2i4,1pe12.5)') 
     &          'e port: ',l,i,j,speede(j)
              endif
c
              if     (i+1.ge.i0+ 1-nbdy .and.
     &                i+1.le.i0+ii+nbdy .and.
     &                j  .ge.j0+ 1-nbdy .and.
     &                j  .le.j0+jj+nbdy      ) then
                iuopn(i-i0+1,j-j0) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)-1
            i = ifport(l)
            call xclget(pline(i),lnport(l), depths,i,j,1,0, 0)
            do i= ifport(l),ilport(l)
              speedn(i) = sqrt(thref/(onem*pline(i)))
              rspedn(i) = 1.0/speedn(i)
              if     (ldebug_latbdt .and. mnproc.eq.1) then
              write(lp,'(a,i2,2i4,1pe12.5)') 
     &          'n port: ',l,i,j,speedn(i)
              endif
c
              if     (i  .ge.i0+ 1-nbdy .and.
     &                i  .le.i0+ii+nbdy .and.
     &                j+1.ge.j0+ 1-nbdy .and.
     &                j+1.le.j0+jj+nbdy      ) then
                ivopn(i-i0,j-j0+1) = 1
              endif
            enddo
c
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            i = ifport(l)
            call xclget(pline(i),lnport(l), depths,i,j,1,0, 0)
            do i= ifport(l),ilport(l)
              speeds(i) = sqrt(thref/(onem*pline(i)))
              rspeds(i) = 1.0/speeds(i)
              if     (ldebug_latbdt .and. mnproc.eq.1) then
              write(lp,'(a,i2,2i4,1pe12.5)') 
     &          's port: ',l,i,j,speeds(i)
              endif
c
              if     (i.ge.i0+ 1-nbdy .and.
     &                i.le.i0+ii+nbdy .and.
     &                j.ge.j0+ 1-nbdy .and.
     &                j.le.j0+jj+nbdy      ) then
                ivopn(i-i0,j-j0) = 1
              endif
            enddo
c
          endif
        enddo  !l=1,nports
        if     (ldebug_latbdt .and. mnproc.eq.1) then
        write(lp,*) 
        call flush(lp)
        endif
c
c       end of initialization
c
        call xcsync(flush_lp)
        return
      endif
c
c --- nested input only required on first barotropic time step.
c
      if     (lll.eq.1) then
        do j= 1,jj
          do i= 1,ii
            if     (maskbc(i,j).eq.1) then
              util1(i,j) = ubnest(i,j,ln0)*wb0+ubnest(i,j,ln1)*wb1
              util2(i,j) = vbnest(i,j,ln0)*wb0+vbnest(i,j,ln1)*wb1
              util3(i,j) = pbnest(i,j,ln0)*wb0+pbnest(i,j,ln1)*wb1
              util4(i,j) = ubpnst(i,j,ln0)*wb0+ubpnst(i,j,ln1)*wb1
              util5(i,j) = vbpnst(i,j,ln0)*wb0+vbpnst(i,j,ln1)*wb1
            endif
          enddo
        enddo
c
        do l= 1,nports
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            j = jfport(l)
            call xclget(plnstw(j),  lnport(l),
     &                  util3(1-nbdy,1-nbdy),   i,  j,  0,1, 0)  ! pbnest
            call xclget(ulnstw(j),  lnport(l),
     &                  util4(1-nbdy,1-nbdy),   i,  j,  0,1, 0)  ! ubpnst
            call xclget(vlnstw(j+1),lnport(l)-1,
     &                  util2(1-nbdy,1-nbdy),   i,  j+1,0,1, 0)  ! vbnest
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)-1
            j = jfport(l)
            call xclget(plnste(j),  lnport(l),
     &                  util3(1-nbdy,1-nbdy),   i,  j,  0,1, 0)  ! pbnest
            call xclget(ulnste(j),  lnport(l),
     &                  util4(1-nbdy,1-nbdy),   i,  j,  0,1, 0)  ! ubpnst
            call xclget(vlnste(j+1),lnport(l)-1,
     &                  util2(1-nbdy,1-nbdy),   i,  j+1,0,1, 0)  ! vbnest
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)-1
            i = ifport(l)
            call xclget(plnstn(i),  lnport(l),
     &                  util3(1-nbdy,1-nbdy),   i,  j,  1,0, 0)  ! pbnest
            call xclget(vlnstn(i),  lnport(l),
     &                  util5(1-nbdy,1-nbdy),   i,  j,  1,0, 0)  ! vbpnst
            call xclget(ulnstn(i+1),lnport(l)-1,
     &                  util1(1-nbdy,1-nbdy),   i+1,j,  1,0, 0)  ! ubnest
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            i = ifport(l)
            call xclget(plnsts(i),  lnport(l),
     &                  util3(1-nbdy,1-nbdy),   i,  j,  1,0, 0)  ! pbnest
            call xclget(vlnsts(i),  lnport(l),
     &                  util5(1-nbdy,1-nbdy),   i,  j,  1,0, 0)  ! vbpnst
            call xclget(ulnsts(i+1),lnport(l)-1,
     &                  util1(1-nbdy,1-nbdy),   i+1,j,  1,0, 0)  ! ubnest
          endif
        enddo  !l=1,nports
      endif !lll.eq.1
c
c --- 'wellposed' treatment of pressure and velocity fields.
c --- alternate order of ports in case corners are open.
c
      if     (mod(lll,2).eq.1) then
        npf =  1
        npl =  nports
        npi =  1
      else
        npf =  nports
        npl =  1
        npi = -1
      endif
      do l= npf,npl,npi
c
        if     (kdport(l).eq.4) then
c
c         western port
c
          i = ifport(l)
          j = jfport(l)
          call xclget(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,  0,1, 0)
          call xclget(uline(j,2),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+2,j,  0,1, 0)
          call xclget(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  0,1, 0)
          do j= jfport(l),jlport(l)
            crs=                   ulnstw(j) +speedw(j)*plnstw(j)
            fin=1.5*uline(j,1)-0.5*uline(j,2)-speedw(j)*pline( j)
            pline(j)  =0.5*(crs-fin)*rspedw(j)
            uline(j,1)=    (crs+fin)-uline(j,1)
          enddo
          j = jfport(l)
          call xclput(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  0,1)
          call xclput(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i,  j,  0,1)  ! normal
          call xclput(vlnstw(j+1),lnport(l)-1,
     &                vbavg(1-nbdy,1-nbdy,n), i,  j+1,0,1)  ! tangential
c
        elseif (kdport(l).eq.3) then
c
c         eastern port
c
          i = ifport(l)-1
          j = jfport(l)
          call xclget(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i,  j,  0,1, 0)
          call xclget(uline(j,2),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i-1,j,  0,1, 0)
          call xclget(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  0,1, 0)
            if     (ldebug_latbdt .and. mnproc.eq.1) then
            j=jlport(l)
            write(lp,'(a,i4,1p2e12.5)') 'e port, uline:',j,uline(j,1:2)
            write(lp,'(a,i4,1p1e12.5)') 'e port, pline:',j,pline(j)
            write(lp,'(a,i4,1p1e12.5)') 'e port, plnst:',j,plnste(j)
            write(lp,'(a,i4,1p1e12.5)') 'e port, ulnst:',j,ulnste(j)
            endif
          do j= jfport(l),jlport(l)
            crs=                   ulnste(j) -speede(j)*plnste(j)
            fin=1.5*uline(j,1)-0.5*uline(j,2)+speede(j)*pline( j)
            pline(j)  =0.5*(fin-crs)*rspede(j)
            uline(j,1)=    (fin+crs)-uline(j,1)
          enddo
            if     (ldebug_latbdt .and. mnproc.eq.1) then
            j=jlport(l)
            write(lp,'(a,i4,1p2e12.5)') 'e port,   crs:',j,crs,fin
            write(lp,'(a,i4,1p1e12.5)') 'e port, pbavg:',j,pline(j)
            write(lp,'(a,i4,1p1e12.5)') 'e port, ubavg:',j,uline(j,1)
            write(lp,'(a,i4,1p1e12.5)') 'e port, vbavg:',j,vlnste(j)
            write(lp,*)
            call flush(lp)
            endif
          j = jfport(l)
          call xclput(pline(j),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  0,1)
          call xclput(uline(j,1),lnport(l),
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,  0,1)  ! normal
          call xclput(vlnste(j+1),lnport(l)-1,
     &                vbavg(1-nbdy,1-nbdy,n), i,  j+1,0,1)  ! tangential
c
        elseif (kdport(l).eq.1) then
c
c         northern port
c
          j = jfport(l)-1
          i = ifport(l)
          call xclget(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j,  1,0, 0)
          call xclget(uline(i,2),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j-1,1,0, 0)
          call xclget(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  1,0, 0)
            if     (ldebug_latbdt .and. mnproc.eq.1) then
            i=ilport(l)
            write(lp,'(a,i4,1p2e12.5)') 'n port, uline:',i,uline(i,1:2)
            write(lp,'(a,i4,1p1e12.5)') 'n port, pline:',i,pline(i)
            write(lp,'(a,i4,1p1e12.5)') 'n port, plnst:',i,plnstn(i)
            write(lp,'(a,i4,1p1e12.5)') 'n port, vlnst:',i,vlnstn(i)
            endif
          do i= ifport(l),ilport(l)
            crs=                   vlnstn(i) -speedn(i)*plnstn(i)
            fin=1.5*uline(i,1)-0.5*uline(i,2)+speedn(i)*pline( i)
            pline(i)  =0.5*(fin-crs)*rspedn(i)
            uline(i,1)=    (fin+crs)-uline(i,1)
          enddo
            if     (ldebug_latbdt .and. mnproc.eq.1) then
            i=ilport(l)
            write(lp,'(a,i4,1p2e12.5)') 'n port,   crs:',i,crs,fin
            write(lp,'(a,i4,1p1e12.5)') 'n port, pbavg:',i,pline(i)
            write(lp,'(a,i4,1p1e12.5)') 'n port, vbavg:',i,uline(i,1)
            write(lp,'(a,i4,1p1e12.5)') 'n port, ubavg:',i,ulnstn(i)
            write(lp,*)
            call flush(lp)
            endif
          i = ifport(l)
          call xclput(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  1,0)
          call xclput(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j+1,1,0)  ! normal
          call xclput(ulnstn(i+1),lnport(l)-1,
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,  1,0)  ! tangential
c
        elseif (kdport(l).eq.2) then
c
c         southern port
c
          j = jfport(l)
          i = ifport(l)
          call xclget(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j+1,1,0, 0)
          call xclget(uline(i,2),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j+2,1,0, 0)
          call xclget(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  1,0, 0)
          do i= ifport(l),ilport(l)
            crs=                   vlnsts(i) +speeds(i)*plnsts(i)
            fin=1.5*uline(i,1)-0.5*uline(i,2)-speeds(i)*pline( i)
            pline(i)  =0.5*(crs-fin)*rspeds(i)
            uline(i,1)=    (crs+fin)-uline(i,1)
          enddo
          i = ifport(l)
          call xclput(pline(i),  lnport(l),
     &                pbavg(1-nbdy,1-nbdy,n), i,  j,  1,0)
          call xclput(uline(i,1),lnport(l),
     &                vbavg(1-nbdy,1-nbdy,n), i,  j,  1,0)  ! normal
          call xclput(ulnsts(i+1),lnport(l)-1,
     &                ubavg(1-nbdy,1-nbdy,n), i+1,j,  1,0)  ! tangential
c
        endif
c
      enddo  !l=1,nports
c
      return
      end
