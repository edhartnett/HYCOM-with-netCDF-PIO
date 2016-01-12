      subroutine zebra(array,idim,ni,nj)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer idim,ni,nj
      real    array(idim,nj)
c
c --- find nice contour interval resulting in 7 to 10 contour lines and
c --- draw contours on line printer through the following set of grid points:
c
c         array( 1,nj) . . . . . . . . .  array(ni,nj)
c              .                               .
c              .                               .          (plot will appear
c              .                               .          on paper as shown,
c              .                               .          j up, i across)
c              .                               .
c         array( 1, 1) . . . . . . . . .  array(ni, 1)
c
c --- ni may be smaller than idim.
c --- thus, plotting of partial arrays is possible.
c
c --- uses micom version of routine via a work array.
c
      integer i,j
      real, allocatable :: aflip(:,:)
c
      if     (mnproc.eq.1) then  ! 1st processor only
        allocate( aflip(nj,ni) )
        do i= 1,ni
          do j= 1,nj
            aflip(nj+1-j,i) = array(i,j)
          enddo
        enddo
        call zebram(aflip,nj,nj,ni)
        deallocate( aflip )
      endif
      return
      end
c
c
      subroutine zebram(array,idim,ii,jj)
c
c     micom version of zebra.
c
c --- find nice contour interval resulting in 7 to 10 contour lines and
c --- draw contours on line printer through the following set of grid points:
c
c         array( 1, 1) . . . . . . . . .  array( 1,jj)
c              .                               .
c              .                               .          (plot will appear
c              .                               .          on paper as shown,
c              .                               .          i down, j across)
c              .                               .
c         array(ii,jj) . . . . . . . . .  array(ii,jj)
c
c --- ii  may be smaller than  idim, the first (row) dimension of 'array'
c --- in the calling program. thus, plotting of partial arrays is possible.
c
      implicit none
c
      integer idim,ii,jj
      real    array(idim,*)
c
      integer        lp
      common/linepr/ lp
      save  /linepr/
c
      integer i,j
      real    amn,amx,contur,q,ratio
c
      real       sqrt2
      parameter (sqrt2=1.414)
c
      amx=-1.e25
      amn= 1.e25
      do 1 i=1,ii
      do 1 j=1,jj
      amx=max(amx,array(i,j))
 1    amn=min(amn,array(i,j))
c
      if (amx.gt.amn) go to 2
      write (lp,100) array(1,1)
 100  format (//' field to be contoured is constant ...',1pe15.5/)
      return
c
 2    contur=(amx-amn)/6.
      q=10.**int(log10(contur))
      if (contur.lt.1.) q=q/10.
      ratio=contur/q
      if (ratio.gt.sqrt2*5.)  contur=q*10.
      if (ratio.le.sqrt2*5.)  contur=q*5.
      if (ratio.le.sqrt2*2.)  contur=q*2.
      if (ratio.le.sqrt2)     contur=q
      write (lp,101) contur,amn,amx
 101  format (' contour interval in plot below is',1pe9.1,
     .        6x,'min/max =',2e11.3/)
      call digplt(array,idim,ii,jj,contur)
c
      return
      end
c
c
      subroutine digplt(array,idim,ii,jj,dec)
c
c --- simulate a contour line plot on the printer
c
      implicit none
c
      integer idim,ii,jj
      real    array(idim,*),dec
c
      integer        lp
      common/linepr/ lp
      save  /linepr/
c
      character*1 digit(130)
      integer     i,ia,j,ja,k,n
      real        dx,dxdy,dy,value,x,xinc,y,yinc
c
c     nchar = number of character increments in 'j' direction
c     ratio = character width / line spacing
c
      integer    nchar
      real       ratio
      parameter (nchar=120, ratio=0.58)
ccc   parameter (nchar= 74, ratio=0.58)
c
      character*1 dig(20)
      data dig/'0',' ','1',' ','2',' ','3',' ','4',' ',
     .         '5',' ','6',' ','7',' ','8',' ','9',' '/
c
      xinc=float(jj-1)/(float(nchar)*ratio)
      yinc=float(jj-1)/ float(nchar)
      k=float(nchar)*ratio*float(ii-1)/float(jj-1)+1.00001
      do 1 i=1,k
      x=1.+float(i-1)*xinc
      ia=min(ii-1,int(x))
      dx=x-float(ia)
      do 2 j=1,nchar+1
      y=1.+float(j-1)*yinc
      ja=min(jj-1,int(y))
      dy=y-float(ja)
      dxdy=dx*dy
      value=array(ia,ja)*(1.-dx-dy+dxdy)
     .     +array(ia+1,ja)*(dx-dxdy)
     .     +array(ia,ja+1)*(dy-dxdy)
     .     +array(ia+1,ja+1)*dxdy
      n=mod(mod(int(2.*value/dec+sign(.5,value)),20)+20,20)+1
      digit(j)=dig(n)
 2    continue
      write (lp,100) 'i',' ',(digit(j),j=1,nchar+1),' ','i'
 1    continue
 100  format(1x,130a1)
      return
      end
