!Uses Laplacian for viscous terms
!Uniform LM in core. P2P-1 scheme
!Separate (block diagonal) preconditioners for pressure and LM terms in the	form GTG and BTB
!Crank-Nicholson
	  implicit double precision (a-h,o-z)
	  parameter (nx=61,ny=61,nz=61,nel=(nx-1)*(ny-1)*(nz-1)/8)
	  parameter (nxp=(nx+1)/2,nyp=(ny+1)/2,nzp=(nz+1)/2)
	  parameter (nvu=nx*ny*nz,ned=(nx-1)/2)
	  parameter (nup=(nx+ny-1)*nz)
	  parameter (nprt=322,nbp=106,nce=50)
!
	  dimension vgco(nvu,8)
!Velocity
	  dimension mgu(nel,27),inu(0:2,0:2,0:2)
	  dimension igiu(nvu)
	  dimension u(nvu),v(nvu),w(nvu)
	  dimension bu(nvu),bv(nvu),bw(nvu)
	  dimension bun(nvu,2),bvn(nvu,2),bwn(nvu,2)
	  dimension u0(nvu),v0(nvu),w0(nvu)
	  dimension iup(nup,2)
	  dimension amu(27,27),alu(27,27),alut(27,27),alx(27,27),aly(27,27),alz(27,27),dul(nvu)
	  dimension ctux(27,27,27),ctuy(27,27,27),ctuz(27,27,27)
	  dimension fi(0:2),dfi(0:2)
!Pressure
	  dimension p(nel,0:3),rpr(nel,0:3),ppr(nel,0:3),zpr(nel,0:3)
	  dimension pgm(27)
      dimension pu(0:3,27),pv(0:3,27),pw(0:3,27)
!Particles
      dimension rot(nprt,3),up(nprt,3)
      dimension qbp(nbp,3)
	  dimension xp(nprt),yp(nprt),zp(nprt)
      dimension fx(nprt),fy(nprt),fz(nprt)
	  dimension iebp(nprt,nbp),bf(nprt,nbp,27)
	  dimension iece(nprt,nce)
	  dimension ac(nprt,nce,3),ab(nprt,nbp,3)
	  dimension rc(nprt,nce,3),rb(nprt,nbp,3)
	  dimension rcm(nprt,3),rbm(nprt,3),rrb(nprt,3)
	  dimension pc(nprt,nce,3),pb(nprt,nbp,3)
	  dimension zc(nprt,nce,3),zb(nprt,nbp,3)
	  dimension au(27),av(27),aw(27)
	  dimension xt(3)
!Shear lagrange multipliers
      dimension ts(ned**2,3),rs(ned**2,3),zs(ned**2,3),ps(ned**2,3)
	  dimension fi1(0:2),fi2(0:2),fl(9),fu1(9),fu2(9)
	  character*3 sim(300)
          integer ios
!
	  character*18 filename
!
!
	  pi=-1.d0
	  pi=dacos(pi)
!
if (nx==61) xpl=15.d0
if (nx==41) xpl=10.d0
	  ypl=xpl
	  zpl=xpl
!
	  radius=1.d0
!
      dt=.001d0
!	  tmax=60.d0
!	  Re=1.d0
      nss=1
      dtf=dt/dfloat(nss)
CALL GETARG(1, filename) 
!filename='prt1.txt'
nol=0
open(3,type='old',file=filename)
do j=1,100000
read(3,*,iostat=ios)
if (ios/=0) exit
nol=nol+1
enddo
close(3)
!print*, nol
icoun=nol/(nprt+1)
!
!Initial particle positions
!
!open(3,type='old',file='315.txt')  !@@@@@@@
!read(3,*) xp,yp,zp
! close(unit=3)
! open(7,type='old',file='315Re1uvw.txt')  !@@@@@@@
!      read(7,*) tinitial
!          read(7,*) xp,yp,zp
!          read(7,*) up
!          read(7,*) u,v,w
!          read(7,*) p
!          read(7,*) ab,ac,ts

!open(12,type='new',file='pit'//filename)
!do iiii=1,icoun
!open(4,type='old',file=filename)
!  read(4,*) t
!         read(4,*) xp,yp,zp
 

!Box
      bs=2.d0*radius/dsqrt(3.d0)
	  nby=dint(ypl/bs)
	  nbz=nby	 !??????
	  nbx=nby  !??????
      dxg=ypl/dfloat(nby)
!
	  dx=xpl/dfloat(nx-1)*2.d0
	  dy=ypl/dfloat(ny-1)*2.d0
	  dz=zpl/dfloat(nz-1)*2.d0
!
do iiii=1,icoun
open(4,type='old',file=filename)
  read(4,*) t
         read(4,*) xp,yp,zp

      xus=modulo(t*xpl,xpl)
	  iedu1=dint(xus/dx)
	  xso=xus/dx-dfloat(iedu1)
	  !fi2(0)=xso-3.d0*xso**2/2.d0+2.d0*xso**3/3.d0
	  !fi2(1)=2.d0*xso-4.d0*xso**3/3.d0
	  !fi2(2)=-xso**2/2.d0+2.d0*xso**3/3.d0
     ! fi1(0)=fi(0)-fi2(0)
      !fi1(1)=fi(1)-fi2(1)
      !fi1(2)=fi(2)-fi2(2)
!
!	  indx=0
!	  do ia=0,2
!	     do ja=0,2
!		    indx=indx+1
!			fl(indx)=fi(ia)*fi(ja) !*dx**2
!			fu1(indx)=fi1(ia)*fi(ja) !*dx**2
!			fu2(indx)=fi2(ia)*fi(ja) !*dx**2
!		 enddo
!	  enddo
!
!Move particles
!
      do k=1,1
	    call forces(nprt,xp,yp,zp,fx,fy,fz,xpl,ypl,zpl,nbx,nby,nbz,radius,dxg,xus,dstmn)
!
!	    up(1:nprt,1)=up(1:nprt,1)+fx*dtf
!	    up(1:nprt,2)=up(1:nprt,2)+fy*dtf
!	    up(1:nprt,3)=up(1:nprt,3)+fz*dtf
!	    yp=modulo(yp+up(1:nprt,2)*dtf,ypl)
!		do ip=1,nprt
!		   zt=zp(ip)+up(ip,3)*dtf
!	       zp(ip)=modulo(zp(ip)+up(ip,3)*dtf,zpl)
!		   disp=(zp(ip)-zt)/zpl*xus
!	       xp(ip)=modulo(xp(ip)+up(ip,1)*dtf+disp,xpl)
!		enddo
      enddo
!
!      ikpp=ikpp+1
 !!     if(ikpp.ge.100) then
!	     ikpp=0
 !        icount=icount+1
!	     filename='sre1'//sim(icount)//'.txt'
		 !write(*,*) t
         !write(*,*) xp,yp,zp
!	  endif
write(*,*) dstmn
!write(12,*) dstmn
enddo
999 close(unit=3)
close(12)
close(4)
end program
!------------------------------------------------------------------------------------------------------
      subroutine forces(nprt,xp,yp,zp,fx,fy,fz,xpl,ypl,zpl,nbx,nby,nbz,a,dx,xu,distmin)
	  implicit double precision (a-h,o-z)
      dimension xp(nprt),yp(nprt),zp(nprt)
	  dimension fx(nprt),fy(nprt),fz(nprt),distm(nprt,nprt)
	  dimension iiap(nprt,3),iap(0:nbx-1,0:nby-1,0:nbz-1)
          fx=0.d0
	  fy=0.d0
	  fz=0.d0
!
	  r=.1d0*a
	  ep=.1d0
! 	  
	  iap=0
	  iapp=0
	  sepmin=1.d0
distm=25.d0
!icountt=0
!Find particle position on grid
	  do ip=1,nprt
	     i=dint(xp(ip)/dx)
		 j=dint(yp(ip)/dx)
		 k=dint(zp(ip)/dx)
		 iiap(ip,1)=i
		 iiap(ip,2)=j
		 iiap(ip,3)=k
		 iap(i,j,k)=ip
	  enddo
	  index=0
	  do i=0,nbx-1
	     do j=0,nby-1
		    do k=0,nbz-1
			   if(iap(i,j,k).gt.0) index=index+1
			enddo
		 enddo
	  enddo
!Find nearest neighbors
      do ip=1,nprt-1
	     ix=iiap(ip,1)
		 iy=iiap(ip,2)
		 iz=iiap(ip,3)
		 do i1=ix-2,ix+2
		    i=modulo(i1,nbx)
		    lx=(i1-i)/nbx
			do j1=iy-2,iy+2
		       j=modulo(j1,nby)
		       ly=(j1-j)/nby
			   do k1=iz-2,iz+2
		          k=modulo(k1,nbz)
		          lz=(k1-k)/nbz
			      m=iap(i,j,k)
                                if (m.gt.ip) then
			         vnx=xp(m)+dfloat(lx)*xpl+dfloat(lz)*xu-xp(ip)
			         vny=yp(m)+dfloat(ly)*ypl-yp(ip)
			         vnz=zp(m)+dfloat(lz)*zpl-zp(ip)
!
		             dist=dsqrt(vnx**2+vny**2+vnz**2)

!icountt=icountt+1
!if (dist.lt.distminm(m,ip)) then
distm(m,ip)=dist
!endif
		         sep=dist-2.d0*a
					 if(sep.lt.sepmin) sepmin=sep
!
			         if(sep.lt.r) then
						ft=(1.d0-sep/r)**2/ep
			            fxa=ft*vnx/dist
			            fya=ft*vny/dist
			            fza=ft*vnz/dist
!
			            fx(ip)=fx(ip)-fxa
			            fy(ip)=fy(ip)-fya
			            fz(ip)=fz(ip)-fza
			            fx(m)=fx(m)+fxa
			            fy(m)=fy(m)+fya
			            fz(m)=fz(m)+fza
				     endif
			      endif
		       enddo
		    enddo
         enddo
	  enddo
!do i=1,nprt
distmin= minval(distm(:,:))
!write(*,*) minval(distm(:,:)),minloc(distm(:,:))
!write(*,*) minval(distm(i,:)),i
!endif
!if (minval(distm(:,i)).ne.distmin) then
!write(*,*) minval(distm(:,i)),i
!endif

!enddo
999  return
      end
!---------------------------------------------------------------------------------------------------------
