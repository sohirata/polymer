module mesh
  implicit none
  integer knr,nr
  parameter (knr=200)
  real(kind=8) ar(knr),rmin,rmax,dr,rhomin
end module mesh

module fact
  implicit none
  real(kind=8) fac(0:80),sgn(-80:80),rttwo,pi
end module fact

program runtrn
  use mesh
  use fact

  implicit none
  integer i,li,mi,lcmax,dim2
  real(kind=8) f(knr),rr(3),aa,ya(knr)
  complex(kind=8), dimension(:,:), allocatable :: tra
  open (6,file='trnout') 

!    initialize some things

  fac(0)=1.0
  sgn(0)=1.0
  do i=1,80
     fac(i)=i*fac(i-1)
     sgn(i)=-sgn(i-1)
     sgn(-i)=sgn(i)
  enddo
  rttwo=sqrt(2.0d0)
  pi=4.0d0*atan(1.0d0)

!   set up the radial mesh, logarithmic coordinates

  nr=100
  rmin=0.001
  rhomin=log(rmin)
  rmax=20.0
  dr=(log(rmax)-log(rmin))/(nr-1)
  do i=1,nr
     ar(i)=rmin*exp((i-1)*dr)
  enddo

!   set up a test case

  li=0
  mi=0
  lcmax=10
  dim2=(lcmax+1)**2-1
  allocate (tra(nr,0:dim2))
  aa=1.0
  do i=1,nr
     f(i)=ar(i)**li*exp(-aa*ar(i))
  enddo
  rr(1)=0.0
  rr(2)=0.0
  rr(3)=1.0

!    run the test case

  call ntrans(f,li,mi,rr,tra,dim2,lcmax)

!    write the real part of the l=0 term

  do i=1,nr
     ya(i)=tra(i,0)
  enddo
  do i=1,nr
    write (6,600) ar(i),f(i),ya(i)
  enddo
600 format(3f10.5)
end program runtrn

subroutine trn(f,rc,li,tx,clz,clzp,nt,lcmax,dim1)

  use mesh
  use fact

!   this subroutine computes the functions given in Eq. (10), apart
!   from the 3j factors

  implicit none

  integer dim1,li,lcmax,ll1,lc2min,lc2max,ic,i,lc1,lc2,ll,ngl,&
       j,k,lbdmax,lbd,ll2,lbdmin,nt,clz(dim1),clzp(dim1)
  real(kind=8) tx(nr,dim1),f(nr),rc,ya(knr),xx,rmn,rmx,rhx,dlt,dy,yy,zz,&
       ss,yb(50),qq1,qq2,qq3,aa
  real(kind=8) plx(100,0:30),tb(knr,0:30),ta(knr,0:10),xgl(50),&
       wgl(50)

  ngl=20
  call gauleg(xgl,wgl,ngl)
  lbdmax=lcmax+li
  if (rc.lt.0.00001d0) then
     nt=1
     clz(1)=li
     clzp(1)=0
     do i=1,nr
        tx(i,1)=f(i)*(2*li+1)
     enddo
     return
  endif

!   get the l=0 factor in the function

  do i=1,nr
     ya(i)=f(i)/ar(i)**li
  enddo

!   compute the functions defined in Eq. (12)

  do i=1,nr

!   table ta contains powers ot ar(i), see Eq. (13)
!   table tb contains function F_lambda, see Eq. (12)

     xx=ar(i)/rc
     ta(i,0)=rc**li
     if (li.gt.0) then
        do ll=1,li
           ta(i,ll)=xx*ta(i,ll-1)
        enddo
     endif
     rmn=min(rc,ar(i))
     rmx=max(rc,ar(i))
     do j=1,ngl   
        xx=rmx+xgl(j)*rmn

!   Here is where the interpolation happens

        rhx=log(xx)
        k=(rhx-rhomin)/dr+1
        k=min(k,nr-3)
        k=max(k,3)
        dlt=rhx-log(ar(k))
        dy=dlt/dr

        yy=(-dy*(dy**2-1.0d0)*(dy-2.0d0)*(dy-3.0d0)*ya(k-2)&
             +5.0d0*dy*(dy-1.0d0)*(dy**2-4.0d0)*(dy-3.0d0)*ya(k-1)&
             -10.0d0*(dy**2-1.0d0)*(dy**2-4.0d0)*(dy-3.0d0)*ya(k)&
             +10.0d0*dy*(dy+1.0d0)*(dy**2-4.0d0)*(dy-3.0d0)*ya(k+1)&
             -5.0d0*dy*(dy**2-1.0d0)*(dy+2.0d0)*(dy-3.0d0)*ya(k+2)&
             +dy*(dy**2-1.0d0)*(dy**2-4.0d0)*ya(k+3))/120.0d0
        yb(j)=xx*yy
     enddo

!    This is where the integration in Eq. (12) happens

     do j=1,ngl
        zz=-xgl(j)+0.5d0*rmn*(1.0d0-xgl(j)**2)/rmx
        plx(j,0)=1.0d0
        if (lbdmax.gt.0) plx(j,1)=zz
        if (lbdmax.gt.1) then
           do lbd=1,lbdmax-1
              plx(j,lbd+1)=((2*lbd+1)*zz*plx(j,lbd)-lbd*plx(j,lbd-1))/(lbd+1)
           enddo
        endif
     enddo
     do lbd=0,lbdmax
        ss=0.0d0
        do j=1,ngl
           ss=ss+wgl(j)*plx(j,lbd)*yb(j)
        enddo
        tb(i,lbd)=0.5d0*(2*lbd+1)*ss/rmx
     enddo
  enddo
  do ic=1,dim1
     do i=1,nr
        tx(i,ic)=0.0d0
     enddo
  enddo
  ic=0
  
!   This carries out the summation implied in Eq. (16)
  do lc1=0,lcmax
     lc2min=abs(lc1-li)
     lc2max=lc1+li
     do lc2=lc2min,lc2max,2
        ic=ic+1
        do ll1=0,li
           ll2=li-ll1
           lbdmin=max(abs(lc1-ll1),abs(lc2-ll2))
           lbdmax=min(lc1+ll1,lc2+ll2)
           if (lbdmin.le.lbdmax) then
              do lbd=lbdmin,lbdmax,2
                 qq1=(fac((lc1+lbd+ll1)/2)*fac(lc1+lbd-ll1))&
                      /(fac((lc1+lbd-ll1)/2)*fac((ll1+lbd-lc1)/2)&
                      *fac((lc1+ll1-lbd)/2)*fac(lc1+lbd+ll1+1))
                 qq3=(fac((lc1+lc2+li)/2)*fac(lc1+lc2-li))&
                      /(fac((lc1+lc2-li)/2)*fac((lc1+li-lc2)/2)&
                      *fac((lc2+li-lc1)/2)*fac(lc1+lc2+li+1))
                 qq2=(fac((lc2+lbd+ll2)/2)*fac(lc2+lbd-ll2))&
                      /(fac((lc2+lbd-ll2)/2)*fac((lc2+ll2-lbd)/2)&
                      *fac((ll2+lbd-lc2)/2)*fac(lc2+lbd+ll2+1))
                 aa=sgn(ll2)*(2*lc1+1)*(2*lc2+1)*qq1*qq2/qq3
                 if (aa.ne.0.0d0) then
                    clz(ic)=lc1
                    clzp(ic)=lc2
                    do i=1,nr
                       tx(i,ic)=tx(i,ic)+aa*ta(i,ll1)*tb(i,lbd)
                    enddo
                 endif
              enddo
           endif
        enddo
     enddo
  enddo
  nt=ic
  return
end subroutine trn

subroutine ntrans(f,li,mi,rr,tra,idim,lcmax)

!   This subroutine carries out the sum on L' and M' to get the 
!   translation functions in complex array TRA 

  use mesh
  use fact
  implicit none

  integer lcmax,idim,li,mi,lcmaxp,ic,iy,lc,lcp,mc,mcp,&
       i,ix,ll2,nt,ll1,lc1,lc2min,lc2max,lc2,dim1
  real(kind=8) f(*),rr(3),rc,rl,thrj
  complex(kind=8) cc,tra(nr,0:idim)
  real(kind=8), dimension (:,:), allocatable :: tx
  complex(kind=8), dimension (:), allocatable :: cylm
  integer, dimension (:), allocatable :: clz,clzp
  lcmaxp=lcmax+li
  ll1=(lcmax+1)**2
  ll2=(lcmaxp+1)**2-1
  rc=sqrt(rr(1)**2+rr(2)**2+rr(3)**2)
  ic=0
  do lc1=0,lcmax+li
     lc2min=abs(lc1-li)
     lc2max=lc1+li
     do lc2=lc2min,lc2max,2
        ic=ic+1
     enddo
  enddo
  dim1=ic
  allocate (tx(nr,dim1),clz(dim1),clzp(dim1),cylm(0:ll2))
  iy=-1
  do lc=0,lcmax
     do mc=-lc,lc
        iy=iy+1
        do i=1,nr
           tra(i,iy)=0.0d0
        enddo
     enddo
  enddo
  if (rc.lt.1.0d-6) then
     iy=li*(li+1)+mi
     do i=1,nr
        tra(i,iy)=f(i)
     enddo
     return
  endif
  call csphar(rr,cylm,lcmaxp)
  call trn(f,rc,li,tx,clz,clzp,nt,lcmax,dim1)
  do ix=1,nt
     lc=clz(ix)
     lcp=clzp(ix)
     do mc=-lc,lc
        mcp=-mc-mi
        ic=lcp*(lcp+1)+mcp
        cc=sgn(mi)*thrj(lc,lcp,li,0,0,0)*thrj(lc,lcp,mi,mc,mcp,mi)*cylm(ic)
        iy=lc*(lc+1)+mc
        do i=1,nr
           tra(i,iy)=tra(i,iy)+cc*tx(i,ix)
        enddo
     enddo
  enddo
  return
600 format (7f10.5)
end subroutine ntrans

subroutine gauleg(x,w,n)
  implicit none
  integer n
  real(kind=8) eps,pi,z,p1,p2,p3,pp,z1,x(*),w(*)
  integer i,m,j
  eps=3.0d-14
  m=(n+1)/2
  pi=4.0d0*atan(1.0d0)
  do i=1,m
     z=cos(pi*(i-0.25d0)/(n+0.5d0))
1    continue
     p1=1.0d0
     p2=0.0d0
     do j=1,n
        p3=p2
        p2=p1
        p1=((2.0d0*j-1)*z*p2-(j-1)*p3)/j
     enddo
     pp=n*(z*p1-p2)/(z*z-1.0d0)
     z1=z
     z=z1-p1/pp
     if (abs(z-z1).gt.eps) go to 1
     x(i)=-z
     x(n+1-i)=z
     w(i)=2.0d0/((1.0d0-z*z)*pp*pp)
     w(n+1-i)=w(i)
  enddo
  return
end subroutine gauleg

subroutine csphar(r,res,lmax) 
  use fact

  implicit none 
  real(kind=8) r(3),x,y,z,dd,phi,cc,ss,al,aa,bb,zz,cs
  integer lmax,ll,l,m,il1,il2,ind,ll2
  complex(kind=8) res(0:*)
  x=r(1) 
  y=r(2) 
  z=r(3) 
  pi=4.0d0*atan(1.0d0)
  rttwo=sqrt(2.0d0)
  dd=sqrt(x*x+y*y+z*z)
  if (dd.lt.1.0d-10) then
     ll=(lmax+1)**2-1 
     do  l=1,ll 
        res(l)=0.0d0 
     end do
     res(0)=1.0d0 
     return
  endif
  if (abs(x).lt.0.00001) then
     phi=0.5d0*pi
     if (y.lt.0.0d0) phi=-phi 
  else
     phi=atan(y/x) 
     if (x.lt.0.0d0) phi=phi+pi 
  endif
  ss=sqrt(x*x+y*y)/dd 
  cc=z/dd
  res(0)=1.0d0 
  if (lmax.eq.0) return
  do l=1,lmax 
     al=l 
     il2=(l+1)**2-1 
     il1=l**2-1
     res(il2)=-ss*sqrt((al-0.5d0)/al)*res(il1) 
     res(il2-1)=cc*sqrt(2.0d0*al-1.0d0)*res(il1)
  end do
  if (lmax.ge.2) then
     do m=0,lmax-2
        if (m.lt.lmax) then
           do l=m+1,lmax-1
              ind=l*(l+1)+m 
              aa=l**2-m**2
              bb=(l+1)**2-m**2
              zz=(2*l+1)*cc*res(ind)-sqrt(aa)*res(ind-2*l) 
              res(ind+2*(l+1))=zz/sqrt(bb) 
           end do
        endif
     end do
  endif
  do l=0,lmax
     ll2=l*(l+1)
     do m=0,l
        cs=sin(m*phi)
        cc=cos(m*phi)
        res(ll2+m)=cmplx(cc,cs)*res(ll2+m)
        res(ll2-m)=sgn(m)*conjg(res(ll2+m))
     enddo
  enddo
  return 
end subroutine csphar

function thrj(l1,l2,l3,m1,m2,m3) 
  use fact

  implicit none
  real(kind=8) thrj,xx,cc,s,ph,t,dlt
  integer l1,l2,l3,m1,m2,m3,lg,it,lgh,itmin,itmax
  thrj=0.0d0 
  if (m1**2.gt.l1**2) return 
  if (m2**2.gt.l2**2) return 
  if (m3**2.gt.l3**2) return 
  if (m1+m2+m3.ne.0) return 
  if (l3.lt.iabs(l1-l2)) return 
  if (l3.gt.l1+l2) return 
  lg=l1+l2+l3
  dlt=sqrt(fac(lg-2*l1)*fac(lg-2*l2)*fac(lg-2*l3)/fac(lg+1))
  if ((m1.eq.0).and.(m2.eq.0)) then
     if (mod(lg,2).eq.1) return
     lgh=lg/2
     thrj=sgn(lgh)*dlt*fac(lgh)/(fac(lgh-l1)*fac(lgh-l2)*fac(lgh-l3))
     return
  endif
  xx=fac(l3+m3)*fac(l3-m3)/(fac(l1+m1)*fac(l1-m1)&
       *fac(l2+m2)*fac(l2-m2)) 
  cc=dlt*sqrt(xx)
  itmin=max0(0,l1-l2+m3) 
  itmax=min0(l3-l2+l1,l3+m3) 
  s=0.0d0 
  ph=1.0d0 
  it=itmin 
  if (mod(it+l2+m2,2).eq.1) ph=-ph 
12 t=ph*fac(l3+l1-m2-it)*fac(l2+m2+it) 
  t=t/(fac(l3+m3-it)*fac(it+l2-l1-m3)*fac(it)*fac(l3-l2+l1-it)) 
  s=s+t 
  it=it+1 
  ph=-ph 
  if (it.le.itmax) go to 12 
  thrj=cc*s 
  return 
end function thrj

