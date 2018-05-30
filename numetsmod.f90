! contains subroutines xflux() and yflux
!----------------------------------------------------------------------!
MODULE nummetmod

  implicit none

CONTAINS

  SUBROUTINE xflux(Q,ny,nxs, flux)
  use initparamod

  implicit none  
  integer,intent(IN):: nxs,ny
  double precision,dimension(ny,nxs)::P
  double precision,dimension(ny,nxs,4),intent(IN)::Q
  double precision,dimension(ny,nxs,4),intent(OUT)::flux

  P = (gamm-1)*(Q(:,:,4)-0.5*(Q(:,:,2)**2+Q(:,:,3)**2)/Q(:,:,1))
  flux(:,:,1) = Q(:,:,2)
  flux(:,:,2) = (Q(:,:,2)**2)/Q(:,:,1)+P
  flux(:,:,3) = Q(:,:,2)*Q(:,:,3)/Q(:,:,1) 
  flux(:,:,4) = (Q(:,:,4)+P)*(Q(:,:,2)/Q(:,:,1))
  RETURN

  END SUBROUTINE


! ===<>==== ===<>==== ===<>==== ===<>==== ===<>==== ===<>==== ===<>====
  SUBROUTINE  yflux(Q,ny,nxs, flux)
  use initparamod

  implicit none
  integer,intent(IN):: ny,nxs  
  double precision,dimension(ny,nxs)::P
  double precision,dimension(ny,nxs,4),intent(IN)::Q
  double precision,dimension(ny,nxs,4),intent(OUT)::flux


  P = (gamm-1)*(Q(:,:,4)-0.5*(Q(:,:,2)**2+Q(:,:,3)**2)/Q(:,:,1))
  flux(:,:,1) = Q(:,:,3)
  flux(:,:,2) = (Q(:,:,2)*Q(:,:,3))/Q(:,:,1)
  flux(:,:,3) = (Q(:,:,3)**2)/Q(:,:,1) + P 
  flux(:,:,4) = (Q(:,:,4)+P)*(Q(:,:,3)/Q(:,:,1))
  RETURN

  END SUBROUTINE 

! ===<>==== ===<>==== ===<>==== ===<>==== ===<>==== ===<>==== ===<>====
  SUBROUTINE  source(Q,gravitymesh,ny,nxs, S)
  use domainparmod

  implicit none
  integer,intent(IN)::ny,nxs
  double precision,dimension(ny,nxs),intent(IN)::gravitymesh
  double precision,dimension(ny,nxs,4),intent(IN)::Q
  double precision,dimension(ny,nxs,4),intent(OUT)::S

  S(:,:,1) = 0*Q(:,:,1)
  S(:,:,2) = 0*Q(:,:,2)
  S(:,:,3) = -Q(:,:,1)*gravitymesh
  S(:,:,4) = -Q(:,:,3)*gravitymesh

  RETURN

  END SUBROUTINE 

! ===<>==== ===<>==== ===<>==== ===<>==== ===<>==== ===<>==== ===<>====
  SUBROUTINE  sourcemnt(ymesh,xmesh,nxs,t,dt, Q)
  use domainparmod
  use initparamod

  implicit none
  integer,intent(IN)::nxs
  double precision,intent(IN)::t,dt
  double precision,dimension(ncy,nxs)::Ff
  double precision,dimension(ncy,nxs),intent(IN)::ymesh,xmesh
  double precision,dimension(ncy,nxs,4),intent(INOUT)::Q
  double precision,dimension(ncy,nxs):: onesm
  double precision,dimension(ncy,nxs)::ys,xs
  double precision::ts
  double precision,parameter:: hk = 8*atan(1.d0)/xdomain
  double precision,parameter:: mwamp = 0.134 ! wave max amplitude [m/s/s]
  double precision,parameter:: ycent = 10000, xcent = 0 ! 
  double precision,parameter:: tcent = 6.d0*60*60 ! 
  double precision,parameter:: sigy = 3000, sigx = sqrt(0.6*xdomain) ! 
  double precision,parameter:: sigt = 1.5*60*60 ! [s]
  double precision,parameter:: sy2sq=-2*(sigy)**2,sx2sq=-2*(sigy)**2;
  double precision,parameter:: st2sq=-2*(sigt)**2

  onesm = ymesh/ymesh
  ys = ymesh-onesm*ycent
  xs = xmesh-onesm*xcent
  ts = t-tcent
!! Mountain Wave Forcing Function [see Heale,Snively]
  Ff = mwamp*Q(:,:,1)*cos(hk*xs) &
        *exp(ys**2/sy2sq + xs**2/sx2sq + ts**2/st2sq);
  Q(:,:,3) = Q(:,:,3)+Ff*dt

  RETURN

  END SUBROUTINE 

END MODULE nummetmod
