! contains subroutine mesh()
!----------------------------------------------------------------------!
MODULE setupmod

  implicit none

CONTAINS

  SUBROUTINE  mesh(xmesh,ymesh,myID,ncx_spc)
  ! =================================================================================
  ! This subroutine uses the step variables dx and dy to create an X and Y mesh of the 
  ! spatial grid centers. 
  ! dx,dy,xdomain,ncx,ncy are imported from domainparmod.f90 these are the universal
  ! spatial parameters. The parallelization scheme divides the system with vertical 
  ! lines. That is, the ydomain remains intact but the xdomain is split into nworkers 
  ! sections. The *_spc postfix on the x variables below represent "specific."
  ! ================================================================================
  use domainparmod

  implicit none

  integer,intent(IN)::myID,ncx_spc

  double precision:: yci,xci,xmin2_spc,xmax2_spc
  integer::iter,jter,ndc
  double precision,dimension(ncx_spc)::xcv
  double precision,dimension(ncy)::ycv
  double precision,dimension(ncy,ncx_spc),intent(OUT)::xmesh,ymesh


  !xmin-3.d0*dx/2 is the minimum x of node 0
  xci=(xmin-3.d0*dx/2)+(ncx_spc-4)*myID*dx

                                                                                              
  ! cell center positions for initial x and y
  !xci = xmin_spc-3.d0*dx/2
  yci = ymin-3.d0*dy/2
  !call sleep(1)
  !write(*,*) 'YCI',yci,myID
  ! create a vector with the position of each cell center
  !call sleep(10)

  do iter = 1,ncx_spc
	xcv(iter) = xci + dx*(iter-1)   
  enddo

  do iter = 1,ncy
	ycv(iter) = yci + dy*(iter-1)   
  enddo




  ! create a mesh with cell center positions
  do iter = 1,ncx_spc
	ymesh(:,iter) = ycv
  enddo
  do iter = 1,ncy
   	xmesh(iter,:) = xcv
  enddo
 

  RETURN
  END SUBROUTINE mesh

! ===<>=== ===<>=== ===<>=== ===<>=== ===<>=== ===<>=== ===<>=== ===<>===

SUBROUTINE windprofile(windtype,ymesh,ncx_spc, wind)

  use domainparmod

  implicit none
  integer,intent(IN)::windtype,ncx_spc
  double precision,dimension(ncy,ncx_spc),intent(IN)::ymesh
  double precision,dimension(ncy,ncx_spc),intent(OUT)::wind
  double precision,dimension(ncy,ncx_spc):: onesm
  onesm = ymesh/ymesh ! create a ones matrix the same size as ymesh
  if(windtype.eq.1) then ! jet stream
    	wind = 40.d0*exp(-(ymesh-100000*onesm)**2/(2*10000**2));
  elseif(windtype.eq.2) then
    	wind = 40.d0*tanh(5/1e5*ymesh-5*onesm);
  elseif(windtype.eq.3) then
    	wind = 20.d0-20.d0*tanh(5/1e5*ymesh-5*onesm);
  elseif(windtype.eq.4) then
  	wind = 40.d0*onesm
  else
	wind = 0.d0*onesm
  endif

  RETURN
END SUBROUTINE windprofile

! ===^<>^=== ===^<>^=== ===^<>^=== ===^<>^=== ===^<>^===

SUBROUTINE pressureprofile(pressuretype,xmesh,ymesh,ncx_spc, pressure)

  use domainparmod

  implicit none
  integer,intent(IN)::pressuretype,ncx_spc
  double precision,dimension(ncy,ncx_spc),intent(IN)::xmesh,ymesh
  double precision,dimension(ncy,ncx_spc),intent(OUT)::pressure
  double precision,dimension(ncy,ncx_spc):: onesm
  double precision,parameter:: amp = 10000,sigX=0,sigY=ydomain/10
  onesm = xmesh/xmesh ! create a ones matrix the same size as ymesh
    	
  if(pressuretype.eq.1) then! Gaussian Pressure Pulse
  	pressure = amp*exp(-((xmesh-xmax*onesm/2)/(1.41*sigX))**2) &
		*exp(-((ymesh-ymax*onesm/2)/(1.41*sigY))**2)	
  elseif(pressuretype.eq.2) then
	pressure=(amp/2)*tanh(xmesh)
  else
	pressure=0.0d0*onesm
  endif

  RETURN
END SUBROUTINE pressureprofile

END MODULE setupmod
