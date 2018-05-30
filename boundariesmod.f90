! contains subroutines bcb(), bct(), and bcs()
! bcs() uses mpi
!----------------------------------------------------------------------!
MODULE boundariesmod

  implicit none

CONTAINS

  SUBROUTINE  bcb(rho,rhoU,rhoV,energy,P0,rho0,wind,ncxs) ! boundary condition closed bottom
  ! use initparamod
  use domainparmod
  use initparamod

  implicit none
  integer,intent(IN)::ncxs
  double precision,dimension(ncy,ncxs),intent(IN):: P0,rho0,wind
  double precision,dimension(ncy,ncxs),intent(INOUT):: rho,rhoU,rhoV,energy
  double precision,dimension(2,ncxs)::tmp1,rhot,wndt,prst

  tmp1(1,:) = rho(3,:); tmp1(2,:) = rho(3,:) ! equivalent to [rho(3,:); rho(3,:)]-------- tmp1
  rhot(1,:) = rho0(3,:); rhot(2,:) = rho0(3,:) ! equivalent to [rho0(3,:); rho0(3,:)]---- rhot
  rho(1:2,:) = rho0(1:2,:)+(tmp1-rhot)&
	*sqrt(rho0(1:2,:)/rhot)

  wndt(1,:) = wind(3,:); wndt(2,:) = wind(3,:) ! :*( equivalent to [wind(3,:); wind(3,:)]---- wndt 
  tmp1(1,:) = rhoU(3,:); tmp1(2,:) = rhoU(3,:) ! equivalent to [rhoU(3,:); rhoU(3,:)]---- tmp1  
  rhoU(1:2,:) = rho0(1:2,:)*wind(1:2,:)&
	+(tmp1-rhot*wndt)&
	*sqrt(rho0(1:2,:)/rhot)

  tmp1(1,:) = rhoV(3,:); tmp1(2,:) = rhoV(3,:) ! equivalent to [rhoV(3,:); rhoV(3,:)]---- tmp1  
  rhoV(1:2,:) = -tmp1*sqrt(rho0(1:2,:)/rhot)
  
  tmp1(1,:) = energy(3,:); tmp1(2,:) = energy(3,:) ! eqvlnt [energy(3,:); energy(3,:)]--- tmp1  
  prst(1,:) = P0(3,:); prst(2,:) = P0(3,:)         ! equivalent to [P0(3,:); P0(3,:)]---- prst 
  energy(1:2,:)=P0(1:2,:)/(gamm-1)& 
	+(tmp1 - prst/(gamm-1) - 0.5*rhot*wndt**2) &
	*sqrt(rho0(1:2,:)/rhot)+(0.5)*rho0(1:2,:)*wind(1:2,:)**2

  return
  end subroutine bcb 

  SUBROUTINE  bct(rho,rhoU,rhoV,energy,P0,rho0,wind,ncxs) ! boundary condition open top
  !use initparamod
  use domainparmod
  use initparamod

  implicit none
  integer,intent(IN)::ncxs
  double precision,dimension(ncy,ncxs),intent(IN):: P0,rho0,wind
  double precision,dimension(ncy,ncxs),intent(INOUT)::rho,rhoU,rhoV,energy
  double precision,dimension(2,ncxs)::tmp1,rhot,wndt,prst ! arrays used to simplify matrix operations

  tmp1(1,:) = rho(ncy-2,:); tmp1(2,:) = rho(ncy-2,:) ! equivalent to [rho(end-2,:); rho(end-2,:)]--- tmp1
  rhot(1,:) = rho0(ncy-2,:); rhot(2,:) = rho0(ncy-2,:) ! eqvalnt to [rho0(end-2,:); rho0(end-2,:)]-- rhot
  rho(ncy-1:ncy,:) = rho0(ncy-1:ncy,:)+(tmp1-rhot)&
	*sqrt(rho0(ncy-1:ncy,:)/rhot) 
  
  tmp1(1,:) = rhoU(ncy-2,:); tmp1(2,:) = rhoU(ncy-2,:)! equivlnt to [rhoU(end-2,:); rhoU(end-2,:)]-- tmp1
  rhoU(ncy-1:ncy,:) = rho0(ncy-1:ncy,:)*wind(ncy-1:ncy,:)&
	+(tmp1-rhot*wndt)*sqrt(rho0(ncy-1:ncy,:)/rhot)

  tmp1(1,:) = rhoV(ncy-2,:); tmp1(2,:) = rhoV(ncy-2,:)! equivlnt to [rhoV(end-2,:); rhoV(end-2,:)]-- tmp1
  rhoV(ncy-1:ncy,:) = tmp1*sqrt(rho0(ncy-1:ncy,:)/rhot)
  
  tmp1(1,:) = energy(ncy-2,:); tmp1(2,:) = energy(ncy-2,:)! eqvl [energy(end-2,:); energy(end-2,:)]- tmp1
  prst(1,:) = P0(ncy-2,:); prst(2,:) = P0(ncy-2,:)     ! equivalent to [P0(end-2,:); P0(end-2,:)]--- prst
  wndt(1,:) = wind(ncy-2,:); wndt(2,:) = wind(ncy-2,:)! equivlnt to [wind(end-2,:); wind(end-2,:)]-- wndt
  energy(ncy-1:ncy,:) = P0(ncy-1:ncy,:)/(gamm-1) &
	+(tmp1-prst/(gamm-1)-(0.5)*rhot*wndt**2) * &
	sqrt(rho0(ncy-1:ncy,:)/rhot) * &
	(0.5)*rho0(ncy-1:ncy,:)*wind(ncy-1:ncy,:)**2 

  return
  end subroutine bct 

! ===<>===<>=== ===<>===<>=== ===<>===<>=== ===<>===<>=== ===<>===<>=== ===<>===<>=== ===<>===<>=== ===<>===<>===

  SUBROUTINE  bcs(rho,rhoU,rhoV,energy,ncxs,nprocs,myID) ! boundary condition periodic sides

  use domainparmod
  use mpi

  implicit none
  integer,intent(IN)::ncxs,nprocs,myID
  integer:: myIDt ! temporary id used for logic operations
  double precision,dimension(ncy,ncxs),intent(INOUT)::rho,rhoU,rhoV,energy


  integer:: mpi_doub_array,msgtag,ierr
  integer:: rnode,lnode,iter,jter,kter,lter
  !integer::status(MPI_STATUS_SIZE)
  integer,dimension(4):: indices

  if(myID.eq.0) then !======================
	lnode = nprocs-1; rnode = 1
  elseif(myID.eq.(nprocs-1)) then
	lnode = myID-1; rnode = 0
  else
	lnode = myID-1; rnode = myID+1
  endif !=============================

  
  ! Documentation for the messy DO,IF,DO loops may be found at the bottom of the page
  myIDt = myID
  DO kter = 1,2 ! <><><> <><><> <><><> <><><> <><><> <><><> <><><> <><><> <><><> <><><> <><><> <><><>
  if(mod(myIDt,2).eq.0) then
  	indices(1:4) = (/ncxs-2,ncxs-3,4,3/)
	do iter = 1,2
	   jter = indices(iter)
	   call MPI_SEND(rho(:,jter),ncy,MPI_DOUBLE_PRECISION,rnode,1,MPI_COMM_WORLD,ierr)
	   call MPI_SEND(rhoU(:,jter),ncy,MPI_DOUBLE_PRECISION,rnode,2,MPI_COMM_WORLD,ierr)
	   call MPI_SEND(rhoV(:,jter),ncy,MPI_DOUBLE_PRECISION,rnode,3,MPI_COMM_WORLD,ierr)
	   call MPI_SEND(energy(:,jter),ncy,MPI_DOUBLE_PRECISION,rnode,4,MPI_COMM_WORLD,ierr)
	enddo ! iter loop
  	do iter = 3,4
	   jter = indices(iter)
	   call MPI_SEND(rho(:,jter),ncy,MPI_DOUBLE_PRECISION,lnode,5,MPI_COMM_WORLD,ierr)
	   call MPI_SEND(rhoU(:,jter),ncy,MPI_DOUBLE_PRECISION,lnode,6,MPI_COMM_WORLD,ierr)
	   call MPI_SEND(rhoV(:,jter),ncy,MPI_DOUBLE_PRECISION,lnode,7,MPI_COMM_WORLD,ierr)
	   call MPI_SEND(energy(:,jter),ncy,MPI_DOUBLE_PRECISION,lnode,8,MPI_COMM_WORLD,ierr)
	enddo
  else! seelseelseelseelseelseelseelseelseelseelseelseelseelseelseelseelseelseelseelseelseelseelseelse

  	indices(1:4) = (/2,1,ncxs,ncxs-1/)	
	do iter = 1,2
	   jter = indices(iter)
	   call MPI_RECV(rho(:,jter),ncy,MPI_DOUBLE_PRECISION,lnode,1,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
	   call MPI_RECV(rhoU(:,jter),ncy,MPI_DOUBLE_PRECISION,lnode,2,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
	   call MPI_RECV(rhoV(:,jter),ncy,MPI_DOUBLE_PRECISION,lnode,3,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
	   call MPI_RECV(energy(:,jter),ncy,MPI_DOUBLE_PRECISION,lnode,4,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
	enddo

	do iter = 3,4
	   jter = indices(iter)
	   call MPI_RECV(rho(:,jter),ncy,MPI_DOUBLE_PRECISION,rnode,5,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
	   call MPI_RECV(rhoU(:,jter),ncy,MPI_DOUBLE_PRECISION,rnode,6,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
	   call MPI_RECV(rhoV(:,jter),ncy,MPI_DOUBLE_PRECISION,rnode,7,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
	   call MPI_RECV(energy(:,jter),ncy,MPI_DOUBLE_PRECISION,rnode,8,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
	enddo

  endif
  myIDt = myIDt+1 ! bunch of bullshit if you ask me (or see below)
  ENDDO !  <><><> <><><> <><><> <><><>  <><><> <><><> <><><> <><><> <><><> <><><> <><><> <><><> <><><>

  return
  end subroutine bcs 

END MODULE boundariesmod

! Maybe I created a problem where there is none, and maybe if there is a problem I approached it the wrong way.
! It's my impression that communication between nodes might cause the program to "stall" if a node is trying to 
! send information to another node, but the other node is also trying to send information. The above do-loops 
! inside an if-block inside a do loop is designed to make sure that senders and receivers are unambigously mixed.
! In the first iteration of the main do-loop, the even workers send information to their neighbors. In the second
! iteration, the odd workers send information to their neighbors.
