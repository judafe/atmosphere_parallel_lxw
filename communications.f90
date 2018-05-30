! contains subroutines 
!----------------------------------------------------------------------!
MODULE communications

  implicit none

CONTAINS

  SUBROUTINE sendQ(dest,nsteps,time,rho,rhoU,rhoV,energy,ncxs,myID)
  use domainparmod
  use mpi

  implicit none  

  integer,intent(IN)::myID,ncxs,nsteps,dest
  double precision,intent(IN)::time
  double precision,dimension(ncy,ncxs),intent(IN)::rho,rhoU,rhoV,energy
  integer::nworkers,start,fin,iter,ierr
  nworkers=(ncx-4)/(ncxs-4)


  if(myID.eq.1) then
     start = 1
     fin = ncxs-2
     call MPI_SEND(nsteps, 1, MPI_INTEGER,dest,myID+100,MPI_COMM_WORLD,ierr)
     call MPI_SEND(time, 1, MPI_DOUBLE_PRECISION,dest,myID+200,MPI_COMM_WORLD,ierr)
  elseif(myID.eq.nworkers) then
     start = 3
     fin = ncxs
  else
     start = 3
     fin = ncxs-2	
  endif
  do iter = start,fin
	call MPI_SEND(rho(:,iter),  ncy, MPI_DOUBLE_PRECISION,dest,myID+1000,MPI_COMM_WORLD,ierr)
	call MPI_SEND(rhoU(:,iter), ncy, MPI_DOUBLE_PRECISION,dest,myID+2000,MPI_COMM_WORLD,ierr)
	call MPI_SEND(rhoV(:,iter), ncy, MPI_DOUBLE_PRECISION,dest,myID+3000,MPI_COMM_WORLD,ierr)
	call MPI_SEND(energy(:,iter),ncy,MPI_DOUBLE_PRECISION,dest,myID+4000,MPI_COMM_WORLD,ierr)
  enddo ! iter loop

  RETURN

  END SUBROUTINE

! ===<>==== ===<>==== ===<>==== ===<>==== ===<>==== ===<>==== ===<>====
  SUBROUTINE  recvQ(nsteps,time,rho,rhoU,rhoV,energy,ncxs)
  use domainparmod
  use mpi

  implicit none  

  double precision,dimension(ncy,ncx),intent(OUT)::rho,rhoU,rhoV,energy
  double precision,intent(IN)::time
  integer,intent(IN)::ncxs,nsteps
  integer::nworkers,sndr,start,fin,iter,ierr
  nworkers=(ncx-4)/(ncxs-4)

  DO sndr=1,nworkers
     if(sndr.eq.1) then
	start = 1
        fin = ncxs-2
        call MPI_RECV(nsteps,1,MPI_INTEGER,sndr,sndr+100,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        call MPI_RECV(time,ncy,MPI_DOUBLE_PRECISION,sndr,sndr+200,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
     elseif(sndr.eq.nworkers) then
        start = 3+(ncxs-4)*(sndr-1)
	fin = ncx	
     else
        start = 3+(ncxs-4)*(sndr-1)
     	fin = start+ncxs-5	
     endif
     !write(*,*) 'recv range',sndr,fin-start
     do iter = start,fin
           call MPI_RECV(rho(:,iter),ncy,MPI_DOUBLE_PRECISION,sndr,sndr+1000,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
           call MPI_RECV(rhoU(:,iter),ncy,MPI_DOUBLE_PRECISION,sndr,sndr+2000,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
           call MPI_RECV(rhoV(:,iter),ncy,MPI_DOUBLE_PRECISION,sndr,sndr+3000,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
           call MPI_RECV(energy(:,iter),ncy,MPI_DOUBLE_PRECISION,sndr,sndr+4000,&
		MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
     enddo
  ENDDO

  END SUBROUTINE 

END MODULE communications











