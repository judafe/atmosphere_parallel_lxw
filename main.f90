!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
!
!  main.f90
!=========================================================================================
! This script initiates mpi, assigns each processor an identification number 
! ("myID") and initiates the two main subroutines of the program, MASTER_PROCESS_MPI(...)
! and WORKER_PROCESS_MPI(...). The Fortran95 function system_clock(...) is
! used to time the duration of the program.
!
!
!=========================================================================================
!
PROGRAM main

  use mpi
  use integratemod
  use domainparmod
  use writeprog


  implicit none

  integer:: nworkers,lastnode,nprocs,ncxs,myid,ierr
  integer(8)::  ti, tf, clock_max, clock_rate
  character (LEN=12):: charID



!! Initialize MPI:
   call MPI_INIT(ierr) 

!! Number rocesses in application:
   call MPI_COMM_SIZE( MPI_COMM_WORLD, Nprocs, ierr )


!! Assign 'myID' to each processor 0 thru Nprocs-1:
   call MPI_COMM_RANK(MPI_COMM_WORLD, myID, ierr)


!! Record Start time of Simulation
   if(myID.eq.0) then
	call system_clock(ti, clock_rate, clock_max)
   endif

!! Calculate the Number of x-Cells Per Integration Processor
   ncxs = (ncx-4)/nprocs+4 

   IF (myID.EQ.0) THEN
	call initialDisplay(1,nprocs)
 	if(nprocs.EQ.1)then
          	write(*,*)'There is only one active processor.'
          	write(*,*)'Try again in Serial'
		goto 999
	endif
   ENDIF
   write(*,*) 'entering integration process',myID
   call INTG_PROCESS_MPI(myID,nprocs,ncxs) 
   write(*,*) 'Integration Processor Finished:',myID,' : ierr = ', ierr
   	
   write(*,*) '==============================================================================',myID

!! Record End time of Simulation
   if(myID.eq.0) then !
	call sleep(1) !! sleep statement is used to organize the terminal output
	call system_clock(tf, clock_rate, clock_max)
	write(*,*) 'Elapsed wall-clock time: ',tf-ti
   endif

999 CONTINUE

!! Finalize MPI
    call MPI_FINALIZE( ierr ) 


END PROGRAM main
