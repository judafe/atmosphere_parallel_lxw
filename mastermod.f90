!  master.f90 contains subroutine master_process_mpi(Nprocs)
!----------------------------------------------------------------------!
!      
!----------------------------------------------------------------------!
MODULE mastermod

  implicit none

CONTAINS

  SUBROUTINE master_process_mpi(nprocs,myID)

  use mpi
  use setupmod

  use initparamod ! initial parameter mod
  use domainparmod ! domain parameter mod
  use mntwavparamod ! mountain wave parameter mod
  use communications
  use writeprog
  implicit none
  
  ! Parallelization Variables 
  	integer,intent(in)::nprocs,myID 
  	integer::nworkers,ncxs,nsteps,recturn,ierr,lastnode
 	integer::status(MPI_STATUS_SIZE)	
  ! Spatial and Temporal Domain variables ===============================
	double precision,dimension(ncy,ncx)::wind
	double precision,dimension(ncy,ncx)::xmesh,ymesh
   	double precision,dimension(ncy,ncx):: pressure,rho,rhoU,rhoV,energy
	! number of x-cell centers per worker	
  	double precision:: dt,time

	double precision,dimension(:,:),allocatable::kinvisc,tdiffus

	integer::iter,jter
	double precision, dimension(:,:),allocatable::onesm
 	integer(8)::  t01,t02,t03,t04,clock_max, clock_rate

	nworkers = nprocs-2
	lastnode = nprocs-1

  	open (unit = 10, file = "pp_rhoA.txt")   
 	open (unit = 20, file = "pp_rhoUA.txt")   
 	open (unit = 30, file = "pp_rhoVA.txt")   
	open (unit = 40, file = "pp_energyA.txt")  
  	open (unit = 11, file = "pp_rhoB.txt")   
 	open (unit = 22, file = "pp_rhoUB.txt")   
 	open (unit = 33, file = "pp_rhoVB.txt")   
	open (unit = 44, file = "pp_energyB.txt")  

   	ncxs = (ncx-4)/nworkers+4
   	dt = dCFL*min(dx,dy)/csos   
	recturn = 2 ! "recording turn", used to track which processor is receiving Q
	nsteps = 0
  	DO WHILE(time.le.tmax.and.nsteps.le.maxsteps) ! ~ TIME INTEGRATION LOOP ~ TIME INTEGRATION LOOP ~.
		IF(mod(nsteps,rrrate).eq.0) THEN ! recording frame rate
			if(myID.eq.0.and.mod(recturn,2).eq.0) then
				call system_clock(t01, clock_rate, clock_max)
				call recvQ(nsteps,time,rho,rhoU,rhoV,energy,ncxs)
				call writetofile(rho,10)
				call writetofile(rhoU,20)
				call writetofile(rhoV,30)
				call writetofile(energy,40)
				write(50,*) time
				call system_clock(t02, clock_rate, clock_max)
				write(*,*) myID,nsteps,time,t02-t01
				!write(*,*) '============================================================================================'
 				!write(*,*) 'Recording Node',myID,'Reocorded Grids at step number:',nsteps,' at simulation Time:',time
				!write(*,*) 'Receive and Record Proces, Elapsed Wallclock time:',t02-t01
				!write(*,*) '============================================================================================'
			elseif(myID.eq.lastnode.and.mod(recturn,2).ne.0) then
				call system_clock(t03, clock_rate, clock_max)
				call  recvQ(nsteps,time,rho,rhoU,rhoV,energy,ncxs)
				call writetofile(rho,11)
				call writetofile(rhoU,22)
				call writetofile(rhoV,33)
				call writetofile(energy,44)
				write(55,*) time
				call system_clock(t04, clock_rate, clock_max)
				write(*,*) myID,nsteps,time,t04-t03
				!write(*,*) '============================================================================================'
 				!write(*,*) 'Recording Node',myID,'Reocorded Grids at step number:',nsteps,' at simulation Time:',time
				!write(*,*) 'Receive and Record Proces, Elapsed Wallclock time:',t04-t03
				!write(*,*) '============================================================================================'
			endif
			recturn = recturn+1
		ENDIF

		nsteps = nsteps+1
		time = nsteps*dt

  	END DO ! ~ TIME INTEGRATION LOOP ~ TIME INTEGRATION LOOP ~ TIME INTEGRATION LOOP ~
  call sleep(1) ! perhaps uneccesary. Used to help organize the terminal/command prompt displayy



   RETURN

   END  SUBROUTINE master_process_mpi

END MODULE mastermod
