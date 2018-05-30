! contains subroutines xflux() and yflux
!----------------------------------------------------------------------!
MODULE writeprog

  implicit none

CONTAINS
  SUBROUTINE initialDisplay(onoff,nworkers)
	use domainparmod
	implicit none
	integer,intent(IN)::onoff,nworkers
	if(onoff.eq.1) then
		write(*,*)	
		write(*,*) '=========================================================================='
        	write(*,*) 'Initiating Dimensionally Split Lax-Wendroff 2-Step Integrator'
		write(*,*)
		write(*,*) 'Number of processors active: ',nworkers
		write(*,*) 'Time Simulated (s):',tmax-tmin
		write(*,*)	
		write(*,*) 'X Domain (m):',xdomain
		write(*,*) 'Y Domain (m):',ydomain
		write(*,*) 'dx (m):',dx
		write(*,*) 'dy (m):',dy
		write(*,*) 'total x-cells',ncx
		write(*,*) 'total y-cells',ncy	
		write(*,*)	
		write(*,*) 'x-cells per worker',(ncx-4)/nworkers+4
   		write(*,*)
   		write(*,*) 'Recording Rate',rrrate
		write(*,*) '=========================================================================='
   		write(*,*)
	endif
  RETURN

  END SUBROUTINE

!!! INITIATE TEXT FILES

  SUBROUTINE openfiles(nprocs,myID)
	
	implicit none
	integer,intent(IN)::nprocs,myID
	character(LEN=12)::charID,rho,rhoU,rhoV,energy
	
	if(nprocs.lt.10) then
	   	write(charID,'(I1)') myID
	elseif(nprocs.lt.100) then
		write(charID,'(I2)') myID
	else
		write(*,*) 'From writeprog.f90:'
		write(*,*) 'STOP KIDDING YOURSELF'
	endif

  	write(rho ,*) trim('rho'//trim(charID)//'.txt')
  	write(rhoU,*) trim('rhoU'//trim(charID)//'.txt')
  	write(rhoV,*) trim('rhoV'//trim(charID)//'.txt')
	write(energy,*) trim('energy'//trim(charID)//'.txt')

	open (unit = myID+100, file = rho)   
 	open (unit = myID+200, file = rhoU)   
 	open (unit = myID+300, file = rhoV)   
 	open (unit = myID+400, file = energy)   
	
  RETURN
  END SUBROUTINE

!!! WRITE TO FILE

  SUBROUTINE writetofile(mesh,fileID,ncxs)
	use domainparmod

	implicit none
	integer,intent(IN)::fileID,ncxs
	double precision,dimension(ncy-4,ncxs-4)::mesh(1:ncy-4,1:ncxs-4)
	integer::iter

	do iter = 1,(ncy-4)	
		write(fileID,*) mesh(iter,:)
	enddo
  RETURN

  END SUBROUTINE

END MODULE writeprog
