MODULE integratemod
!
!
! ! The integration process was copied from Dr. J.Snively's Lax
!

contains

   SUBROUTINE INTG_PROCESS_MPI(myID,nprocs,ncxs)

   use mpi
   use setupmod

   use initparamod ! initial parameter mod
   use domainparmod ! domain parameter mod
   use nummetmod ! numerical methods mod
   use boundariesmod   
   use communications
   use writeprog

   implicit none
   integer,intent(IN)::myID,nprocs,ncxs
   double precision,dimension(ncy,ncxs):: onesm
   double precision,dimension(ncy,ncxs):: xmesh,ymesh,wind
   double precision,dimension(ncy,ncxs):: P0,rho0,kinvisc,tdiffus,pressure,velU,velV
   double precision,dimension(ncy-1,ncxs-1)::gravitym ! gravity mesh
   double precision,dimension(ncy,ncxs,4):: Q
   double precision,dimension(ncy,ncxs,2):: rhotemp,kinvite
   double precision,dimension(ncy-1,ncxs-1,4):: Qh ! used for half step
   double precision,dimension(ncy,ncxs,4)::flux
   double precision,dimension(ncy-1,ncxs-1,4)::fluxh,Sh
   double precision,dimension(ncy-2,ncxs-2,4)::S
   double precision::time,dt
   integer::iter,nsteps,nframe,ierr

   integer,parameter:: snapshotframe = 0, snapshotID = 0
   double precision:: difmax,difCFL
   integer:: ndtr,dtiter,dtr
   integer(8)::  tir, tfr, clock_max, clock_rate

!! Open Output Files
   call openfiles(nprocs,myID)

!! Create X and Y coordinate mesh
   call mesh(xmesh,ymesh,myID,ncxs)
   onesm = ymesh/ymesh 

!! Initial Pressure Mesh ====
   P0   = prssures*exp(-ymesh/H)
!! Initial Density Mesh ====
   rho0 = densitys*exp(-ymesh/H)

!! Gravity Mesh ==<>===<>==<>===<>==
   gravitym(:,1) = (P0(2:ncy,1)-P0(1:ncy-1,1))/&
	(-0.5*dy*(rho0(2:ncy,1)+rho0(1:ncy-1,1)))
   do iter = 2,ncxs-1
   	gravitym(:,iter) = gravitym(:,1)
   enddo !      ==<>===<>==<>===<>==

!! Generate Initial Wind Mesh. 'windtype' allows user to pick between
!!        several different wind profiles. See domainparmod.f90
   call windprofile(windtype,ymesh,ncxs, wind) 

!! Generate Initial Pressure Mesh. 'pressuretype' allows user to pick from
!!        several different pressure profiles. See domainparmod.f90
   call pressureprofile(pressuretype,xmesh,ymesh,ncxs, pressure)

!! Initiate Conservation Variables. Used in integration
   Q(:,:,1) = rho0
   Q(:,:,2) = rho0*wind
   Q(:,:,3) =  0.d0*onesm
   Q(:,:,4) = (pressure+P0)/(gamm-1)+(0.5)*rho0*(wind**2) 

 !! write Q to file.========= WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ
	!if(mod(nsteps,rrrate).eq.0) then ! if recording condtion is met
		!if(myID.eq.0) call system_clock(tir, clock_rate, clock_max)
		!do iter = 1,4
		!	call writetofile(Q(3:ncy-2,3:ncxs-2,iter), myID+(iter*100),ncxs)
		!enddo
			!call writetofile(ymesh(3:ncy-2,3:ncxs-2), myID+(4*100),ncxs)
	!elseif(myID.eq.0) then
		!write(*,*) nsteps, time
	!endif
  !!  WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ

!! Initiate Conservation Variables used in Half Step
   pressure = 0.d0*onesm(2:ncy,2:ncxs)
   Qh(:,:,1)   = 0.d0*onesm(2:ncy,2:ncxs)
   Qh(:,:,2)  = 0.d0*onesm(2:ncy,2:ncxs)
   Qh(:,:,3)  = 0.d0*onesm(2:ncy,2:ncxs)
   Qh(:,:,4) = 0.d0*onesm(2:ncy,2:ncxs)    



   call bcb(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),P0,rho0,wind,ncxs) ! closed bottom
   call bct(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),P0,rho0,wind,ncxs) ! open top
   call bcs(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),ncxs,nprocs,myID) ! bc, periodic sides


! ==============Mesh Viscosity Coefficients==================================!
   kinvisc = (1.3e-5*ymesh/ymesh)/rho0
   kinvite(:,:,1) = kinvisc; kinvite(:,:,2) = kinvisc
   tdiffus = kinvisc/prandtl

! === Initial Time Step ================================================!
  dt = dCFL*min(dx,dy)/csos ! changes within loop
  nsteps=0
  time = 0
  nframe = 1;

  if(myID.eq.0) then ! write to terminal
	write(*,*) '======================================================================'
	write(*,*) '   Step No.       Simulated Time     Frame Number       Recording time'
  endif   ! end write to terminal

  DO WHILE(time.le.tmax.and.nsteps.le.maxsteps) ! ~ TIME INTEGRATION LOOP ~ TIME INTEGRATION LOOP ~

	call xflux(Q,ncy,ncxs, flux)
	Qh(2:ncy-1,2:ncxs-1,:) = &
		(0.5)*(Q(2:ncy-1,2:ncxs-1,:)+Q(2:ncy-1,3:ncxs,:))&
		-(0.5*dt/dx)*(flux(2:ncy-1,3:ncxs,:)-flux(2:ncy-1,2:ncxs-1,:))

             
	call xflux(Qh,ncy-1,ncxs-1, fluxh)
	Q(2:ncy-1,2:ncxs-1,:) = &
		Q(2:ncy-1,2:ncxs-1,:) - &
		(dt/dx)*(fluxh(2:ncy-1,2:ncxs-1,:)-fluxh(2:ncy-1,1:ncxs-2,:))


	!! Boundary Conditions
	call bcb(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),P0,rho0,wind,ncxs) ! boundary condition
	call bct(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),P0,rho0,wind,ncxs) ! bc, open top
	call bcs(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),ncxs,nprocs,myID) ! bc, periodic sides

        !! Y Split ============================================================
	call yflux(Q,ncy,ncxs, flux)
	  ! gravity "source"; half-step 1
	call source((0.5*Q(2:ncy-1,2:ncxs-1,:)+0.5*Q(3:ncy,2:ncxs-1,:)),gravitym(2:ncy-1,2:ncxs-1),ncy-2,ncxs-2,S)

	Qh(2:ncy-1,2:ncxs-1,:) = &
  		 (0.5)*(Q(2:ncy-1,2:ncxs-1,:)+Q(3:ncy,2:ncxs-1,:)) &
		-(0.5*dt/dy)*(flux(3:ncy,2:ncxs-1,:)-flux(2:ncy-1,2:ncxs-1,:)) &
		+(0.5*dt)*S


	call yflux(Qh,ncy-1,ncxs-1, fluxh)
	  ! gravity "source"; half-step 2	
	call source(Qh,gravitym,ncy-1,ncxs-1,Sh)
	Q(2:ncy-1,2:ncxs-1,:) = &
		Q(2:ncy-1,2:ncxs-1,:) &
		-(dt/dy)*(fluxh(2:ncy-1,2:ncxs-1,:)-fluxh(1:ncy-2,2:ncxs-1,:)) &
		+ 0.5*dt*(Sh(2:ncy-1,2:ncxs-1,:)+Sh(1:ncy-2,2:ncxs-1,:))

	!! Boundary Conditions
	call bcb(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),P0,rho0,wind,ncxs) ! boundary condition, closed bottom
	call bct(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),P0,rho0,wind,ncxs) ! bc, open top
	call bcs(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),ncxs,nprocs,myID) ! bc, periodic sides


  !! --Viscosity         --~---~--~---~--~---~--~---~--~---~--~---~--~---~--~---~--~---~--~---~--~---~--~---~
	IF(viscosity01.eq.1) then
	   !write(*,*) 'viscosity is on'
	   difmax = maxval(kinvisc)
	   difCFL = 0.25
	   ndtr = max(dt*difmax/(dxymin**2)/difCFL,1.d0)


		DO dtiter=1,ndtr
		dtr = dt/ndtr
		Q(:,:,4) = Q(:,:,4)-0.5*(Q(:,:,2)**2+Q(:,:,3)**2)/Q(:,:,1)

		! calculate velocities
		rhotemp(:,:,1) = Q(:,:,1); rhotemp(:,:,2) = Q(:,:,1) 
		Q(:,:,2:3) = Q(:,:,2:3)/rhotemp

		Q(2:ncy-1,2:ncxs-1,2:3) = rhotemp(2:ncy-1,2:ncxs-1,:) &
			*(Q(2:ncy-1,2:ncxs-1,2:3)+kinvite(2:ncy-1,2:ncxs-1,:)*(dtr/dx**2) &
			*(Q(2:ncy-1,3:ncxs,2:3)-2*Q(2:ncy-1,2:ncxs-1,2:3)+Q(2:ncy,1:ncxs-2,2:3)))
		
		! re-calculating energy
		Q(2:ncy-1,2:ncxs-1,4) = Q(2:ncy-1,2:ncxs-1,4) &
			+ 0.5*(Q(2:ncy-1,2:ncxs-1,2)**2+Q(2:ncy-1,2:ncxs-1,3)**2)/Q(2:ncy-1,2:ncxs-1,1)
		!! Boundary Conditions
		call bcb(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),P0,rho0,wind,ncxs) ! boundary condition, closed bottom
		call bct(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),P0,rho0,wind,ncxs) ! bc, open top
		call bcs(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),ncxs,nprocs,myID) ! bc, periodic sides

		Q(:,:,4) = Q(:,:,4)-0.5*(Q(:,:,2)**2+Q(:,:,3)**2)/Q(:,:,1)

		! calculating velocity
		rhotemp(:,:,1) = Q(:,:,1); rhotemp(:,:,2) = Q(:,:,1) 
		Q(:,:,2:3) = Q(:,:,2:3)/rhotemp
		Q(2:ncy-1,2:ncxs-1,2:3) = rhotemp(2:ncy-1,2:ncxs-1,:) &
			*(Q(2:ncy-1,2:ncxs-1,2:3)+kinvite(2:ncy-1,2:ncxs-1,:)*(dtr/dy**2) &
			*(Q(3:ncy,2:ncxs-1,2:3)-2*Q(2:ncy-1,2:ncxs-1,2:3)+Q(2:ncy,1:ncxs-2,2:3)))

		! recalculating energy
		Q(2:ncy-1,2:ncxs-1,4) = Q(2:ncy-1,2:ncxs-1,4) &
			+ 0.5*(Q(2:ncy-1,2:ncxs-1,2)**2+Q(2:ncy-1,2:ncxs-1,3)**2)/Q(2:ncy-1,2:ncxs-1,1)
		!! Boundary Conditions
		call bcb(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),P0,rho0,wind,ncxs) ! boundary condition, closed bottom
		call bct(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),P0,rho0,wind,ncxs) ! bc, open top
		call bcs(Q(:,:,1),Q(:,:,2),Q(:,:,3),Q(:,:,4),ncxs,nprocs,myID) ! bc, periodic sides

		ENDDO
	ENDIF ! viscosity01 		
  !! --End Viscosity      --~---~--~---~--~---~--~---~--~---~--~---~--~---~--~---~--~---~--~---~--~---~--~---~

  !! Mountain Wave Source Function -----MWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMW
	if(mntsource01.eq.1) then
  		call sourcemnt(ymesh,xmesh,ncxs,time,dt, Q)
	endif
  !!--End Mountain Wave Source Function -----MWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWMWM


 !! write Q to file.========= WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ
	if(mod(nsteps,rrrate).eq.0) then ! if recording condtion is met
		if(myID.eq.0) call system_clock(tir, clock_rate, clock_max)
		do iter = 1,4
			call writetofile(Q(3:ncy-2,3:ncxs-2,iter), myID+(iter*100),ncxs)
		enddo
			!call writetofile(ymesh(3:ncy-2,3:ncxs-2), myID+(4*100),ncxs)
		if(myID.eq.0) then
			call system_clock(tfr, clock_rate, clock_max)
			write(*,*)  nsteps, time, nframe,tfr-tir
			nframe = nframe+1
		endif
	elseif(myID.eq.0) then
		!write(*,*) nsteps, time
	endif
  !!  WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ WQ


	nsteps = nsteps+1
	time = nsteps*dt


  END DO ! ~ TIME INTEGRATION LOOP ~ TIME INTEGRATION LOOP ~ TIME INTEGRATION LOOP ~
  


  RETURN

  END SUBROUTINE INTG_PROCESS_MPI

END MODULE  integratemod


	! Used primarily for troubleshooting. May be used to take a "snap shot" of Q at any chosen 
	! step. Snapshot parameters at the top of this mod.
  !	if(myID.eq.snapshotID.and.nsteps.eq.snapshotframe) then
 !		open (unit = myID+10, file = "rho_ss.csv")   
 !		open (unit = myID+20, file = "rhoU_ss.csv")   
 !		open (unit = myID+30, file = "rhoV_ss.csv")   
 !		open (unit = myID+40, file = "energy_ss.csv")   
!		do iter = 1,ncy
!			write(myID+10,*) rho(iter,:)
!			write(myID+20,*) rhoU(iter,:)
!			write(myID+30,*) rhoV(iter,:)
!			write(myID+40,*) energy(iter,:)
!		enddo
!		write(*,*) myID,' finished writing'
 ! 		close(myID+10)
 !  		close(myID+20)
!   		close(myID+30)
!  		close(myID+40)  
!  	endif 
