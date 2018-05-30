!initparamod.f90

MODULE initparamod

  implicit none


   double precision,parameter:: gravitys = 9.8d0 ! gravity scalar     
   double precision,parameter:: gasconst = 287.d0 ! gass constant 'R'
   double precision,parameter:: prssures = 1.0e5 ! pressure scalar
   double precision,parameter:: densitys = 1.2 ! density scalar
   double precision,parameter:: gamm = 1.4d0 ! gamma ratio
   double precision,parameter:: csos = sqrt(gamm*prssures/densitys)! speed of sound
   double precision,parameter:: H = prssures/(densitys*gravitys) ! scaled eight
 

   integer,parameter:: windtype = 1 ! wind profile command (see below)
   integer,parameter:: pressuretype = 0 ! pressure profile command (see below)
   integer,parameter:: viscosity01 = 1 ! viscosity switch
   integer,parameter:: mntsource01 = 1 ! mount. wave source switch
!==================================================
!--------  Scenario Switches ---------------------
!==================================================
!wind type
!  1 --> jet stream near 100 km altitude
!  2 --> sharp decrease near 100 km, 40 m/s to 0 m/s
!  3 --> sharp decrease near 100 km, 40 m/s to -40 m/s
!  4 --> constant wind at 40 m/s
!  else> no wind
!
!pressure type
!
!  1 --> gaussian pressure pulse
!  2 --> tanh vertical split
!  else> no pressure perturbation
!
! viscosity 0/1
!  1 --> turns off viscosity
  
   

END MODULE initparamod
