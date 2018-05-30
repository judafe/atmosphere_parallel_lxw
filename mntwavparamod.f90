!mntwavparamod.f90

MODULE mntwavparamod

  use domainparmod
  implicit none

   double precision,parameter:: amplitude = 0.34 ! wave amplitude [m/s^2]     
   double precision,parameter:: sigy = 3000, sigx = sqrt(0.6*xdomain) ! wave sigma [m]
   double precision,parameter:: sigt = 1.5*60.d0*60.d0 ! wave sigma [s]
   double precision,parameter:: xw_cent = 0, yw_cent = 10000 ! spatial wave center
   double precision,parameter:: tw_cent = 6.d0*60*60 ! temporal wave center


END MODULE mntwavparamod

