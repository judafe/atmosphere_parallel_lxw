! domainparmod.f90

MODULE domainparmod

  implicit none

   ! RESTRICTION: mod(xdomain/dx,nworkers) = 0
   double precision,parameter:: Xmin = -100000, Xmax = 100000
   double precision,parameter:: Ymin = 0, Ymax = 160000
   double precision,parameter:: dx = 200, dy = 200
   double precision,parameter:: dxymin = min(dx,dy)
   double precision,parameter:: xdomain = Xmax - Xmin
   double precision,parameter:: ydomain = Ymax - Ymin
   integer,parameter:: ncx = xdomain/dx+4 ! number of cell centers (I in matlab code)
   integer,parameter:: ncy = ydomain/dy+4 ! number of cell centers (J in matlab code)

   double precision,parameter:: tmin = 0, tmax = 36000!2*60*60 ! 
   integer,parameter:: maxsteps = 1000000
   integer,parameter:: rrrate = 40 ! reciprocal record rate

   double precision,parameter:: dCFL = 0.8d0 ! desired CFL number

   double precision,parameter:: dynvisc = 1.3e-5 ! dynamic viscosity
   double precision,parameter::prandtl = 0.7d0 ! Prandtl number



END MODULE domainparmod


