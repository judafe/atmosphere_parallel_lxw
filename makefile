#---------------------- Makefile for parallel -------------------------#
	##usage:  make compile ; make run  or  make pbs
#
#Nprocs=5 #Nprocs=Nworkers+1 must match data file (line 1)
#
#----------------------------------------------------------------------#
#============================= set MPI, compiler ======================#
# In Vega first, load appropriate modules for a C/Fortran compiler and 
# one of the mpi-wrappers (openmpi/mpich/mvapich2 etc)   
# Check: module list, module avail 
# 1. GNU-MPI (load one of following set) 
#    module load mvapich2/gcc/64/2.2rc1  
#    module load mpich/ge/gcc/64/3.2
#    module load openmpi/gcc/64/1.10.3  
# 2. INTEL-MPI (load one of following set) 
#    module load intel/compiler/64/2017/17.0.2 intel/mpi/64/2017/2.174 
#    module load intel/compiler/64/2017/17.0.2 mvapich2/intel/64/2.2rc1 
#    module load intel/compiler/64/2017/17.0.2 openmpi/intel/64/1.10.3 
# 3. CRAY-INTEL-MPI (load all of the following set)  
#    module load craype-broadwell    
#    module load PrgEnv-cray/1.0.2   
#    module unload mvapich2_cce/2.2rc1.0.3   
#    module load craype-network-infiniband
#    module load intel/compiler/64/2017/17.0.2
#    module load intel/mpi/64/2017/2.174
#    module load cray-impi/1.1.4  
##-----> set appropriate compiler_wrapper: mpif77 mpif90 mpicc mpic++
  COMP   = mpifort
##-----> set appropriate extension: f90 f  c  cpp
  EXT    = f90
  LFLAGs = 
#for C:  LFLAGs = -lm
##-------------------------- for all:
#  FLAGs  =  -fopenmp 
#  FLAGs  = -g  
# FLAGs  = -O3 
# FLAGs  = -fast 
#=========================== set source code  =========================#
##--------------->set names for your PROGram and std I/O files: 
  PROG = lxw2s.x 
  INPUT  = ./dat 
  OUTPUT = o.out 
# ##--------------------> set code components: List of object files
CODE_o = domainparmod.o initparamod.o setup.o numetsmod.o writeprog.o\
	mntwavparamod.o boundariesmod.o communications.o \
	integratemod.o main.o
#======================= create executable: make compile ============# 
compile: $(CODE_o) 
	$(COMP) $(FLAGs)  $(CODE_o)  -o $(PROG)  $(LFLAGs)
	@echo " >>> compiled on `hostname -s` with  $(COMP) <<<" 
	rm -f *.o 

$(CODE_o):%.o: %.$(EXT)
	$(COMP) $(FLAGs) -c $< -o $@
#======================= execute: make run | make pbs ==========# 
run: 
	mpirun  -np 2  $(PROG)  


gplot:
	gnuplot> load 'lab8.plt'
	rm -f load


pbs: 
	@ vi PBS_script 
	qsub PBS_script 

clean:
	rm -f *.o *~ *.mod

wipe:
	rm -f *.out out.* fort.* *.x *.err core* *.o *.mod

