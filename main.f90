PROGRAM main
  !
  !use linop, only: Flin
  !use foo  , only: mpi_test
  !
  IMPLICIT NONE
  !
  !include "dp.h"
  !include "param.h"
  !
#ifdef do_mpi
  !Loads mpi module
  include 'mpif.h'
  !INTEGER, PARAMETER :: MPI_OWNREAL=MPI_REAL
#endif
  INTEGER :: rang, nb_procs
#ifdef do_mpi
  INTEGER :: code
#endif
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%---MPI---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef do_mpi
  ! MPI code
  CALL MPI_INIT(code)
  ! Determines number of processes
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
  ! Dtermines rank of process
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
#else
  ! No MPI code
  rang=0
  nb_procs=1
#endif
  !
#ifdef do_mpi
  print *, MPI_COMM_WORLD
#endif
  !
  WRITE(UNIT=*,FMT="(A,I2)") "Hello, I'm processor number", rang
  !
END PROGRAM main
