MODULE foo
  !
  IMPLICIT NONE
  !
  include "dp.h"
  include "param.h"
#ifdef do_mpi
  ! Loads mpi module
  include 'mpif.h'
  include 'dp_mpi.h'
#endif
  !
  PRIVATE
  PUBLIC zadd, zaddn, test_mpi
CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE test_mpi(a,b,comm,nz,nx)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz,nx),INTENT(out) :: a
    REAL(4),DIMENSION(nz,nx),INTENT(in)  :: b
    INTEGER, INTENT(in) :: comm,nz,nx
    !Local variables
    INTEGER :: rang,code
    !
    code = 17
    !
#ifdef do_mpi
    !INTEGER :: code
    CALL MPI_COMM_RANK(comm,rang,code)
#else
    rang=0
#endif
    !
    WRITE(UNIT=*,FMT="(A,I10,A,I10,A,I10,A,F5.3)")&
         "rang=",rang," comm=",comm, "code=",code," tab(1,1)=",b(1,1)
    !
#ifdef do_mpi
    CALL MPI_ALLREDUCE(b,a,nz*nx,MPI_OWNREAL,MPI_SUM,comm,code)
#else
    a = b
#endif
    !
  END SUBROUTINE test_mpi
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE zaddn(s,a,b,nz,nx)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz,nx),INTENT(out) :: s
    REAL(4),DIMENSION(nz,nx),INTENT(in)  :: a,b
    INTEGER                 ,INTENT(in)  :: nz,nx
    !
    s = a + b
  END SUBROUTINE zaddn
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE zadd(s,a,b)
    IMPLICIT NONE
    !Parameters
    REAL(4),INTENT(out) :: s
    REAL(4),INTENT(in)  :: a,b
    !
    s = a + b
  END SUBROUTINE zadd
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE foo
