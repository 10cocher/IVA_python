MODULE outdata
  !
  IMPLICIT NONE
  !
  include "dp.h"
  !
  PRIVATE
  PUBLIC outmatrix, outmatrix3
  !
CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE outmatrix(M,i,j,filename,comm)
    IMPLICIT NONE
    !Parameters
    INTEGER                ,INTENT(in)  :: i,j,comm
    REAL(4),DIMENSION(i,j),INTENT(in)  :: M
    CHARACTER(len=*)       ,INTENT(in)  :: filename
    !Local variables
    INTEGER :: ios,urang,rang
#ifdef do_mpi
    INTEGER :: code
#endif
    !
#ifdef do_mpi
    CALL MPI_COMM_RANK(comm,rang,code)
#else
    rang = 0 + 0*comm
#endif
    urang=20+rang
    !
    OPEN(UNIT=urang,IOSTAT=ios,FILE="out/"//filename//".dat",&
         STATUS="replace",ACCESS="stream",FORM="unformatted",&
         POSITION="rewind",ACTION="write")
    IF (ios.NE.0) WRITE(UNIT=*,FMT="(A,A)")"Pb open2 ",filename
    WRITE(UNIT=urang,IOSTAT=ios) M
    IF (ios.NE.0) WRITE(UNIT=*,FMT="(A,A)")"Pb write2 ",filename
    CLOSE(UNIT=urang,IOSTAT=ios)
    IF (ios.NE.0) WRITE(UNIT=*,FMT="(A,A)")"Pb close2 ",filename
  END SUBROUTINE outmatrix
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE outmatrix3(M,nz,nx,nh,filename,comm)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz,nx,nh),INTENT(in)  :: M
    INTEGER                     ,INTENT(in)  :: nz,nx,nh,comm
    CHARACTER(len=*)            ,INTENT(in)  :: filename
    !Local variables
    INTEGER :: ios,urang,rang
#ifdef do_mpi
    INTEGER :: code
#endif
    !
#ifdef do_mpi
    CALL MPI_COMM_RANK(comm,rang,code)
#else
    rang = 0 + 0*comm
#endif
    urang = 20 + rang
    !
    OPEN(UNIT=urang,IOSTAT=ios,FILE="out/"//filename//".dat",&
         STATUS="replace",ACCESS="stream",FORM="unformatted",&
         POSITION="rewind",ACTION="write")
    IF (ios.NE.0) WRITE(UNIT=*,FMT="(A,A)"), "/!\ Pb open3 ",filename
    WRITE(UNIT=urang,IOSTAT=ios) M
    IF (ios.NE.0) WRITE(UNIT=*,FMT="(A,A)"), "/!\ Pb write3 ",filename
    CLOSE(UNIT=urang,IOSTAT=ios)
    IF (ios.NE.0) WRITE(UNIT=*,FMT="(A,A)"), "/!\ Pb close3 ",filename
  END SUBROUTINE outmatrix3
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE outdata
