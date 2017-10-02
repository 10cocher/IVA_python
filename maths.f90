MODULE maths
  !
  PRIVATE
  PUBLIC SmoothSin, deftapZ, deftapX, deftapX2, deftapXdisym
  !
  include "dp.h"
  include "param.h"
  !
CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE REAL(4) FUNCTION SmoothSin(x)
    IMPLICIT NONE
    !
    REAL(4),INTENT(in) :: x
    !
    !REAL(4) :: pi
    !
    !pi=3.141592653589793_4
    !
    SmoothSin=( sin(pi*x/2._4) )**2
  END FUNCTION SmoothSin
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE deftapZ(tapZ,iz1,iz2,izb,nz)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz),INTENT(out) :: tapZ
    INTEGER               ,INTENT(in)  :: iz1,iz2,izb,nz
    !Local variables
    INTEGER :: iz
    !
    tapZ(:)=1._4
    tapZ(1:iz1)=0._4
    DO iz=iz1+1,iz2-1
       tapZ(iz)=SmoothSin(REAL(iz-iz1,4)/REAL(iz2-iz1-1,4))
    END DO
    DO iz=1,izb-1
       tapZ(nz-iz+1)=SmoothSin(REAL(iz-1,4)/REAL(izb-1,4))       
    END DO
  END SUBROUTINE deftapZ
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE deftapX(tapX,ixt,nxx)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nxx),INTENT(out) :: tapX
    INTEGER                ,INTENT(in)  :: ixt,nxx
    !Local variables
    REAL(4) :: val
    INTEGER :: ix
    !
    tapX(:)=1._4
    IF (ixt.GT.nxx) RETURN
    IF (nxx.GT.1) THEN
       DO ix=1,ixt-1
          val=SmoothSin(REAL(ix-1,4)/REAL(ixt-1,4))
          tapX(ix)=val
          tapX(nxx-ix+1)=val
       END DO
    END IF
    !
  END SUBROUTINE deftapX
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE deftapXdisym(tapX,ixt,ixt2,nxx)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nxx),INTENT(out) :: tapX
    INTEGER                ,INTENT(in)  :: ixt,ixt2,nxx
    !Local variables
    REAL(4) :: val
    INTEGER :: ix
    !
    tapX(:)=1._4
    IF (ixt.GT.nxx) RETURN
    IF (nxx.GT.1) THEN
       DO ix=1,ixt-1
          val=SmoothSin(REAL(ix-1,4)/REAL(ixt-1,4))
          tapX(ix)=val
          tapX(ixt2-ix+1)=val
       END DO
       DO ix=ixt2,nxx
          tapX(ix)=0._4
       END DO
    END IF
    !
  END SUBROUTINE deftapXdisym
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE deftapX2(tapX,ixt,nxx)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nxx),INTENT(out) :: tapX
    INTEGER                ,INTENT(in)  :: ixt,nxx
    !Local variables
    REAL(4) :: val
    INTEGER :: ix
    !
    tapX(:)=1._4
    DO ix=1,ixt
       val=SmoothSin(REAL(ix,4)/REAL(ixt,4))
       tapX(ix)=val
       tapX(nxx-ix+1)=val
    END DO
  END SUBROUTINE deftapX2
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE maths
