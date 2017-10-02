MODULE intderiv
  !
  !
  IMPLICIT NONE
  !
  !include "dp.h"
  !include "param.h"
  !
  !dp = 4
  !
  PRIVATE
  PUBLIC dev2x,dev1x,int1x,int1xAdj,int2x,&
       fliptime,&
       int1x_1D,int1xAdj_1D,dev2x_1D,dev1x_1D,dev2x2_1D
  !
CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE dev2x(Q,P,dt,ntt,nz,noff)
    !Parameters
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(out):: Q
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(in) :: P
    REAL(4)                       ,INTENT(in) :: dt
    INTEGER                        ,INTENT(in) :: ntt,nz,noff
    !Local variables
    REAL(4) :: temp
    temp=1._4/(dt**2)
    Q(2:ntt-1,:,:)=(P(3:ntt,:,:)+P(1:ntt-2,:,:)-2._4*P(2:ntt-1,:,:))*temp
    Q(1  ,:,:)=Q(2    ,:,:)
    Q(ntt,:,:)=Q(ntt-1,:,:)
  END SUBROUTINE dev2x
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE dev1x(Q,P,dt,ntt,nz,noff)
    !Parameters
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(out):: Q
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(in) :: P
    REAL(4)                       ,INTENT(in) :: dt
    INTEGER                        ,INTENT(in) :: ntt,nz,noff
    !Local variables
    REAL(4) :: temp,temp2
    temp=1._4/(2._4*dt)
    temp2=1._4/(dt)
    Q(2:ntt-1,:,:)=(P(3:ntt,:,:)-P(1:ntt-2,:,:))*temp
    Q(1  ,:,:)=(P(2,:,:)-P(1,:,:))*temp2
    Q(ntt,:,:)=(P(ntt,:,:)-P(ntt-1,:,:))*temp2
  END SUBROUTINE dev1x
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE int1x(Q,P,dt,ntt,nz,noff)
    !Parameters
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(out) :: Q
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(in)  :: P
    REAL(4)                       ,INTENT(in)  :: dt
    INTEGER                        ,INTENT(in)  :: noff,nz,ntt
    !Local variables
    REAL(4) :: temp
    INTEGER  :: it,iz,ix
    !
    temp=dt/24._4

    DO ix=1,noff
       DO iz=1,nz
          Q(1,iz,ix)=0._4
          DO it=2,ntt-1
             Q(it,iz,ix)=Q(it-1,iz,ix)&
                  +temp*(P(it+1,iz,ix)+22._4*P(it,iz,ix)+P(it-1,iz,ix))
          END DO
          Q(ntt,iz,ix)=Q(ntt-1,iz,ix)&
               +temp*(23._4*P(ntt,iz,ix)+P(ntt-1,iz,ix))
       END DO
    END DO
  END SUBROUTINE int1x
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE int1xAdj(Q,P,dt,ntt,nz,noff)
    !Parameters
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(out) :: Q
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(in)  :: P
    REAL(4)                       ,INTENT(in)  :: dt
    INTEGER                        ,INTENT(in)  :: noff,nz,ntt
    !Local variables
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: S,R
    !
    ALLOCATE(S(ntt,nz,noff))
    S(:,:,:)=0._4
    ! Flip in time
    CALL FlipTime(S,P,ntt,nz,noff)
    ! Integerate
    ALLOCATE(R(ntt,nz,noff))
    CALL int1x(R,S,dt,ntt,nz,noff)
    ! Flip in time again
    DEALLOCATE(S)
    CALL FlipTime(Q,R,ntt,nz,noff)
    DEALLOCATE(R)
    !
  END SUBROUTINE int1xAdj
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE int2x(Q,P,dt,ntt,nz,noff)
    !Parameters
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(out) :: Q
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(in)  :: P
    REAL(4)                       ,INTENT(in)  :: dt
    INTEGER                        ,INTENT(in)  :: noff,nz,ntt
    !Local variables
    REAL(4),DIMENSION(ntt,nz,noff) :: temp
    !
    CALL int1x(temp,P,dt,ntt,nz,noff)
    CALL int1x(Q,temp,dt,ntt,nz,noff)
  END SUBROUTINE int2x
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE FlipTime(Q,P,ntt,nz,noff)
    !Parameters
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(out) :: Q
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(in)  :: P
    INTEGER                        ,INTENT(in)  :: noff,nz,ntt
    !Local variables
    INTEGER :: it
    DO it=1,ntt
       Q(it,:,:)=P(ntt+1-it,:,:)
    END DO
    !
  END SUBROUTINE FlipTime
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE dev2x_1D(Q,P,dt,ntt,nz)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz),INTENT(out):: Q
    REAL(4),DIMENSION(ntt,nz),INTENT(in) :: P
    REAL(4)              ,INTENT(in) :: dt
    INTEGER               ,INTENT(in) :: ntt,nz
    !Local variables
    REAL(4) :: temp
    temp=1._4/(dt**2)
    Q(2:ntt-1,:)=(P(3:ntt,:)+P(1:ntt-2,:)-2._4*P(2:ntt-1,:))*temp
    Q(1  ,:)=Q(2    ,:)
    Q(ntt,:)=Q(ntt-1,:)
  END SUBROUTINE dev2x_1D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE dev2x2_1D(Q,P,dz,ntt,nz)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz),INTENT(out):: Q
    REAL(4),DIMENSION(ntt,nz),INTENT(in) :: P
    REAL(4)              ,INTENT(in) :: dz
    INTEGER               ,INTENT(in) :: ntt,nz
    !Local variables
    REAL(4) :: temp
    temp=1._4/(dz**2)
    Q(:,2:nz-1)=(P(:,3:nz)+P(:,1:nz-2)-2._4*P(:,2:nz-1))*temp
    Q(:,1)=Q(:,2)
    Q(:,nz)=Q(:,nz-1)
  END SUBROUTINE dev2x2_1D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE dev1x_1D(Q,P,dt,ntt,nz)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz),INTENT(out):: Q
    REAL(4),DIMENSION(ntt,nz),INTENT(in) :: P
    REAL(4)              ,INTENT(in) :: dt
    INTEGER               ,INTENT(in) :: ntt,nz
    !Local variables
    REAL(4) :: temp,temp2
    temp=1._4/(2._4*dt)
    temp2=1._4/(dt)
    Q(2:ntt-1,:)=(P(3:ntt,:)-P(1:ntt-2,:))*temp
    Q(1  ,:)=(P(2,:)-P(1,:))*temp2
    Q(ntt,:)=(P(ntt,:)-P(ntt-1,:))*temp2
  END SUBROUTINE dev1x_1D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE int1x_1D(Q,P,dt,ntt,nz)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz),INTENT(out):: Q
    REAL(4),DIMENSION(ntt,nz),INTENT(in) :: P
    REAL(4)              ,INTENT(in) :: dt
    INTEGER               ,INTENT(in) :: ntt,nz
    !Local variables
    REAL(4),DIMENSION(ntt,nz) :: S
    INTEGER :: it
    !
    S=0.5_4*P*dt
    !
    Q(1,:)=0._4
    Q(2,:)=S(1,:)+S(2,:)
    DO it=3,ntt
       Q(it,:)=Q(it-1,:)+S(it-1,:)+S(it,:)
    END DO
  END SUBROUTINE int1x_1D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE int1xAdj_1D(Q,P,dt,ntt,nz)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz),INTENT(out):: Q
    REAL(4),DIMENSION(ntt,nz),INTENT(in) :: P
    REAL(4)              ,INTENT(in) :: dt
    INTEGER               ,INTENT(in) :: ntt,nz
    !Local variables
    REAL(4),DIMENSION(ntt,nz) :: S1,S2
    INTEGER :: it
    !
    DO it=1,ntt
       S1(it,:)=P(ntt+1-it,:)
    END DO
    CALL int1x_1D(S2,S1,dt,ntt,nz)
    DO it=1,ntt
       Q(it,:)=S2(ntt+1-it,:)
    END DO
  END SUBROUTINE int1xAdj_1D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE intderiv
