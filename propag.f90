MODULE propag
  USE intderiv, ONLY : dev2x,dev1x,int1x,int1xAdj,int2x
  USE findiff, ONLY : findiff1D,findiff2D
  !
  IMPLICIT NONE
  !
  !include "dp.h"
  !include "param.h"
  !
  PRIVATE
  PUBLIC calcS0,calcPropag,calcKP1
CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE calcKP1(K,P0,xiin,cfv,dh,dt,xizero,order,backforw,mode1D2D,&
       nx,nz,nh,ntt) !!PURE
    !Parameters
    REAL(4),DIMENSION(ntt,nz,nx),INTENT(out) :: K
    REAL(4),DIMENSION(ntt,nz,nx),INTENT(in)  :: P0
    REAL(4),DIMENSION(nz,nx,nh) ,INTENT(in)  :: xiin
    REAL(4),DIMENSION(nz,nx)    ,INTENT(in)  :: cfv
    REAL(4)                     ,INTENT(in)  :: dh,dt
    INTEGER                      ,INTENT(in)  :: nx,nz,nh,ntt,order,mode1D2D
    LOGICAL                      ,INTENT(in)  :: xizero,backforw
    !Local variables
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: P0der,xi
    INTEGER :: ih,nh2,xmin,xmax,si,ix,iz,it
    !
    IF (xizero) THEN
       K(:,:,:)=0._4
       RETURN
    END IF
    !
    ALLOCATE(P0der(ntt,nz,nx),xi(nz,nx,nh))
    !
    DO ih=1,nh
       xi(:,:,ih)=dh*cfv*xiin(:,:,ih)
    END DO
    ! determines values for integer si = +1 si = -1
    !  forward propag xi(z, x-h, h) ...etc
    ! backward propag xi(z, x+h, h) ...etc
    IF (backforw) THEN
       si=-1
    ELSE
       si=1
    END IF
    !
    SELECT CASE (order)
    CASE (-11)
       CALL int1xAdj(P0der,P0,dt,ntt,nz,nx)
    CASE (-2)
       CALL int2x(P0der,P0,dt,ntt,nz,nx)
    CASE (-1)
       CALL int1x(P0der,P0,dt,ntt,nz,nx)
    CASE (0)
       P0der=P0
    CASE (1)
       CALL dev1x(P0der,P0,dt,ntt,nz,nx)
    CASE (2)
       CALL dev2x(P0der,P0,dt,ntt,nz,nx)
    CASE DEFAULT
       P0der(:,:,:)=0._4 ! pseudo help -> would look strange
    END SELECT
    !
    K(:,:,:)=0._4
    SELECT CASE (mode1D2D)
    CASE (1) ! no need to do fancy loops in 1D
       DO it = 1,ntt
          K(it,:,:) = P0der(it,:,:)*xi(:,:,1)
       END DO
    CASE (2)
       nh2=(nh-1)/2
       !
       DO ih=-nh2,nh2
          xmin=max(1,1-si*2*ih)
          xmax=min(nx,nx-si*2*ih)
          DO ix=xmin,xmax
             DO iz=1,nz
                K(:,iz,ix)=K(:,iz,ix)+&
                     P0der(:,iz,ix+si*2*ih)*xi(iz,ix+si*ih,ih+nh2+1)
             END DO
          END DO
       END DO
    END SELECT
    !
    DEALLOCATE(P0der,xi)
    !-------------------------------
  END SUBROUTINE calcKP1
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !SUBROUTINE calcKP1(K,P0,xi,cfv,dt,nz,ntt,order)
  !IMPLICIT NONE
  !!Parameters
  !REAL(4),DIMENSION(ntt,nz),INTENT(out) :: K
  !REAL(4),DIMENSION(ntt,nz),INTENT(in)  :: P0
  !REAL(4),DIMENSION(nz)    ,INTENT(in)  :: xi,cfv
  !REAL(4)                  ,INTENT(in)  :: dt
  !INTEGER                  ,INTENT(in)  :: nz,ntt,order
  !!Local variables
  !REAL(4),DIMENSION(ntt,nz) :: P0der
  !REAL(4),DIMENSION(ntt,nz) :: xih
  !!
  !SELECT CASE (order)
  !CASE (-11)
  !   CALL int1xAdj(P0der,P0,dt,ntt,nz)
  !   !P0der=P0
  !   !CALL dev1x(P0der,P0,dt,ntt,nz)
  !   !P0der=-P0der
  !   !CALL dev2x(P0der,P0,dt,ntt,nz)
  !CASE (-1)
  !   CALL int1x_1D(P0der,P0,dt,ntt,nz)
  !   !P0der=P0
  !   !CALL dev1x(P0der,P0,dt,ntt,nz)
  !   !CALL dev2x(P0der,P0,dt,ntt,nz)
  !CASE (2)
  !   CALL dev2x(P0der,P0,dt,ntt,nz)
  !CASE DEFAULT
  !   P0der(:,:)=0._4
  !   WRITE(UNIT=*,FMT="(A)") "error in calcKP1"
  !END SELECT
  !!
  !xih=SPREAD(source=xi,dim=1,ncopies=ntt)
  !K=xih*P0der*SPREAD(source=cfv,dim=1,ncopies=ntt)
  !END SUBROUTINE calcKP1
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE calcPropag(P,S,v0,vmin,vmax,dx,dz,dt,backforw,mode1D2D,&
       npml,nx,nz,ntt)
    !Parameters
    REAL(4),DIMENSION(ntt,nz,nx),INTENT(out) :: P
    REAL(4),DIMENSION(ntt,nz,nx),INTENT(in)  :: S
    REAL(4),DIMENSION(nz,nx)    ,INTENT(in)  :: v0
    REAL(4)                     ,INTENT(in)  :: dx,dz,dt,vmin,vmax
    INTEGER                     ,INTENT(in)  :: npml,nx,nz,ntt,mode1D2D
    LOGICAL                     ,INTENT(in)  :: backforw
    !Local variables
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: Q,K
    !
    IF (backforw) THEN ! forward propagation
       SELECT CASE (mode1D2D)
       CASE (1) ! 1D propagation
          CALL findiff1D(P,S,v0,dz,dt,nz,ntt)
       CASE (2) ! 2D propagation
          CALL findiff2D(P,S,v0,vmin,vmax,dx,dz,dt,npml,nx,nz,ntt)
       END SELECT
    ELSE ! backward propagation
       !
       ALLOCATE(Q(ntt,nz,nx),K(ntt,nz,nx))
       Q(:,:,:)=S(ntt:1:-1,:,:)
       !
       SELECT CASE (mode1D2D)
       CASE (1) ! 1D propagation
          CALL findiff1D(K,Q,v0,dz,dt,nz,ntt)
       CASE (2) ! 2D propagation
          CALL findiff2D(K,Q,v0,vmin,vmax,dx,dz,dt,npml,nx,nz,ntt)
          P(:,:,:)=K(ntt:1:-1,:,:)
       END SELECT
       DEALLOCATE(Q,K)
    END IF
  END SUBROUTINE calcPropag
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE calcS0(S0,src,zsrc,is,dz,dx,noff,ntt,nz,nxx,nsrc)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz,nxx),INTENT(out) :: S0
    REAL(4),DIMENSION(nsrc)      ,INTENT(in)  :: src
    REAL(4)                      ,INTENT(in)  :: dz,dx,zsrc
    INTEGER                      ,INTENT(in)  :: ntt,nz,nxx,nsrc,is,noff
    !Local variables
    INTEGER :: ispos,maxoff,izsrc
    !
    maxoff = (noff-1)/2
    ispos  = MIN(is, maxoff+1)
    izsrc  = FLOOR(zsrc/dz) + 1
    !
    S0(:,:,:)=0._4
    S0(1:nsrc,izsrc,ispos)=src/(dx*dz)
  END SUBROUTINE calcS0
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE propag
