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
  PUBLIC calcPropag, calcKP1, calcQ, calcS0, calcS0W, calcKL3, calcKL3W,&
       DerivZ, DerivX
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
  !PURE SUBROUTINE calcQ_1D(Q,P,L,orderQ,dt,nh,ntt,nz,nx)
  !IMPLICIT NONE
  !!Parameters
  !REAL(4),DIMENSION(nz,nx,nh) ,INTENT(out) :: Q
  !REAL(4),DIMENSION(ntt,nz,nx),INTENT(in)  :: P,L
  !REAL(4)                     ,INTENT(in)  :: dt
  !INTEGER                      ,INTENT(in)  :: nx,nz,nh,ntt,orderQ
  !!Local variables
  !REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: Pder
  !!INTEGER :: ih,nh2,xmin,xmax,lmin,pmin,ix,nxx,iz
  !!
  !!nh2=(nh-1)/2
  !!
  !ALLOCATE(Pder(ntt,nz,nx))
  !!
  !SELECT CASE (orderQ)
  !CASE (-1)
  !   CALL int1x(Pder,P,dt,ntt,nz,nx)
  !CASE (2)
  !   CALL dev2x(Pder,P,dt,ntt,nz,nx)
  !CASE DEFAULT
  !   Pder(:,:,:) = 0._4
  !END SELECT
  !!
  !Q(:,:,:)=0._4
  !!
  !Q(:,:,1)=dt*SUM(Pder*L,dim=1)
  !!
  !DEALLOCATE(Pder)
  !!
  !END SUBROUTINE calcQ_1D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE calcQ(Q,P,L,orderQ,dt,nh,ntt,nz,nx)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz,nx,nh) ,INTENT(out) :: Q
    REAL(4),DIMENSION(ntt,nz,nx),INTENT(in)  :: P,L
    REAL(4)                     ,INTENT(in)  :: dt
    INTEGER                      ,INTENT(in)  :: nx,nz,nh,ntt,orderQ
    !Local variables
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: Pder
    INTEGER :: ih,nh2,xmin,xmax,lmin,pmin,ix,nxx,iz
    !
    nh2=(nh-1)/2
    !
    ALLOCATE(Pder(ntt,nz,nx))
    !
    SELECT CASE (orderQ)
    CASE (-1)
       CALL int1x(Pder,P,dt,ntt,nz,nx)
    CASE (2)
       CALL dev2x(Pder,P,dt,ntt,nz,nx)
    CASE DEFAULT
       Pder(:,:,:) = 0._4
    END SELECT
    !
    Q(:,:,:)=0._4
    !
    DO ih=-nh2,nh2
       xmin=1 +abs(ih)
       xmax=nx-abs(ih)
       nxx=xmax-xmin+1
       pmin=max(1 ,1 -2*ih)
       lmin=max(1 ,1 +2*ih)
       DO ix=0,nxx-1
          DO iz=1,nz
             Q(iz,xmin+ix,ih+nh2+1)=dt*&
                  SUM(Pder(:,iz,pmin+ix)*L(:,iz,lmin+ix),dim=1)
          END DO
       END DO
    END DO
    !
    DEALLOCATE(Pder)
    !
  END SUBROUTINE calcQ
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
       !
       Q(:,:,:)=S(ntt:1:-1,:,:)
       !
       SELECT CASE (mode1D2D)
       CASE (1) ! 1D propagation
          CALL findiff1D(K,Q,v0,dz,dt,nz,ntt)
       CASE (2) ! 2D propagation
          CALL findiff2D(K,Q,v0,vmin,vmax,dx,dz,dt,npml,nx,nz,ntt)
       END SELECT
       !
       P(:,:,:)=K(ntt:1:-1,:,:)
       !
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
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE calcS0W(S0W,srcdcnv,zsrc,is,dz,dx,noff,ntt,nz,nxx,nsrc)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz,nxx),INTENT(out) :: S0W
    REAL(4),DIMENSION(nsrc)      ,INTENT(in)  :: srcdcnv
    REAL(4),INTENT(in) :: dz,dx,zsrc
    INTEGER ,INTENT(in) :: ntt,nz,nxx,nsrc,is,noff
    !Local variables
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: S0
    !
    ALLOCATE(S0(ntt,nz,nxx))
    CALL calcS0(S0,srcdcnv,zsrc,is,dz,dx,noff,ntt,nz,nxx,nsrc)
    CALL DerivZ(S0W,S0,zsrc,dz,ntt,nz,nxx)
    DEALLOCATE(S0)
  END SUBROUTINE calcS0W
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE calcKL3(K,S,dz,ntt,nz,nx)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz,nx),INTENT(out) :: K
    REAL(4),DIMENSION(ntt,nz,nx),INTENT(in)  :: S
    REAL(4)                     ,INTENT(in)  :: dz
    INTEGER                      ,INTENT(in)  :: nx,nz,ntt
    !Local variables
    K=S/dz
  END SUBROUTINE calcKL3
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE calcKL3W(K,S,zsrc,xr,dz,dx,MethAcq2,ntt,nz,nx)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz,nx),INTENT(out) :: K
    REAL(4),DIMENSION(ntt,nz,nx),INTENT(in)  :: S
    REAL(4)                     ,INTENT(in)  :: dz,dx,zsrc
    INTEGER                      ,INTENT(in)  :: nx,nz,ntt,xr,MethAcq2
    !Local variables
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: Q
    !
    ALLOCATE(Q(ntt,nz,nx))
    Q=S/dz
    IF (MethAcq2.EQ.0) THEN ! surface acquisition
       CALL DerivZ(K,Q,zsrc,dz,ntt,nz,nx)
    ELSE ! VSP acquisition
       CALL DerivX(K,Q,xr,dx,ntt,nz,nx)
    END IF
    DEALLOCATE(Q)
  END SUBROUTINE calcKL3W
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE DerivZ(Q,P,zsrc,dz,ntt,nz,noff)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(out) :: Q
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(in)  :: P
    REAL(4)                       ,INTENT(in)  :: dz,zsrc
    INTEGER                        ,INTENT(in)  :: ntt,nz,noff
    !Local variables
    INTEGER :: izsrc
    izsrc  = FLOOR(zsrc/dz) + 1
    !
    Q(:,:,:) = 0._4
    Q(:,izsrc-1,:) = -P(:,izsrc,:)/(2._4*dz)
    Q(:,izsrc+1,:) =  P(:,izsrc,:)/(2._4*dz)
    !Q(:,izsrc,:)=P(:,izsrc,:)
  END SUBROUTINE DerivZ
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE DerivX(Q,P,xr,dx,ntt,nz,noff)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(out) :: Q
    REAL(4),DIMENSION(ntt,nz,noff),INTENT(in)  :: P
    REAL(4)                       ,INTENT(in)  :: dx
    INTEGER                        ,INTENT(in)  :: ntt,nz,noff,xr
    !Local variables
    !
    Q(:,:,:)=0._4
    Q(:,:,xr-1)=-P(:,:,xr)/(2._4*dx)
    Q(:,:,xr+1)=P(:,:,xr)/(2._4*dx)
    Q=-Q
  END SUBROUTINE DerivX
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE propag
