MODULE modelspace
  USE maths!, ONLY: SmoothSin, deftapZ, deftapX, deftapX2, deftapXdisym
  USE intderiv, ONLY: dev2x,dev1x,int1x,int2x
  !
#ifdef do_mpi
  include 'mpif.h'
#endif
  !
  PRIVATE
  PUBLIC taperModel, inttolog3, defcfv, VelocityBounds
  !
CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE taperModel(X0,ptap,MethTap,MethTap2,dz,dx,dh,nz,nxx,nh)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz,nxx,nh),INTENT(inout) :: X0
    REAL(4),DIMENSION(5)        ,INTENT(in)    :: ptap
    REAL(4),INTENT(in) :: dz,dx,dh
    INTEGER ,INTENT(in) :: nz,nxx,nh,MethTap,MethTap2
    !Local variables
    REAL(4),DIMENSION(nz) :: tapZ
    REAL(4),DIMENSION(nxx):: tapX
    REAL(4),DIMENSION(nh) :: tapH
    INTEGER :: izt1,izt2,izb,ixt,iht,iz,ix,ih
    LOGICAL :: DoTapZ,DoTapX,DoTapH
    !
    izt1=FLOOR(ptap(1)/dz)
    izt2=FLOOR(ptap(2)/dz)
    izb=FLOOR(ptap(3)/dz)
    ixt=FLOOR(ptap(4)/dx)
    iht=FLOOR(ptap(5)/dh)
    !
    SELECT CASE (MethTap)
    CASE (0) ! no taper
       X0=X0
    CASE (1) ! edge taper
       tapZ(:)=1._4
       tapX(:)=1._4
       tapH(:)=1._4
       CALL IntToLog3(DoTapZ,DoTapX,DoTapH,MethTap2)
       !-------------------------------
       ! ztaper
       IF (DoTapZ) CALL deftapZ(tapZ,izt1,izt2,izb,nz)
       !-------------------------------
       ! xtaper
       IF (DoTapX) CALL deftapX(tapX,ixt,nxx)
       !-------------------------------
       ! htaper
       IF (DoTapH) CALL deftapX(tapH,iht,nh)
       !-------------------------------
       ! fills the (z,x,h) array
       DO ih=1,nh
          DO ix=1,nxx
             DO iz=1,nz
                X0(iz,ix,ih)=X0(iz,ix,ih)*tapZ(iz)*tapX(ix)*tapH(ih)
             END DO
          END DO
       END DO
       !-------------------------------
    END SELECT
  END SUBROUTINE taperModel
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE IntToLog3(a,b,c,S)
    IMPLICIT NONE
    !Parameters
    LOGICAL,INTENT(out) :: a,b,c
    INTEGER,INTENT(in)  :: S
    !Local variables
    SELECT CASE (S)
    CASE (0)
       a=.FALSE. ; b=.FALSE. ; c=.FALSE.
    CASE (1)
       a=.FALSE. ; b=.FALSE. ; c=.TRUE.
    CASE (2)
       a=.FALSE. ; b=.TRUE.  ; c=.FALSE.
    CASE (3)
       a=.FALSE. ; b=.TRUE.  ; c=.TRUE.
    CASE (4)
       a=.TRUE.  ; b=.FALSE. ; c=.FALSE.
    CASE (5)
       a=.TRUE.  ; b=.FALSE. ; c=.TRUE.
    CASE (6)
       a=.TRUE.  ; b=.TRUE.  ; c=.FALSE.
    CASE DEFAULT
       a=.TRUE.  ; b=.TRUE.  ; c=.TRUE.
    END SELECT
    !
  END SUBROUTINE IntToLog3
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE defcfv(cfv,v,MethXi,nz,nx)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz,nx),INTENT(out) :: cfv
    REAL(4),DIMENSION(nz,nx),INTENT(in)  :: v
    INTEGER                  ,INTENT(in)  :: nz,nx,MethXi
    !Local variables
    SELECT CASE (MethXi)
    CASE (0)
       cfv=4._4/(v**2)
    CASE (1)
       cfv=1._4
    END SELECT
  END SUBROUTINE defcfv
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE VelocityBounds(v,vmin,vmax,comm,nz,nx)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz,nx),INTENT(in) :: v
    REAL(4)                 ,INTENT(in) :: vmin,vmax
    INTEGER                  ,INTENT(in) :: nz,nx
    !Local variables
    LOGICAL,DIMENSION(nz,nx) :: okmin,okmax
    INTEGER                  :: rang,cpt,iz,ix,comm
#ifdef do_mpi
    INTEGER :: code
#endif
    !-----------------------------------
#ifdef do_mpi
    CALL MPI_COMM_RANK(comm,rang,code)
#else
    rang = 0 + 0*comm
#endif
    !
    okmin=(v.GE.vmin)
    okmax=(v.LE.vmax)
    !
    IF (ANY(.NOT.okmin) .OR. ANY(.NOT.okmax)) THEN
       IF (rang.EQ.0) WRITE(UNIT=*,FMT="(A)") "error: velocity out of bounds !!"
       cpt=0
       DO ix=1,nx
          DO iz=1,nz
             IF (.NOT. okmin(iz,ix) .AND. cpt.LE.10) THEN
                cpt=cpt+1
                WRITE(UNIT=*,FMT="(A,I4,A,I4,A,F6.1,A,F6.1,A)")&
                     "v(",iz,",",ix,")=",v(iz,ix)," m/s -- vmin=",vmin
             END IF
             IF (.NOT. okmax(iz,ix) .AND. cpt.LE.10) THEN
                cpt=cpt+1
                WRITE(UNIT=*,FMT="(A,I4,A,I4,A,F6.1,A,F6.1,A)")&
                     "v(",iz,",",ix,")=",v(iz,ix)," m/s -- vmax=",vmax
             END IF
          END DO
       END DO
#ifdef do_mpi
       CALL MPI_BARRIER(comm,code)
       CALL MPI_ABORT(comm,123,code)
#else
       !STOP
#endif
    END IF
  END SUBROUTINE VelocityBounds
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE modelspace
