MODULE linop
  !
  use propag, only : calcPropag, calcS0, calcKP1
  use modelspace, only : defcfv, VelocityBounds, taperModel
  use dataspace, only : srcminmax, projectDobs, projectD, taperData
  use outdata!, only !: outmatrix, outmatrix3
  !
  private
  public Flin, UnPMeths
  !
CONTAINS
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Flin(P,Pzero,xi,v,src,rOK,Mtap,ptap,zsrc,dx,dh,dz,dt,vmin,vmax,&
       PMeths,xizero,stap,slist,mode1D2D,comm,npml,noff,ntt,&
       nrcv,ns,nz,nx,nh,nsrc)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nrcv,ns),INTENT(out):: P
    LOGICAL                        ,INTENT(out):: Pzero
    REAL(4),DIMENSION(nz,nx,nh)   ,INTENT(in) :: xi
    REAL(4),DIMENSION(nz,nx)      ,INTENT(in) :: v
    REAL(4),DIMENSION(nsrc)       ,INTENT(in) :: src
    REAL(4),DIMENSION(5)          ,INTENT(in) :: ptap
    REAL(4),DIMENSION(5)          ,INTENT(in) :: Mtap
    REAL(4),DIMENSION(ns)         ,INTENT(in) :: stap
    INTEGER ,DIMENSION(ns)         ,INTENT(in) :: slist
    INTEGER ,DIMENSION(9)          ,INTENT(in) :: PMeths
    LOGICAL ,DIMENSION(nrcv,ns)    ,INTENT(in) :: rOK
    !CHARACTER(len=6)               ,INTENT(in) :: folder
    REAL(4),INTENT(in) :: zsrc,dx,dh,dz,dt,vmin,vmax
    INTEGER ,INTENT(in) :: ntt,nz,noff,nh,nx,ns,npml,nsrc,nrcv,comm,mode1D2D
    LOGICAL ,INTENT(in) :: xizero
    !Local varibles
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: P0,Q
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: X0
    REAL(4),DIMENSION(:,:)  ,ALLOCATABLE :: Aobs
    REAL(4),DIMENSION(:,:)  ,ALLOCATABLE :: cfv,V0,V2
    LOGICAL ,DIMENSION(:)    ,ALLOCATABLE :: posr
    !CHARACTER(len=6) :: suffix
    INTEGER :: is,iis,xmin,xmax,offmin,offmax,nxx,xr,nrcv2,rcvmin,rcvmax,&
         ispos
    INTEGER :: MethJxi,MethXi,MethTap,MethTap2,MethMR,MethTapD,MethQ,MethDcnv,&
         MethAcq2
    INTEGER :: rang,orderK2
    LOGICAL :: forward
    !
#ifdef do_mpi
    INTEGER :: code
    CALL MPI_COMM_RANK(comm,rang,code)
#else
    rang= 0 + 0*comm
#endif
    !
    P(:,:,:)=0._4 ! just to avoid "unused argument" warning
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IF (xizero) THEN
       Pzero=.TRUE. ; RETURN
    ELSE
       Pzero=.FALSE.
    END IF
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL UnPMeths(PMeths,MethJxi,MethXi,MethTap,MethTap2,MethMR,MethTapD,&
         MethQ,MethDcnv,MethAcq2)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL VelocityBounds(v,vmin,vmax,comm,nz,nx)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Initialization
    ALLOCATE(cfv(nz,nx))
    CALL defcfv(cfv,v,MethXi,nz,nx)
    forward = .TRUE.
    orderK2 = 2
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Loop over sources
    DO iis=1,ns
       is=slist(iis)
       WRITE(UNIT=*,FMT="(4(A,I2))")&
            "proc=",rang," iis=",iis,"/",ns," is=",is
       CALL srcminmax(xmin,xmax,offmin,offmax,nxx,nrcv2,rcvmin,rcvmax,ispos,&
            nz,nx,noff,is,MethAcq2,comm)
       !-------------------------------
       !--- Allocations ---------------
       !-------------------------------
       ALLOCATE(P0(ntt,nz,nxx),Q(ntt,nz,nxx),&
            Aobs(ntt,nrcv2),X0(nz,nxx,nh),V0(nz,nxx),V2(nz,nxx),posr(nrcv2))
       !
       X0=xi(:,xmin:xmax,:)
       V0=v(:,xmin:xmax)
       V2=cfv(:,xmin:xmax)
       !
       CALL taperModel(X0,ptap,MethTap,MethTap2,dz,dx,dh,nz,nxx,nh)
       !--------------------------------------------
       !----- Compute P0 ---------------------------
       CALL calcS0(Q,src,zsrc,is,dz,dx,noff,ntt,nz,nxx,nsrc)
       CALL calcPropag(P0,Q,V0,vmin,vmax,dx,dz,dt,forward,mode1D2D,&
            npml,nxx,nz,ntt)
       !--------------------------------------------
       !----- Compute P1 ---------------------------
       CALL calcKP1(Q,P0,X0,V2,dh,dt,xizero,orderK2,forward,mode1D2D,&
            nxx,nz,nh,ntt)
       CALL calcPropag(P0,Q,V0,vmin,vmax,dx,dz,dt,forward,mode1D2D,&
            npml,nxx,nz,ntt)
       !--------------------------------------------
       !----------- Projection Dobs-----------------
       posr=rOK(rcvmin:rcvmax,iis)
       CALL projectDobs(Aobs,P0,posr,zsrc,xr,MethAcq2,dz,ntt,nz,nxx,nrcv2)
       !--------------------------------------------
       !-------------------- taper -----------------
       CALL taperData(Aobs,posr,Mtap,ispos,iis,dx,dz,dt,stap,MethMR,MethTapD,&
            MethAcq2,nsrc,ntt,nrcv2,ns)
       !Aobs=Aobs!*stap(iis)
       !--------------------------------------------
       !--------------- fill big table -------------
       P(:,rcvmin:rcvmax,iis)=Aobs
       !--------------------------------------------
       DEALLOCATE(P0,Q,Aobs,X0,V0,V2,posr)
    END DO
    !--------------------------------------------
    DEALLOCATE(cfv)
  END SUBROUTINE Flin
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE UnPMeths(PMeths,MethJxi,MethXi,MethTap,MethTap2,MethMR,&
       MethTapD,MethQ,MethDcnv,MethAcq2)
    IMPLICIT NONE
    !Parameters
    INTEGER    ,INTENT(out):: MethJxi,MethXi,MethTap,MethTap2,MethMR,MethTapD,&
         MethQ,MethDcnv,MethAcq2
    INTEGER,DIMENSION(9),INTENT(in) :: PMeths
    !
    MethJxi =PMeths(1)
    MethXi  =PMeths(2)
    MethTap =PMeths(3)
    MethTap2=PMeths(4)
    MethMR  =PMeths(5)
    MethTapD=PMeths(6)
    MethQ   =PMeths(7)
    MethDcnv=PMeths(8)
    MethAcq2=PMeths(9)
  END SUBROUTINE UnPMeths
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
END MODULE linop
