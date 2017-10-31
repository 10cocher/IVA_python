MODULE linop
  !
  use propag, only : calcPropag, calcKP1, calcQ, calcS0, calcS0W, calcKL3, calcKL3W
  use modelspace, only : defcfv, defcfvW, VelocityBounds, taperModel, Wmodel, WmodelT
  use dataspace, only : srcminmax, projectDobs, projectD, taperData
  !use outdata!, only !: outmatrix, outmatrix3
  !
  implicit none
  !
#ifdef do_mpi
  include 'mpif.h'
#endif
  !
  private
  public Flin, Fadj, Fdag,&
       ReduceZXH, UnPMeths
  !
contains
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Flin(P,Pzero,xi,xizero,v,src,rOK,Mtap,ptap,zsrc,dx,dh,dz,dt,vmin,vmax,&
       PMeths,stap,slist,mode1D2D,comm,npml,noff,ntt,&
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
    rang = 0 + 0*comm
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
       !WRITE(UNIT=*,FMT="(4(A,I2))")&
       !"proc=",rang," iis=",iis,"/",ns," is=",is
       CALL srcminmax(xmin,xmax,offmin,offmax,nxx,nrcv2,rcvmin,rcvmax,ispos,&
            nz,nx,noff,is,MethAcq2,comm)
       !-------------------------------
       !--- Allocations ---------------
       !-------------------------------
       ALLOCATE(P0(ntt,nz,nxx),Q(ntt,nz,nxx),&
            Aobs(ntt,nrcv2),X0(nz,nxx,nh),V0(nz,nxx),V2(nz,nxx),posr(nrcv2))
       !
       X0   = xi(:,xmin:xmax,:)
       V0   = v(:,xmin:xmax)
       V2   = cfv(:,xmin:xmax)
       posr = rOK(rcvmin:rcvmax,iis)
       !--------------------------------------
       !--- taper reflectivity ---------------
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
       CALL projectDobs(Aobs,P0,posr,zsrc,xr,MethAcq2,dz,ntt,nz,nxx,nrcv2)
       !--------------------------------------------
       !---------------- taper data ----------------
       CALL taperData(Aobs,posr,Mtap,ispos,iis,dx,dz,dt,stap,MethMR,MethTapD,&
            MethAcq2,nsrc,ntt,nrcv2,ns,.FALSE.)
       !--------------------------------------------
       !------ fill big (output) array -------------
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
  SUBROUTINE Fadj(xi,xizero,P,Pzero,v,src,rOK,Mtap,ptap,zsrc,dx,dh,dz,dt,&
       vmin,vmax,PMeths,stap,slist,mode1D2D,comm,npml,noff,nh,&
       ntt,nz,nx,ns,nsrc,nrcv)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz,nx,nh)   ,INTENT(out):: xi
    LOGICAL                        ,INTENT(out):: xizero
    REAL(4),DIMENSION(ntt,nrcv,ns),INTENT(in) :: P
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
    LOGICAL ,INTENT(in) :: Pzero
    !Local varibles
    REAL(4),DIMENSION(:,:,:)  ,ALLOCATABLE :: P0,L1,Q0,tmpZXHS
    REAL(4),DIMENSION(:,:,:)  ,ALLOCATABLE :: X0
    REAL(4),DIMENSION(:,:)    ,ALLOCATABLE :: Aobs
    REAL(4),DIMENSION(:,:)    ,ALLOCATABLE :: cfv,V0,V2
    LOGICAL ,DIMENSION(:)      ,ALLOCATABLE :: posr
    INTEGER :: is,iis,xmin,xmax,offmin,offmax,nxx,xr,nrcv2,rcvmin,rcvmax,&
         ispos
    INTEGER :: MethJxi,MethXi,MethTap,MethTap2,MethMR,MethTapD,MethQ,MethDcnv,&
         MethAcq2
    INTEGER :: rang,orderQ2,ih
    LOGICAL :: forward,backward
    !
#ifdef do_mpi
    INTEGER :: code
    CALL MPI_COMM_RANK(comm,rang,code)
#else
    rang = 0 + 0*comm
#endif
    xi(:,:,:)=0._4
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IF (Pzero) THEN
       xizero=.TRUE. ; RETURN
    ELSE
       xizero=.FALSE.
    END IF
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL UnPMeths(PMeths,MethJxi,MethXi,MethTap,MethTap2,MethMR,MethTapD,&
         MethQ,MethDcnv,MethAcq2)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL VelocityBounds(v,vmin,vmax,comm,nz,nx)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Initialization
    ALLOCATE(tmpZXHS(nz,nx,nh),cfv(nz,nx))
    CALL defcfv(cfv,v,MethXi,nz,nx)
    tmpZXHS(:,:,:) = 0._4
    forward  = .TRUE.
    backward = .FALSE.
    orderQ2  = 2
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Loop over sources
    DO iis=1,ns
       is=slist(iis)
       CALL srcminmax(xmin,xmax,offmin,offmax,nxx,nrcv2,rcvmin,rcvmax,ispos,&
            nz,nx,noff,is,MethAcq2,comm)
       !-------------------------------------------------------
       !--- allocations ---------------------------------------
       ALLOCATE(P0(ntt,nz,nxx),L1(ntt,nz,nxx),Q0(ntt,nz,nxx),&
            Aobs(ntt,nrcv2),X0(nz,nxx,nh),V0(nz,nxx),V2(nz,nxx),posr(nrcv2))
       !
       V0   = v(:,xmin:xmax)
       V2   = cfv(:,xmin:xmax)
       Aobs = P(:,rcvmin:rcvmax,iis)
       posr = rOK(rcvmin:rcvmax,iis)
       !-------------------------------------------------------
       !--- taper data ----------------------------------------
       CALL taperData(Aobs,posr,Mtap,ispos,iis,dx,dz,dt,stap,MethMR,MethTapD,&
            MethAcq2,nsrc,ntt,nrcv2,ns,.FALSE.)
       !-------------------------------------------------------
       !--- project Dobs space -> D space ---------------------
       CALL projectD(L1,Aobs,posr,zsrc,xr,MethAcq2,dz,nz,nxx,ntt,nrcv2)
       !-------------------------------------------------------
       !--- backpropagate data (L1 wavefield) -----------------
       CALL calcKL3(Q0,L1,dz,ntt,nz,nxx)
       CALL calcPropag(L1,Q0,V0,vmin,vmax,dx,dz,dt,backward,mode1D2D,&
            npml,nxx,nz,ntt)
       !-------------------------------------------------------
       !--- forward propagate source wavelet (P0 wavefield) ---
       CALL calcS0(Q0,src,zsrc,is,dz,dx,noff,ntt,nz,nxx,nsrc)
       CALL calcPropag(P0,Q0,V0,vmin,vmax,dx,dz,dt,forward,mode1D2D,&
            npml,nxx,nz,ntt)
       !-------------------------------------------------------
       !--- Compute correlation between P0 and L1 wavefields --
       CALL calcQ(X0,P0,L1,orderQ2,dt,nh,ntt,nz,nxx)
       !-------------------------------------------------------
       !--- taper reflectivity models and apply coeff fct of v
       CALL taperModel(X0,ptap,MethTap,MethTap2,dz,dx,dh,nz,nxx,nh)
       DO ih=1,nh
          X0(:,:,ih)=X0(:,:,ih)*V2
       END DO
       !-------------------------------------------------------
       !--- adds the results into tmpZXHS ---------------------
       tmpZXHS(:,xmin:xmax,:)=tmpZXHS(:,xmin:xmax,:)+X0!*stap(iis)
       !-------------------------------------------------------
       !--- Deallocations -------------------------------------
       DEALLOCATE(P0,L1,Q0,Aobs,X0,V0,V2,posr)
    END DO
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !--------------Reduce over processors-------------------
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL ReduceZXH(xi,tmpZXHS,dx,comm,nz,nx,nh)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !--------------Deallocations----------------------------
    DEALLOCATE(tmpZXHS,cfv)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  END SUBROUTINE Fadj
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Fdag(xi,xizero,P,Pzero,v,srcdcnv,rOK,Mtap,ptap,zsrc,dx,dh,dz,dt,ds,&
       vmin,vmax,PMeths,stap,slist,mode1D2D,comm,npml,noff,nh,&
       ntt,nz,nx,ns,nsrc,nrcv)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz,nx,nh)   ,INTENT(out):: xi
    LOGICAL                        ,INTENT(out):: xizero
    REAL(4),DIMENSION(ntt,nrcv,ns),INTENT(in) :: P
    REAL(4),DIMENSION(nz,nx)      ,INTENT(in) :: v
    REAL(4),DIMENSION(nsrc)       ,INTENT(in) :: srcdcnv
    REAL(4),DIMENSION(5)          ,INTENT(in) :: ptap
    REAL(4),DIMENSION(5)          ,INTENT(in) :: Mtap
    REAL(4),DIMENSION(ns)         ,INTENT(in) :: stap
    INTEGER ,DIMENSION(ns)         ,INTENT(in) :: slist
    INTEGER ,DIMENSION(9)          ,INTENT(in) :: PMeths
    LOGICAL ,DIMENSION(nrcv,ns)    ,INTENT(in) :: rOK
    !CHARACTER(len=6)               ,INTENT(in) :: folder
    REAL(4),INTENT(in) :: zsrc,dx,dh,dz,dt,ds,vmin,vmax
    INTEGER ,INTENT(in) :: ntt,nz,noff,nh,nx,ns,npml,nsrc,nrcv,comm,mode1D2D
    LOGICAL ,INTENT(in) :: Pzero
    !Local varibles
    REAL(4),DIMENSION(:,:,:)  ,ALLOCATABLE :: P0,L1,Q0,tmpZXHS
    REAL(4),DIMENSION(:,:,:)  ,ALLOCATABLE :: X0,Y0
    REAL(4),DIMENSION(:,:)    ,ALLOCATABLE :: Aobs
    REAL(4),DIMENSION(:,:)    ,ALLOCATABLE :: cfvW,V0,V4
    LOGICAL ,DIMENSION(:)      ,ALLOCATABLE :: posr
    INTEGER :: is,iis,xmin,xmax,offmin,offmax,nxx,xr,nrcv2,rcvmin,rcvmax,&
         ispos
    INTEGER :: MethJxi,MethXi,MethTap,MethTap2,MethMR,MethTapD,MethQ,MethDcnv,&
         MethAcq2
    INTEGER :: rang,orderQ1
    LOGICAL :: forward,backward
    !
#ifdef do_mpi
    INTEGER :: code
    CALL MPI_COMM_RANK(comm,rang,code)
#else
    rang = 0 + 0*comm
#endif
    xi(:,:,:)=0._4
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IF (Pzero) THEN
       xizero=.TRUE. ; RETURN
    ELSE
       xizero=.FALSE.
    END IF
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL UnPMeths(PMeths,MethJxi,MethXi,MethTap,MethTap2,MethMR,MethTapD,&
         MethQ,MethDcnv,MethAcq2)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL VelocityBounds(v,vmin,vmax,comm,nz,nx)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Initialization
    ALLOCATE(tmpZXHS(nz,nx,nh),cfvW(nz,nx))
    CALL defcfvW(cfvW,v,MethXi,mode1D2D,nz,nx)
    tmpZXHS(:,:,:) = 0._4
    forward  = .TRUE.
    backward = .FALSE.
    orderQ1  = -1
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Loop over sources
    DO iis=1,ns
       is=slist(iis)
       CALL srcminmax(xmin,xmax,offmin,offmax,nxx,nrcv2,rcvmin,rcvmax,ispos,&
            nz,nx,noff,is,MethAcq2,comm)
       !-------------------------------------------------------
       !--- allocations ---------------------------------------
       ALLOCATE(P0(ntt,nz,nxx),L1(ntt,nz,nxx),Q0(ntt,nz,nxx),&
            Aobs(ntt,nrcv2),X0(nz,nxx,nh),Y0(nz,nxx,nh),&
            V0(nz,nxx),V4(nz,nxx),posr(nrcv2))
       !
       V0   = v(:,xmin:xmax)
       V4   = cfvW(:,xmin:xmax)
       Aobs = P(:,rcvmin:rcvmax,iis)
       posr = rOK(rcvmin:rcvmax,iis)
       !-------------------------------------------------------
       !--- taper data ----------------------------------------
       CALL taperData(Aobs,posr,Mtap,ispos,iis,dx,dz,dt,stap,MethMR,MethTapD,&
            MethAcq2,nsrc,ntt,nrcv2,ns,.TRUE.)
       !-------------------------------------------------------
       !--- project Dobs space -> D space ---------------------
       CALL projectD(L1,Aobs,posr,zsrc,xr,MethAcq2,dz,nz,nxx,ntt,nrcv2)
       !-------------------------------------------------------
       !--- backpropagate data (L1 wavefield) -----------------
       CALL calcKL3W(Q0,L1,zsrc,xr,dz,dx,MethAcq2,ntt,nz,nxx)
       CALL calcPropag(L1,Q0,V0,vmin,vmax,dx,dz,dt,backward,mode1D2D,&
            npml,nxx,nz,ntt)
       !-------------------------------------------------------
       !--- forward propagate source wavelet (P0 wavefield) ---
       CALL calcS0W(Q0,srcdcnv,zsrc,is,dz,dx,noff,ntt,nz,nxx,nsrc)
       CALL calcPropag(P0,Q0,V0,vmin,vmax,dx,dz,dt,forward,mode1D2D,&
            npml,nxx,nz,ntt)
       !-------------------------------------------------------
       !--- Compute correlation between P0 and L1 wavefields --
       CALL calcQ(Y0,P0,L1,orderQ1,dt,nh,ntt,nz,nxx)
       !-------------------------------------------------------
       !--- taper reflectivity models and apply coeff fct of v
       CALL Wmodel(X0,Y0,V4,ptap,MethTap,MethTap2,dz,dx,dh,nz,nxx,nh)
       !-------------------------------------------------------
       !--- adds the results into tmpZXHS ---------------------
       tmpZXHS(:,xmin:xmax,:)=tmpZXHS(:,xmin:xmax,:)+X0!*stap(iis)
       !-------------------------------------------------------
       !--- Deallocations -------------------------------------
       DEALLOCATE(P0,L1,Q0,Aobs,X0,Y0,V0,V4,posr)
    END DO
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !--------------Reduce over processors-------------------
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL ReduceZXH(xi,tmpZXHS,dx,comm,nz,nx,nh)
    xi = xi*ds
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !--------------Deallocations----------------------------
    DEALLOCATE(tmpZXHS,cfvW)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  END SUBROUTINE Fdag
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
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE ReduceZXH(zxh,tmpZXHS,dx,comm,nz,nx,nh)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nz,nx,nh),INTENT(out):: zxh
    REAL(4),DIMENSION(nz,nx,nh),INTENT(in) :: tmpZXHS
    REAL(4)                    ,INTENT(in) :: dx
    INTEGER                     ,INTENT(in) :: nz,nx,nh,comm
    !Local variables
    INTEGER :: code
    !
#ifdef do_mpi
    CALL MPI_ALLREDUCE(tmpZXHS,zxh,nz*nx*nh,MPI_REAL,MPI_SUM,comm,code)
#else
    zxh=tmpZXHS
    code = comm ! just to avoid "unused variable comm" while compiling...
#endif
    zxh=zxh*dx
  END SUBROUTINE ReduceZXH
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !
END MODULE linop
