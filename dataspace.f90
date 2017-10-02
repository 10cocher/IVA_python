MODULE dataspace
  !
  USE maths!, ONLY: SmoothSin, deftapX, deftapX2, deftapXdisym
  !
  PRIVATE
  PUBLIC projectDobs, projectD, taperData, srcminmax,&
       defMR2, taperLine, talpha
  !
CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE projectDobs(Pobs,P,rOK,zsrc,ixrcv,MethAcq2,dz,&
       ntt,nz,nxx,nrcv2)
    ! Operator M of the equations,adjoint of projectD
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nrcv2) ,INTENT(out) :: Pobs
    REAL(4),DIMENSION(ntt,nz,nxx),INTENT(in)  :: P
    LOGICAL ,DIMENSION(nrcv2)     ,INTENT(in)  :: rOK
    REAL(4),INTENT(in) :: dz,zsrc
    INTEGER ,INTENT(in) :: ntt,nz,nxx,nrcv2,MethAcq2,ixrcv
    !Local variables
    INTEGER :: ircv,izsrc
    !
    Pobs(:,:)=0._4
    !
    SELECT CASE (MethAcq2)
    CASE (0) ! surface acquisition
       izsrc = FLOOR(zsrc/dz) + 1
       DO ircv=1,nrcv2 ! nrcv2==nxx
          IF (rOK(ircv)) Pobs(:,ircv)=P(:,izsrc,ircv)
       END DO
    CASE (1) ! VSP acqusitiion
       !ixrcv = FLOOR(xrcv/dx) + 1
       DO ircv=1,nrcv2 ! nrcv2==nz
          IF (rOK(ircv)) Pobs(:,ircv)=P(:,ircv,ixrcv)
       END DO
    CASE DEFAULT
       Pobs(:,:)=0._4
       !IF (rang.EQ.0) WRITE(UNIT=*,FMT="(A)")&
       !     "MethAcq2 not defined in projectDobs"
    END SELECT
    !
  END SUBROUTINE projectDobs
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE projectD(P,Pobs,rOK,zsrc,ixrcv,MethAcq2,dz,nz,nxx,ntt,nrcv2)
    ! Operator M^T of the equations, adjoint of projectDobs
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz,nxx),INTENT(out) :: P
    REAL(4),DIMENSION(ntt,nrcv2) ,INTENT(in)  :: Pobs
    LOGICAL ,DIMENSION(nrcv2)     ,INTENT(in)  :: rOK
    REAL(4),INTENT(in) :: zsrc,dz
    INTEGER ,INTENT(in) :: ntt,nz,nxx,nrcv2,MethAcq2,ixrcv
    !Local variables
    INTEGER :: ircv,izsrc
    !#ifdef do_mpi
    !INTEGER :: code
    !#endif
    !
    !#ifdef do_mpi
    !CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
    !#else
    !rang=0
    !#endif
    !
    P(:,:,:)=0._4
    !
    SELECT CASE (MethAcq2)
    CASE (0) ! surface acquisition
       izsrc = FLOOR(zsrc/dz) + 1
       DO ircv=1,nrcv2 ! nrcv2=nxx
          IF (rOK(ircv)) P(:,izsrc,ircv)=Pobs(:,ircv)
       END DO
    CASE (1) ! VSP acqusitiion
       !ixrcv = FLOOR(xrcv/dx) + 1
       DO ircv=1,nrcv2
          IF (rOK(ircv)) P(:,ircv,ixrcv)=Pobs(:,ircv)
       END DO
    CASE DEFAULT
       P(:,:,:)=0._4
       !IF (rang.EQ.0) WRITE(UNIT=*,FMT="(A)")&
       !     "MethAcq2 not defined in projectD"
    END SELECT
    !
  END SUBROUTINE projectD
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE taperData(P,rOK,Mtap,ispos,iis,dx,dz,dt,stap,&
       MethMR,MethTapD,MethAcq2,nsrc,ntt,nrcv2,ns)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nrcv2),INTENT(inout) :: P
    REAL(4),DIMENSION(5)        ,INTENT(in)    :: Mtap
    REAL(4),DIMENSION(ns)       ,INTENT(in)    :: stap
    LOGICAL ,DIMENSION(nrcv2)    ,INTENT(in)    :: rOK
    REAL(4),INTENT(in) :: dx,dz,dt
    INTEGER ,INTENT(in) :: ntt,nsrc,nrcv2,MethMR,MethTapD,MethAcq2,iis,ispos,ns
    !Local variables
    REAL(4),DIMENSION(:,:),ALLOCATABLE :: MR
    !
    ALLOCATE(MR(ntt,nrcv2))
    !
    CALL defMR2(MR,rOK,Mtap,ispos,dx,dz,dt,nsrc,ntt,nrcv2,&
         MethMR,MethTapD,MethAcq2)
    !
    P=P*MR*stap(iis)
    !
    DEALLOCATE(MR)
    !
  END SUBROUTINE taperData
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE defMR2(MR,rOK,Mtap,ispos,dx,dz,dt,nsrc,ntt,nrcv2,MethMR,MethTapD,&
       MethAcq2)
    ! normally nrcv2=nz or nrcv2=nxx
    IMPLICIT NONE    
    !Parameters
    REAL(4),DIMENSION(ntt,nrcv2),INTENT(out) :: MR
    REAL(4),DIMENSION(5)        ,INTENT(in)  :: Mtap
    LOGICAL ,DIMENSION(nrcv2)    ,INTENT(in)  :: rOK
    REAL(4),INTENT(in) :: dx,dz,dt
    INTEGER ,INTENT(in) :: ntt,nsrc,nrcv2,MethMR,MethTapD,MethAcq2,ispos
    !Local variables
    REAL(4),DIMENSION(:),ALLOCATABLE :: tax,tapT,taprcv
    REAL(4):: expMR,t0,at,Tapod2time,tmax
    INTEGER :: nsrc2,it,Tapod,Xapod,ircv,Tapod2!,rang
    LOGICAL :: TaperTime,TaperRcv,DisymTime
    !
    !#ifdef do_mpi
    !INTEGER :: code
    !#endif
    !!
    !#ifdef do_mpi
    !CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
    !#else
    !rang=0
    !#endif
    !
    nsrc2=(nsrc-1)/2
    DisymTime=.FALSE.
    Tapod2time=0.7_4
    tmax=REAL(ntt-nsrc2,4)*dt
    IF (Tapod2time.GT.tmax .AND. DisymTime) THEN
       !IF (rang.EQ.0) THEN
       !WRITE(UNIT=*,FMT="(A)") "Problem for Tapod2 (Tapod2>tmax) in defMR2."
       !END IF
    END IF
    Tapod2=FLOOR(Tapod2time/dt)+nsrc2
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SELECT CASE (MethTapD)
    CASE (0)
       TaperTime=.FALSE.
       TaperRcv =.FALSE.
    CASE (1)
       TaperTime=.TRUE.
       TaperRcv =.FALSE.
    CASE (2)
       TaperTime=.FALSE.
       TaperRcv =.TRUE.
    CASE (3)
       TaperTime=.TRUE.
       TaperRcv =.TRUE.
       !CASE DEFAULT
       !IF (rang.EQ.0) WRITE(UNIT=*,FMT="(A)") "MethTapD not defined in defMR2"
    END SELECT
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    expMR=Mtap(1)
    Tapod=FLOOR(Mtap(2)/dt)
    SELECT CASE (MethAcq2)
    CASE (0) ! surface acquisition
       Xapod=FLOOR(Mtap(3)/dx)
    CASE (1) ! VSP acquisition
       Xapod=FLOOR(Mtap(3)/dz)
    CASE DEFAULT
       MR(:,:)=0._4
       RETURN
       !IF (rang.EQ.0) WRITE(UNIT=*,FMT="(A)") "MethAcq2 not defined in defMR2"
    END SELECT
    t0=Mtap(4)
    at=Mtap(5)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MR(:,:)=0._4
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! First step : muliply data by the time if requested
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ALLOCATE(tax(ntt))
    !
    SELECT CASE (MethMR)
    CASE (0) ! M=1
       MR(:,:)=1._4
    CASE (1) ! M=t^alpha
       tax(1:nsrc2)=0._4
       tax(nsrc2+1:ntt)=(/ ( (REAL(it-1,4)*dt)**(expMR),it=1,ntt-nsrc2 )  /)
       DO ircv=1,nrcv2
          MR(:,ircv)=tax
       END DO
    CASE (2) ! M=(t-t(offset))^alpha
       ! we choose t(offset) = t0 + at*abs(offset)
       DO ircv=1,nrcv2
          CALL talpha(tax,expMR,t0,at,dx,dt,ircv,ispos,ntt,nsrc2)
          MR(:,ircv)=tax
       END DO
    CASE DEFAULT
       MR(:,:)=0._4
       RETURN
       !IF (rang.EQ.0) WRITE(UNIT=*,FMT="(A)") "MethMR not defined in defMR2"
    END SELECT
    !
    DEALLOCATE(tax)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Second step : time taper (beginning and end of recording)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ALLOCATE(tapT(ntt))
    !
    IF (TaperTime) THEN
       IF (DisymTime) THEN
          CALL deftapXdisym(tapT,Tapod,Tapod2,ntt)
       ELSE
          CALL deftapX(tapT,Tapod,ntt)
       END IF
    ELSE
       tapT(:)=1._4
    END IF
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Third step : receiver taper
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ALLOCATE(taprcv(nrcv2))
    !
    IF (TaperRcv) THEN
       !-------------------------------
       ! data tapered
       CALL taperLine(taprcv,rOK,nrcv2,Xapod)
       !-------------------------------
    ELSE
       !-------------------------------
       ! data not tapered
       DO ircv=1,nrcv2
          IF (rOK(ircv)) THEN
             taprcv(ircv)=1._4
          ELSE
             taprcv(ircv)=0._4
          END IF
       END DO
       !-------------------------------
    END IF
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Last step : applies tapers to MR
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DO ircv=1,nrcv2
       DO it=1,ntt
          MR(it,ircv)=MR(it,ircv)*taprcv(ircv)*tapT(it)
       END DO
    END DO
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DEALLOCATE(tapT,taprcv)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  END SUBROUTINE defMR2
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE taperLine(tap,ok,nxx,Xapod)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(nxx),INTENT(out) :: tap
    LOGICAL ,DIMENSION(nxx),INTENT(in)  :: ok
    INTEGER                ,INTENT(in)  :: Xapod,nxx
    !Local variables
    INTEGER :: rmin,rmax,nr,ix,Xwindow
    LOGICAL :: up,down
    !
    tap(:)=0._4
    DO ix=1,nxx
       !-------------------------------
       ! defines up and down
       IF (ix.EQ.1) THEN
          up=ok(ix)
       ELSE
          up=.NOT.ok(ix-1) .AND. ok(ix)
       END IF
       IF (ix.EQ.nxx) THEN
          down=ok(ix)
       ELSE
          down=ok(ix) .AND. .NOT. ok(ix+1)
       END IF
       !
       IF (up) rmin=ix
       IF (down) THEN
          rmax=ix
          nr=rmax-rmin+1
          IF (nr.LT.5) THEN
             tap(rmin:rmax)=1._4
          END IF
          Xwindow=min(Xapod,max(2,nr/3))
          CALL deftapX2(tap(rmin:rmax),Xwindow,nr)
       END IF
    END DO
  END SUBROUTINE taperLine
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE talpha(tapT,expMR,t0,at,dx,dt,ircv,ispos,ntt,nsrc2)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt),INTENT(out) :: tapT
    REAL(4),INTENT(in) :: expMR,t0,at,dx,dt
    INTEGER ,INTENT(in) :: ntt,nsrc2,ispos,ircv
    !Local variables
    REAL(4):: xrcv,tau,seuil,den
    INTEGER :: it,itau,ifin
    !
    xrcv=dx*ABS(REAL(ispos-ircv,4))
    !tau=t0+at*xrcv
    tau=sqrt( t0**2 + (at*xrcv)**2  )
    itau=nsrc2+1+FLOOR(tau/dt)
    !
    IF (itau.GE.ntt) THEN
       tapT(:)=0._4
       RETURN
    END IF

    !
    tapT(1:itau)=0._4
    tapT(itau+1:ntt)=(/ ( (REAL(it-itau,4)*dt)**(expMR),it=itau+1,ntt)  /)
    !
    seuil=10._4
    !
    IF (expMR.LT.0._4) THEN
       ifin=ntt
       DO it=itau+1,ntt
          IF (tapT(it).LT.seuil) THEN
             ifin=it
             EXIT
          END IF
       END DO
       den=1._4/REAL(ifin-itau,4)
       DO it=itau,ifin
          IF (tapT(it).GT.seuil)&
               tapT(it)=seuil*SmoothSin(den*REAL(it-itau,4))
       END DO
    END IF
    !
  END SUBROUTINE talpha
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE srcminmax(xmin,xmax,offmin,offmax,nxx,nrcv2,rcvmin,rcvmax,&
       ispos,nz,nx,noff,is,MethAcq2,comm)
    !Parameters
    INTEGER,INTENT(out) :: xmin,xmax,offmin,offmax,nxx,nrcv2,rcvmin,rcvmax,&
         ispos
    INTEGER,INTENT(in)  :: nz,nx,noff,is,MethAcq2,comm
    !Local variables
    INTEGER :: maxoff,miloff,rang
    !
#ifdef do_mpi
    INTEGER :: code
    CALL MPI_COMM_RANK(comm,rang,code)
#else
    rang = 0 + 0*comm
#endif
    !
    maxoff=(noff-1)/2
    miloff=maxoff+1
    !
    !-------------------------------
    ! source position
    ispos=MIN(is,maxoff+1)
    !-------------------------------
    xmin=max(1,is-maxoff)
    xmax=min(nx,is+maxoff)
    offmin=max(1   ,miloff+1 -is)
    offmax=min(noff,miloff+nx-is)
    nxx=xmax-xmin+1
    !-------------------------------
    SELECT CASE (MethAcq2)
    CASE (0) ! surface acquisition
       xr=0
       nrcv2=nxx
       rcvmin=offmin
       rcvmax=offmax
    CASE (1) ! VSP acquisition
       IF(xmax.LT.nx) THEN
          WRITE(UNIT=*,FMT="(A,I4)")&
               "pb VSP: xr not in the range of offsets for is=",is
       END IF
       !xr=nx-is-1+maxoff
       xr=nxx-1
       !
       nrcv2=nz
       rcvmin=1
       rcvmax=nz
    CASE DEFAULT
       xr=0
       nrcv2=0
       rcvmin=0
       rcvmax=0
       IF (rang.EQ.0) THEN
          WRITE(UNIT=*,FMT="(A)") "MethAcq2 not defined in srcminmax"
       END IF
    END SELECT
    !-------------------------------
  END SUBROUTINE srcminmax
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE dataspace
