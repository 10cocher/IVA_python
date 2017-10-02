!#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
MODULE findiff
  !
  USE intderiv, ONLY : int1x
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC findiff1D,findiff2D,&
       extendmodel,defnsig,defparamcpml2d,defparamcpml
CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE findiff1D(P,S,v,dz,dt,nz,ntt)
    IMPLICIT NONE
    !Parameters
    REAL(4),DIMENSION(ntt,nz),INTENT(out) :: P
    REAL(4),DIMENSION(ntt,nz),INTENT(in)  :: S
    REAL(4),DIMENSION(nz)    ,INTENT(in)  :: v
    REAL(4)                  ,INTENT(in)  :: dz,dt
    INTEGER                  ,INTENT(in)  :: nz,ntt
    !Local variables
    REAL(4),DIMENSION(:,:),ALLOCATABLE :: Sext,Pext
    REAL(4),DIMENSION(:)  ,ALLOCATABLE :: vext
    INTEGER                            :: it,iz,izm,izp
    !
    ALLOCATE(Sext(ntt,5*nz),Pext(ntt,5*nz),vext(5*nz))
    !
    Sext(:,:)=0._4
    Sext(:,2*nz+1:3*nz)=S
    !
    Pext(:,:)=0._4
    !
    vext(2*nz+1:3*nz)=v
    vext(1:2*nz)=v(1)
    vext(3*nz+1:5*nz)=v(nz);
    !
    DO it=2,ntt-1
       DO iz=1,5*nz
          izm=max(1,iz-1)
          izp=min(5*nz,iz+1)
          !
          Pext(it+1,iz)=2._4*Pext(it,iz)-Pext(it-1,iz)&
               +(vext(iz)*dt/dz)**2* (Pext(it,izm)&
               -2._4*Pext(it,iz)+Pext(it,izp))&
               +Sext(it,iz)*(dt*vext(iz))**2
       END DO
    END DO
    P=Pext(:,2*nz+1:3*nz)
    !
    DEALLOCATE(Sext,Pext,vext)
    !
  END SUBROUTINE findiff1D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE findiff2D(P,S,v0,vmin,vmax,dx,dz,dt,npml,nx,nz,nt)
    !Parameters
    REAL(4),DIMENSION(nt,nz,nx),INTENT(out) :: P
    REAL(4),DIMENSION(nt,nz,nx),INTENT(in)  :: S
    REAL(4),DIMENSION(nz,nx)   ,INTENT(in)  :: v0
    REAL(4)                    ,INTENT(in)  :: dx,dz,dt,vmin,vmax
    INTEGER                    ,INTENT(in)  :: nx,nz,nt,npml
    ! Internal variables
    REAL(4),DIMENSION(:,:,:),ALLOCATABLE :: Q
    REAL(4),DIMENSION(:,:),ALLOCATABLE :: uz, ux, vz, vx, pp, m0e
    REAL(4),DIMENSION(:,:),ALLOCATABLE :: uzn, uxn, vzn, vxn
    REAL(4),DIMENSION(:,:),ALLOCATABLE :: dpdz, dvzdz, Qe, dpdx, dvxdx
    REAL(4),DIMENSION(:,:),ALLOCATABLE :: m0
    INTEGER :: it,iz,ix,nze,nxe,nzext,nxext
    !
    REAL(4),DIMENSION(:)  ,ALLOCATABLE :: sigpmlza,sigpmlzb,sigpmlxa,sigpmlxb
    INTEGER :: nzsig,nxsig
    !-----------------------------------
    ! for the PML
    CALL defnsig(nzsig,nxsig,nz,nx,npml)
    ALLOCATE(sigpmlza(nzsig),sigpmlzb(nzsig-1),&
         sigpmlxa(nxsig),sigpmlxb(nxsig-1))
    CALL defparamcpml2d(sigpmlza,sigpmlzb,sigpmlxa,sigpmlxb,vmin,vmax,npml,&
         dz,dx,nz,nx,nzsig,nxsig)
    !-----------------------------------
    ! creates m0 from v0
    ALLOCATE(m0(nz,nx))
    m0 = 1._4 / (v0**2)
    !-----------------------------------
    ! Integration in time of the source along the third dimension
    !CALL defintegrate(Q,S,dt,nx,nz,nt)
    ALLOCATE(Q(nt,nz,nx))
    CALL int1x(Q,S,dt,nt,nz,nx)
    !
    ! Extend the model (for PML)
    !
    nxext = int((nxsig-nx)/2)
    nxe   = nx + 2*nxext
    nzext = int((nzsig-nz)/2)
    nze   = nz + 2*nzext
    !
    ! Allocate memory
    !
    allocate(uz(nze,nxe),vz(nze-1,nxe),uzn(nze,nxe), &
         vzn(nze-1,nxe),pp(nze,nxe),dpdz(nze-1,nxe),dvzdz(nze-1,nxe), &
         Qe(nze,nxe))
    !
    uz  = 0._4
    vz  = 0._4
    uzn = 0._4
    vzn = 0._4
    pp  = 0._4
    !
    allocate(ux(nze,nxe),vx(nze,nxe-1),uxn(nze,nxe),vxn(nze,nxe-1), &
         dpdx(nze,nxe-1),dvxdx(nze,nxe-1))
    ux  = 0._4
    vx  = 0._4
    uxn = 0._4
    vxn = 0._4
    !
    ! Extend the model
    !
    allocate(m0e(nze,nxe))
    !
    CALL extendmodel(m0e,m0,nzext,nxext,nz,nx)
    !
    !
    !------------------------------------------------------------------
    !
    P(:,:,:) = 0._4
    !
    ! 2d propagation
    !
    do it = 1, nt-1
       !WRITE(UNIT=*,FMT="(A3,I3,A1,I3)"),"it=",it,"/",nt-1
       !
       ! Derivative in space
       !
       dpdz(:,:) = 0._4
       !
       do ix = 1, nxe
          do iz = 1, nze-1
             dpdz(iz,ix) = (pp(iz+1,ix)-pp(iz,ix))/dz
          end do
       end do
       !
       ! Derivative in space
       !
       dpdx = 0._4
       !
       do ix = 1, nxe-1
          do iz = 1, nze
             dpdx(iz,ix) = (pp(iz,ix+1)-pp(iz,ix))/dx
          end do
       end do
       !
       ! Integration in time
       !
       do ix = 1, nxe
          do iz = 1, nze-1
             vzn(iz,ix) = vz(iz,ix) + dpdz(iz,ix) * dt - &
                  sigpmlzb(iz) * vz(iz,ix) * dt
          end do
       end do
       !
       do ix = 1, nxe-1
          do iz = 1, nze
             vxn(iz,ix) = vx(iz,ix) + dpdx(iz,ix) * dt - &
                  sigpmlxb(ix) * vx(iz,ix) * dt
          end do
       end do
       !
       ! Derivative in time (for the next time)
       !
       dvzdz = 0._4
       !
       do ix = 1, nxe
          do iz = 2, nze-1
             dvzdz(iz,ix) = (vzn(iz,ix) - vzn(iz-1,ix))/dz
          end do
       end do
       !
       dvxdx = 0._4
       !
       do ix = 2, nxe-1
          do iz = 1, nze
             dvxdx(iz,ix) = (vxn(iz,ix) - vxn(iz,ix-1)) / dx
          end do
       end do
       !
       ! Selection of the source term
       !
       Qe = 0._4
       !
       do ix = 1, nx
          do iz = 1, nz
             Qe(iz+nzext,ix+nxext) = Q(it,iz,ix)
          end do
       end do
       !
       ! Integration in time
       !
       do ix = 2, nxe-1
          do iz = 2, nze-1
             uzn(iz,ix) = uz(iz,ix) + &
                  (dvzdz(iz,ix) + Qe(iz,ix)) * dt/m0e(iz,ix) - &
                  sigpmlza(iz)*uz(iz,ix)*dt
             uxn(iz,ix) = ux(iz,ix) + &
                  (dvxdx(iz,ix) + Qe(iz,ix)) * dt/m0e(iz,ix) - &
                  sigpmlxa(ix)*ux(iz,ix)*dt
          end do
       end do
       !
       ! Total pressure
       !
       pp = uxn + uzn
       !
       do ix = 1, nx
          do iz = 1, nz
             P(it+1,iz,ix) = pp(iz+nzext,ix+nxext)
          end do
       end do
       !
       ! Update
       !
       ux = uxn
       uz = uzn
       vx = vxn
       vz = vzn
       !
    end do
    !
    P=0.5_4*P
    !
    !------------------------------------------------------------------
    !
    ! Deallocate
    !
    deallocate(uz,vz,uzn,vzn,pp,dpdz,dvzdz,Qe)
    deallocate(ux,vx,uxn,vxn,dpdx,dvxdx,m0e)
    !
    DEALLOCATE(Q,m0,sigpmlza,sigpmlzb,sigpmlxa,sigpmlxb)
    !
    !------------------------------------------------------------------
    !
  end subroutine findiff2D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE extendmodel(m0e,m0,nzext,nxext,nz,nx)
    !Parameters
    REAL(4),DIMENSION(nz+2*nzext,nx+2*nxext),INTENT(out) :: m0e
    REAL(4),DIMENSION(nz,nx)                ,INTENT(in)  :: m0
    INTEGER                                 ,INTENT(in) :: nzext,nxext,nz,nx
    !Local variables
    INTEGER :: iz, ix
    !
    m0e=0._4
    ! Inside
    do ix = 1, nx
       do iz = 1, nz
          m0e(iz+nzext,ix+nxext) = m0(iz,ix)
       end do
    end do
    !
    ! Edges
    !
    do ix = 1, nx
       do iz = 1, nzext
          m0e(iz,ix+nxext) = m0(1,ix)
          m0e(iz+nz+nzext,ix+nxext) = m0(nz,ix)
       end do
    end do
    !
    do ix= 1, nxext
       do iz = 1, nz
          m0e(iz+nzext,ix) = m0(iz,1)
          m0e(iz+nzext,ix+nx+nxext) = m0(iz,nx)
       end do
    end do
    !
    ! Corners
    !
    do ix = 1, nxext
       do iz = 1, nzext
          m0e(iz,ix) = m0(1,1)
          m0e(iz,ix+nx+nxext) = m0(1,nx)
          m0e(iz+nz+nzext,ix) = m0(nz,1)
          m0e(iz+nz+nzext,ix+nx+nxext) = m0(nz,nx)
       end do
    end do
    !
  END SUBROUTINE extendmodel
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE defnsig(nzsig,nxsig,nz,nx,npml)
    IMPLICIT NONE
    !Parameters
    INTEGER,INTENT(out) :: nzsig,nxsig
    INTEGER,INTENT(in)  :: nz,nx,npml
    !
    nzsig=nz+2*npml
    nxsig=nx+2*npml
  END SUBROUTINE defnsig
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE defparamcpml2d(sigpmlza,sigpmlzb,sigpmlxa,sigpmlxb,vmin,vmax,npml,&
       dz,dx,nz,nx,nzsig,nxsig)
    !Parameters
    REAL(4),DIMENSION(nzsig)  ,INTENT(out) :: sigpmlza
    REAL(4),DIMENSION(nzsig-1),INTENT(out) :: sigpmlzb
    REAL(4),DIMENSION(nxsig)  ,INTENT(out) :: sigpmlxa
    REAL(4),DIMENSION(nxsig-1),INTENT(out) :: sigpmlxb
    REAL(4)                   ,INTENT(in)  :: vmin,vmax,dz,dx
    INTEGER                   ,INTENT(in)  :: nz,nx,npml,nzsig,nxsig
    !Local variables
    REAL(4) :: vm
    !
    vm=(vmin+vmax)/2._4
    !
    CALL defparamcpml(vm,npml,sigpmlza,dz,nz,nzsig)
    IF (nz.GT.1) THEN
       CALL defparamcpml(vm,npml,sigpmlzb,dz,nz-1,nzsig-1)
    ELSE
       sigpmlzb = sigpmlza ! Not used
    END IF
    !
    IF (nx.GT.1) THEN
       call defparamcpml(vm,npml,sigpmlxa,dx,nx,nxsig)
       call defparamcpml(vm,npml,sigpmlxb,dx,nx-1,nxsig-1)
    ELSE
       sigpmlxa = 0._4 ! Not used
       sigpmlxb = 0._4 ! Not used
    END IF
  END SUBROUTINE defparamcpml2d
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PURE SUBROUTINE defparamcpml(vm,npml,sig,dz,nz,nzsig)
    !Parameters
    REAL(4),DIMENSION(nzsig),INTENT(out) :: sig
    REAL(4)                 ,INTENT(in)  :: vm,dz
    INTEGER                 ,INTENT(in)  :: npml,nz,nzsig
    !Internal variables
    REAL(4)  :: L, R, vpml, val
    INTEGER  :: ind
    !
    L    = REAL(npml-1,4)*dz
    R    = 5.0_4
    vpml = 3._4 * vm / 2._4 / L * R
    !
    sig(:) = 0._4
    !
    DO ind = 1, npml
       val = vpml*(REAL(ind,4)/REAL(npml,4))**2
       sig(npml - ind + 1)  = val
       sig(nz + npml + ind) = val
    END DO
    !
  END SUBROUTINE defparamcpml
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE findiff
