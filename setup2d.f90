    module setup2d
    implicit none

    ! MPI related variable
    integer :: size,rank,ierror,NPROC

    ! file names
    character(64) :: NAMECEL,NAMEVRT,NAMEVRT8,NAMEBND,METIS_CELL,NAMERESTART

    ! connectivity and grid variables
    ! global
    integer :: NCELLGLOB,NVERTGLOB,NBOUNDEFINED
    integer,allocatable :: IVGCELLG(:,:),IVGBOUN(:,:)
    character(8),allocatable :: IBOUNTYPE(:)
    double precision,allocatable :: XVG(:),YVG(:)
    integer,allocatable :: PROC_NCELL(:)

    ! local
    integer :: NCELL,NVERT,NFACE,NCELLTOT,NBOUN
    integer :: NINLET,NOUTLET,NSYMP,NWALL,NCYCL,NCYCLOC,NCYCREM,NPAIR,NFREE
    integer :: NPROCINT,NINTER
    integer :: MAXPROCINT,MAXCYCREM
    integer,allocatable :: IBFINL(:),IBFOUT(:),IBFSYMP(:),IBFWAL(:),IBFCYC(:),&
    IBFPROC(:),IF2IBPROC(:),IBFCYCLOC(:),IBFCYCREM(:),IF2IBCYCREM(:),IBFINTER(:),&
    IBFREE(:)
    integer,allocatable :: BOUNFACEINDEX(:)
    integer,allocatable :: IC2PROC(:)
    integer,allocatable :: ICG2IC(:),IC2ICG(:),IV2IVG(:),IVG2IV(:)
    integer,allocatable :: IVCELL(:,:),IC2F(:,:),IF2V(:,:),IF2C(:,:)
    double precision,allocatable :: XV(:),YV(:)
    integer,allocatable :: ICVERT(:),ICVSTA(:)
    integer,dimension(4,2) :: IVFACE
    integer,allocatable :: PROCSBUF(:),PROCRBUF(:,:),PROCINT2PROC(:),PROCINT2F_PROC(:),&
    CYCREM2PROC(:),CYCREM2F_PROC(:)
    integer,allocatable :: CYC2FPAIR(:,:)
    double precision,allocatable :: CYCRSBUF(:),CYCRRBUF(:,:)
    
    !---------------------coordinate related--------------------
    double precision,allocatable,dimension(:,:,:,:)::XXsolu,XXfluxi,XXfluxj
    double precision,allocatable,dimension(:,:,:) :: Jac
    double precision,allocatable,dimension(:,:,:,:,:) :: S1,S2
    double precision,allocatable,dimension(:,:,:,:) :: dmdxs,dmdxf1,dmdxf2
    integer :: CURVE_WALL
    double precision,allocatable,dimension(:,:,:) :: xxf_cell

    ! SD method
    integer :: N,numv
    double precision,allocatable,dimension(:) :: Xs,Xf
    double precision,allocatable,dimension(:,:) :: Lmat,Mmat

    integer,allocatable,dimension(:,:) :: iface2fp,jface2fp

    !-----------conserved variables and flux variables-------------
    double precision,allocatable,dimension(:,:,:,:) :: Q,F1,G2,Fv1,Gv2,resid
    double precision,allocatable,dimension(:,:,:,:) :: Qvfi,Qvfj
    double precision,allocatable,dimension(:,:,:,:,:) :: nablaQs,nablaQvfi,nablaQvfj
    double precision,allocatable,dimension(:,:,:) :: Qfl_p,Qfr_p,Qfl_c,Qfr_c,Qfl_p2,Qfr_p2,&
    Qfl_c2,Qfr_c2

    ! Boundary
    double precision :: DXCYCL,DYCYCL
    DOUBLE PRECISION,PARAMETER :: tolCYC = 1d-5
    double precision,allocatable :: dedxi(:,:,:),dedeta(:,:,:)

    !--------------------const parameter------------------------
    double precision,parameter :: pi = atan(1.d0)*4.d0
    double precision,parameter :: gam = 5.d0/3
    double precision,parameter :: Rair = 287.04d0
    double precision,parameter :: lambda = -2.d0/3.d0
    double precision,parameter :: prandt = 0.72d0
    double precision :: mu,eta0

    integer :: vismode

    double precision :: divBglob_report

    !-----------------------Time Marching-----------------------
    double precision :: ctime,dt
    integer :: k_stage,MAXITER,iter

    !-------------------divergence cleaning--------------------
    double precision :: c_h,alpha,umax,umax_g

    !-------------------artificial viscosity--------------------
    integer :: ARTIV
    double precision,allocatable :: muAV(:,:),muAVfi(:,:,:),muAVfj(:,:,:),SM_IND(:,:)
    double precision,allocatable :: Xs_low(:),Xf_low(:)
    double precision :: eigvmax,eigvmax_g,epsilon0
    double precision,allocatable :: plotAV(:,:)

    !---------------------postprocessing------------------------
    integer :: NRE,nwrite,plotnodes,plotcells
    integer,allocatable :: gnumcell(:,:,:),connec_c2n(:,:),&
    gnode_factor(:),gnumnode(:,:,:),isfacenod(:),&
    gfacenode(:,:),PROC_plotnodes(:),PROC_plotcells(:)
    double precision,allocatable :: plotX(:,:),plotQ(:,:)

    !-------------------------restart---------------------------
    integer :: restart,restart_ord

    !------------------------free-stream------------------------
    double precision :: rinf,pinf,uinf,vinf,Bxinf,Byinf

    end module setup2d