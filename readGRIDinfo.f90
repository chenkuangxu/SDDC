    SUBROUTINE READ_CELL_DATA

    use setup2d
    IMPLICIT NONE
    include 'mpif.h'

    integer :: I,IC,K,ICTYPE,ICKEY,IOERR

    !go through the file to get total number of cells

    if(rank.eq.0) then
        WRITE(*,*) 'reading cells ...'    
        I = 0
        IOERR = 0
        open(20,file=NAMECEL)
        do while(IOERR.EQ.0)
            read(20,*,iostat=IOERR)
            I=I+1
        end do
        close(20)

        NCELLGLOB = I-1
        write(*,'(A,I6)') 'number of global cells=', NCELLGLOB
    endif

    call MPI_BCAST(NCELLGLOB,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)

    allocate(IVGCELLG(NCELLGLOB,4))

    !go through the file again to read data
    open(20,file=NAMECEL)
    do I=1,NCELLGLOB
      read(20,*) IC,(IVGCELLG(I,K),K=1,4),ICTYPE,ICKEY
    end do
    close(20)

    END SUBROUTINE READ_CELL_DATA

    ! read XVG YVG
    SUBROUTINE READ_VRT_DATA

    use setup2d
    IMPLICIT NONE
    include 'mpif.h'
    
    integer :: I,IV,IOERR
    
    !go through the file to get total number of vert
    
    if(rank.eq.0) then
        
        write(*,*) 'reading vertices ...'
    
        I = 0
        IOERR = 0
        open(10,file=NAMEVRT)
        do while(IOERR.EQ.0)
            read(10,*,iostat=IOERR)
            I=I+1
        end do
        close(10)
    
        NVERTGLOB = I-1
        write(*,'(A,I6)') ' number of global vertices = ',NVERTGLOB
    
    endif
    
    call MPI_BCAST(NVERTGLOB,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
    
    allocate(XVG(NVERTGLOB))
    allocate(YVG(NVERTGLOB))
    
    !go throught the file again to get the data
    
    open(10,file=NAMEVRT)
    
    do I=1, NVERTGLOB
        read(10,*) IV,XVG(I),YVG(I)
! 1000    FORMAT(I9,6X,3(D20.12,1X))
    end do
    
    close(10)
    
    END SUBROUTINE READ_VRT_DATA

    SUBROUTINE READ_VRT8_DATA
    use setup2d
    IMPLICIT NONE
    include 'mpif.h'
        
    integer :: ICG,K,IC
    double precision,allocatable :: XV8G(:,:),YV8G(:,:)
        
    ALLOCATE(XV8G(NCELLGLOB,8),YV8G(NCELLGLOB,8))
        
    open(12,FILE=NAMEVRT8)
    do ICG=1,NCELLGLOB
        do K=1,8
            read(12,*) XV8G(ICG,K),YV8G(ICG,K)
        end do
    end do
        
    ALLOCATE(xxf_cell(2,8,NCELL))
        
    do ICG=1,NCELLGLOB
        if (ICG2IC(ICG).NE.0) then
            IC = ICG2IC(ICG)
            do K=1,8
                xxf_cell(1,K,IC) = XV8G(ICG,K)
                xxf_cell(2,K,IC) = YV8G(ICG,K)
            end do
        end if
    end do
        
    DEALLOCATE(XV8G,YV8G)
        
    END SUBROUTINE READ_VRT8_DATA

    SUBROUTINE READ_METIS
    use setup2d
    IMPLICIT NONE
    
    include 'mpif.h'
    
    integer :: I
        
    allocate(IC2PROC(NCELLGLOB))
    allocate(PROC_NCELL(NPROC))
    
    NCELL = 0
    
    if (NPROC.NE.1) then
    
        open(13,FILE=METIS_CELL)
        do I = 1,NCELLGLOB
            READ(13,*) IC2PROC(I)
            if (IC2PROC(I).EQ.rank) NCELL = NCELL + 1
        end do
        CLOSE(13)
    
    else if (NPROC.EQ.1) then
    
        IC2PROC(:) = 0
    
        NCELL = NCELLGLOB
    
    end if
    
    print *,'processor',rank,'has',NCELL,'local cells'
    
    CALL MPI_ALLGATHER(NCELL,1,MPI_INTEGER,PROC_NCELL,1,MPI_INTEGER,MPI_COMM_WORLD,IERROR)
    
    END SUBROUTINE READ_METIS


    SUBROUTINE GLOB2LOCAL
    use setup2d
    IMPLICIT NONE
    
    include 'mpif.h'
    integer :: IC,IV,IVNEW,ICNEW,I
    integer :: IVG,IVL
    integer,allocatable :: IV2IVG_tmp(:)
    double precision,allocatable :: XV_tmp(:),YV_tmp(:)
    
    allocate(IVCELL(NCELL,4))
    allocate(ICG2IC(NCELLGLOB))
    allocate(IC2ICG(NCELL))
    allocate(IVG2IV(NVERTGLOB))
    
    allocate(XV_tmp(4*NCELL))
    allocate(YV_tmp(4*NCELL))
    allocate(IV2IVG_tmp(4*NCELL))
    
    ICG2IC = 0
    IVNEW  = 0
    ICNEW  = 0
    IVG2IV = 0
    
    do IC = 1,NCELLGLOB
            
        if (IC2PROC(IC).EQ.rank) then
              
            ICNEW = ICNEW + 1
            ICG2IC(IC) = ICNEW
            IC2ICG(ICNEW) = IC
    
            do IV = 1,4
                IVG = IVGCELLG(IC,IV)
                if (IVG2IV(IVG).EQ.0) then
                    IVNEW = IVNEW + 1
                    IVL   = IVNEW
                    IVG2IV(IVG) = IVL
                    IV2IVG_tmp(IVL) = IVG
                    XV_tmp(IVL) = XVG(IVG)
                    YV_tmp(IVL) = YVG(IVG)
                else
                    IVL = IVG2IV(IVG)
                end if
    
                IVCELL(ICNEW,IV) = IVL
    
            end do
    
        end if
    end do
    
    NVERT = IVNEW
    allocate(XV(NVERT))
    allocate(YV(NVERT))
    allocate(IV2IVG(NVERT))
    
    XV(1:NVERT) = XV_tmp(1:NVERT)
    YV(1:NVERT) = YV_tmp(1:NVERT)
    IV2IVG(1:NVERT) = IV2IVG_tmp(1:NVERT)
    
    deallocate(XVG)
    deallocate(YVG)
    deallocate(IC2PROC)
    deallocate(IVGCELLG)
        
    ! check a specific local vrt
    ! IV = 10
    ! print *,'this is rank',rank
    ! IVG = IV2IVG(IV)
    ! print *,'local index',IV
    ! print *,'global index',IVG
    ! print *,'coordinate',XV(IV),YV(IV)
    ! print *,IV2IVG
       
    END SUBROUTINE GLOB2LOCAL
    
    SUBROUTINE READ_BND_DATA
    use setup2d
    IMPLICIT NONE
    include 'mpif.h'
    
    INTEGER :: I,IOERR,IINLET,IOUTLET,ISYMP,IWALL,ICYCL,IFREE
    INTEGER :: IB,IPATCH,IREGION
    
    if (rank.eq.0) then
        write(*,*) 'reading boundaries ...'
        I = 0
        IOERR = 0
        open(11,FILE=NAMEBND)
        do while (IOERR.EQ.0)
            read(11,*,iostat=IOERR)
            I=I+1
        end do
        close(11)
    
        NBOUNDEFINED=I-1
        write(*,'(A,I6)') 'number of global boundaries =',NBOUNDEFINED
    end if
        
    call MPI_BCAST(NBOUNDEFINED,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
    
    ALLOCATE(IVGBOUN(NBOUNDEFINED,4),IBOUNTYPE(NBOUNDEFINED))
    
    IINLET = 0
    IOUTLET = 0
    ISYMP = 0
    IWALL = 0
    ICYCL = 0
    IFREE = 0
    
    open(11,FILE=NAMEBND)
    do I = 1,NBOUNDEFINED
        read(11,*) IB,IVGBOUN(I,1:2),IREGION,&
        IPATCH,IBOUNTYPE(I)
        if (IBOUNTYPE(I).EQ.'INLE') IINLET=IINLET+1
        if (IBOUNTYPE(I).EQ.'OUTL') IOUTLET=IOUTLET+1
        if (IBOUNTYPE(I).EQ.'SYMP') ISYMP=ISYMP+1
        if (IBOUNTYPE(I).EQ.'WALL') IWALL=IWALL+1
        if (IBOUNTYPE(I).EQ.'CYCL') ICYCL=ICYCL+1
        if (IBOUNTYPE(I).EQ.'FREE') IFREE=IFREE+1
    end do
    close(11)
        
    if (rank.eq.0) then
        write(*,*) '--------------------------------'
        write(*,*) 'INFO FROM READ BND DATA'
        write(*,*) 'Number of inlet boundary faces=',IINLET
        write(*,*) 'Number of outlet boundary faces=',IOUTLET
        write(*,*) 'Number of symmetry boundary faces=',ISYMP
        write(*,*) 'Number of wall boundary faces=',IWALL
        write(*,*) 'Number of cyclic faces=',ICYCL
        write(*,*) 'Number of free boundary faces=',IFREE
    endif
            
    END SUBROUTINE READ_BND_DATA