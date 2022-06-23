    SUBROUTINE GETIVCELL_PROC

    use setup2d
    IMPLICIT NONE

    include 'mpif.h'
        
    ! F2PROC(i) is the processor number connected to boundary face i
    ! PROCINT2F_PROC(IB) is the cell number IC2 of processor F2PROC(i) connected to face IB

    INTEGER :: I,K,IC1,IC2,M,K1,ks,js,is,IV
    INTEGER :: TAG,DEST,SOURCE,MLENGTH
    INTEGER :: REQHANDLE(NPROCINT),ISTAT(MPI_STATUS_SIZE,NPROCINT)

    INTEGER :: VRBUF(4,NPROCINT),VSBUF(4,NPROCINT)
    ! Message length
    MLENGTH = 4

    ! ----------------------------------------------------------------------

    ! POST NON-BLOCKING RECEIVES

    DO I=1,NPROCINT ! Loop over the faces with type '    PROC'
        
        K = IBFPROC(I) ! Get face index
        ! CELL CONTAINING THE INFORMATION TO RECEIVE
        ! TAG IS CELL IC2 of Current Processor
        TAG = K
        ! PROCESSOR TO RECEIVE FROM
        SOURCE = PROCINT2PROC(I)
        CALL MPI_IRECV(VRBUF(1,I),MLENGTH,MPI_INTEGER,&
        SOURCE,TAG,MPI_COMM_WORLD,REQHANDLE(I),IERROR)
    END DO


    ! SEND THE DATA

    DO I=1,NPROCINT ! Loop over the faces with type '    PROC'
        
        K = IBFPROC(I) ! Get face index
        ! CELL CONTAINING THE INFORMATION TO SEND
        IC1=IF2C(K,1)
        ! TAG IS FACE K of PROCESSOR F2PROC(I)
        TAG = PROCINT2F_PROC(I)
        ! PROCESSOR TO SEND TO
        DEST = PROCINT2PROC(I)
        ! PACK THE SENDING BUFFER

        M=1
        DO IV=1,4
            VSBUF(M,I) = IV2IVG(IVCELL(IC1,IV))
            M = M+1
        END DO

        CALL MPI_SEND(VSBUF(1,I),MLENGTH,MPI_INTEGER,DEST,TAG,&
        MPI_COMM_WORLD,IERROR)

    END DO

    ! WAIT FOR THE COMPLETION OF THE NON_BLOCKING RECEIVES

    CALL MPI_WAITALL(NPROCINT,REQHANDLE,ISTAT,IERROR)

    ! UNPACK THE RECEIVED BUFFER

    DO I=1,NPROCINT
        
        K = IBFPROC(I)
        IC2 = IF2C(K,2)

        M=1
        DO IV=1,4
            IVCELL(IC2,IV) = IVG2IV(VRBUF(M,I))
            ! NOTE THAT SOME OF THE IVCELL WILL BE ZERO
            M = M+1
        END DO
    END DO
  
    END SUBROUTINE GETIVCELL_PROC