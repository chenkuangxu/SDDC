    subroutine MATCHCYCREM
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: MLENGTH,count,I,TAG,M,IB,K,match_count,match_count_all,IP
    integer :: REQHANDLE(NPROC-1),ISTAT(MPI_STATUS_SIZE,NPROC-1),KLOC
    double precision :: XFA,YFA,XFLOCAL,YFLOCAL
    double precision :: epsX,epsY

    CALL MPI_ALLREDUCE(NCYCREM,MAXCYCREM,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)

    MAXCYCREM = MAXCYCREM + 1

    ALLOCATE(CYCRSBUF(MAXCYCREM*3))
    ALLOCATE(CYCRRBUF(MAXCYCREM*3,NPROC))
    ALLOCATE(CYCREM2PROC(NCYCREM),CYCREM2F_PROC(NCYCREM))

    CYCRSBUF = 0.d0
    CYCRRBUF = 0.d0

    MLENGTH = 3*MAXCYCREM

    count = 0

    DO I=1,NPROC
        IF (I.EQ.RANK+1) THEN
            CYCLE
        ELSE
            count = count + 1
            TAG = I*NPROC+RANK+1
            CALL MPI_IRECV(CYCRRBUF(1,I),MLENGTH,MPI_DOUBLE_PRECISION,I-1,&
            TAG,MPI_COMM_WORLD,REQHANDLE(count),IERROR)
        END IF
    END DO

    DO I=1,NPROC
        IF (I.EQ.RANK+1) THEN
            CYCLE
        ELSE
            M = 1
            DO IB=1,NCYCREM
                K = IBFCYCREM(IB)
                CYCRSBUF(M) = dble(K)
                CYCRSBUF(M+1) = ( XV(IF2V(K,1)) + XV(IF2V(K,2)) ) / 2
                CYCRSBUF(M+2) = ( YV(IF2V(K,1)) + YV(IF2V(K,2)) ) / 2
                M = M+3
            END DO
        END IF

        TAG = (RANK+1)*NPROC+I
        CALL MPI_SEND(CYCRSBUF,MLENGTH,MPI_DOUBLE_PRECISION,I-1,TAG,&
        MPI_COMM_WORLD,IERROR)
    END DO

    CALL MPI_WAITALL(NPROC-1,REQHANDLE,ISTAT,IERROR)

    DO IB=1,NCYCREM
        CYCREM2PROC(IB) = -1
        CYCREM2F_PROC(IB) = -1
    END DO

    match_count_all = 0
    DO I=1,NPROC
        match_count = 0
        IF (I.EQ.RANK+1) THEN
            CYCLE
        END IF
        M=1
        DO WHILE (CYCRRBUF(M,I).NE.0.d0)

            K = NINT(CYCRRBUF(M,I))
            XFA = CYCRRBUF(M+1,I)
            YFA = CYCRRBUF(M+2,I)
            M = M+3
            
            DO IP=1,NCYCREM
                KLOC = IBFCYCREM(IP)
                XFLOCAL = ( XV(IF2V(KLOC,1)) + XV(IF2V(KLOC,2)) ) / 2
                YFLOCAL = ( YV(IF2V(KLOC,1)) + YV(IF2V(KLOC,2)) ) / 2
                

                epsX = abs(abs(XFLOCAL - XFA) - DXCYCL)
                epsY = abs(abs(YFLOCAL - YFA) - DYCYCL)

                IF ( (epsX.lt.tolCYC.and.abs(YFA-YFLOCAL).lt.tolCYC).OR.& 
                (epsY.lt.tolCYC.and.abs(XFA-XFLOCAL).lt.tolCYC) ) THEN
                    CYCREM2PROC(IP) = I-1
                    CYCREM2F_PROC(IP) = K
                    match_count = match_count + 1
                END IF
            END DO
        END DO

        WRITE(*,400) 'processor',rank,'has matched',match_count,&
        'remote cyclic faces from rank',I-1
400     FORMAT(A,I6,1X,A,I8,1X,A,I6)
        match_count_all = match_count_all + match_count

    END DO  ! LOOP OVER FACES RECEIVED FROM PROCESSOR I

    end subroutine MATCHCYCREM