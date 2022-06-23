    subroutine MATCHPROC
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: MLENGTH,TAG
    integer :: I,M,K,IB,IFA,KLOC,IV1,IV2,IDIF11,IDIF21,IDIF13,IDIF23
    integer :: REQHANDLE(NPROC-1),ISTAT(MPI_STATUS_SIZE,NPROC-1)
    integer :: IVB1,IVB2,ICSTA1,ICEND1,ICSTA3,ICEND3,ICC,K1,K3,ICC1,ICC3
    integer :: count,match_count,match_count_all

    CALL MPI_ALLREDUCE(NPROCINT,MAXPROCINT,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)

    MAXPROCINT = MAXPROCINT + 1

    allocate(PROCSBUF(MAXPROCINT*3))
    allocate(PROCRBUF(MAXPROCINT*3,NPROC))
    allocate(PROCINT2PROC(NPROCINT),PROCINT2F_PROC(NPROCINT))

    PROCSBUF = 0
    PROCRBUF = 0

    MLENGTH = MAXPROCINT*3

    count = 0

    DO I=1,NPROC
        IF (I.EQ.RANK+1) THEN
            CYCLE
        ELSE
            count = count + 1
            TAG = I*NPROC+RANK+1
            CALL MPI_IRECV(PROCRBUF(1,I),MLENGTH,MPI_INTEGER,I-1,&
            TAG,MPI_COMM_WORLD,REQHANDLE(count),IERROR)
        END IF
    END DO

    DO I=1,NPROC
        IF (I.EQ.RANK+1) THEN
            CYCLE
        ELSE
            M = 1
            DO IB=1,NPROCINT
                K = IBFPROC(IB)
                PROCSBUF(M) = K 
                PROCSBUF(M+1) = IV2IVG(IF2V(K,1))
                PROCSBUF(M+2) = IV2IVG(IF2V(K,2))
                M=M+3
            END DO
        END IF
  
        TAG = (rank+1)*NPROC+I
        CALL MPI_SEND(PROCSBUF,MLENGTH,MPI_INTEGER,I-1,TAG,&
                      MPI_COMM_WORLD,IERROR)
    END DO

    CALL MPI_WAITALL(NPROC-1,REQHANDLE,ISTAT,IERROR)

    DO IB=1,NPROCINT
        PROCINT2PROC(IB) = -1
        PROCINT2F_PROC(IB) = -1
    END DO

    print *,'NPROCINT=',NPROCINT,'rank=',rank

    match_count_all = 0
    DO I=1,NPROC
        match_count = 0
        IF (I.EQ.RANK+1) THEN
            CYCLE
        END IF
        M=1
        DO 40 WHILE (PROCRBUF(M,I).NE.0)

            K = PROCRBUF(M,I)
            IVB1 = IVG2IV(PROCRBUF(M+1,I))
            IVB2 = IVG2IV(PROCRBUF(M+2,I))

            M = M+3

            IF ((IVB1.NE.0).AND.(IVB2.NE.0)) THEN

                ICSTA1=ICVSTA(IVB1)
                ICEND1=ICVSTA(IVB1+1)-1
                ICSTA3=ICVSTA(IVB2)
                ICEND3=ICVSTA(IVB2+1)-1
                ICC = 0

                ! Looping over all cells touching vertex IVB1
                DO K1=ICSTA1,ICEND1
                    ICC1=ICVERT(K1)
                ! Looping over all cells touching vertex IVB2
                DO K3=ICSTA3,ICEND3
                    ICC3=ICVERT(K3)

                IF (ICC1.EQ.ICC3) THEN
                    ICC=ICC1
                    ! Meaning we found the common cell
                    GOTO 15
                END IF
                !
                END DO
                END DO

                IF (ICC.eq.0) THEN
                    WRITE(*,*) 'No common cell ---- PROBLEM!!!'
                    CYCLE
                END IF 

15              DO IFA=1,4
                    IV1=IVCELL(ICC,IVFACE(IFA,1))
                    IV2=IVCELL(ICC,IVFACE(IFA,2))

                    IDIF11=IV1-IVB1
                    IDIF21=IV2-IVB1
                    IDIF13=IV1-IVB2
                    IDIF23=IV2-IVB2
                    IF ((IDIF11.EQ.0.OR.IDIF21.EQ.0) &
                    .AND. &
                   (IDIF13.EQ.0.OR.IDIF23.EQ.0)) GOTO 41
                END DO
                
                WRITE(*,*) 'Problem: Could not find common face'

41              KLOC=IC2F(ICC,IFA)
                IB = IF2IBPROC(KLOC)
                PROCINT2PROC(IB) = I-1
                PROCINT2F_PROC(IB) = K
                match_count= match_count+1

            END IF

40      CONTINUE

        WRITE(*,500) 'processor',rank,'has matched',match_count,&
        'processor interfaces from rank',I-1
500     FORMAT(A,I6,1X,A,I8,1X,A,I6)
        match_count_all = match_count_all + match_count 
    
    END DO
                
!     WRITE(*,200) 'processor',rank,'has matched',match_count_all,&
!     'processor interfaces'
! 200 FORMAT(A,I6,1X,A,I8,1X,A,I6)

    deallocate(PROCSBUF,PROCRBUF)

    ! IB = 1
    ! IF (RANK.EQ.0) THEN
    !     K = PROCINT2F_PROC(IB)
    !     PRINT *,RANK,PROCINT2PROC(IB),PROCINT2F_PROC(IB),IBFPROC(IB)
    ! END IF
    ! IF (RANK.EQ.2) THEN
    !     IB = IF2IBPROC(51)
    !     PRINT *,RANK,PROCINT2PROC(IB),PROCINT2F_PROC(IB),IBFPROC(IB)
    ! END IF

    end subroutine MATCHPROC