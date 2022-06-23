    SUBROUTINE GETIVCELL_CYC
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: I,K,TAG,SOURCE,MLENGTH,IC1,IC2,M,IV,IV2,DEST,KLOC
    double precision :: CRBUF(4*2,NCYCREM),CSBUF(4*2),XVR,YVR,&
    XVLOCAL,YVLOCAL,epsX,epsY
    integer :: REQHANDLE(NCYCREM),ISTAT(MPI_STATUS_SIZE,NCYCREM)

    ! MESSAGE LENGTH
    MLENGTH = 4*2

    DO I=1,NCYCREM
        K = IBFCYCREM(I)
        TAG = K
        SOURCE = CYCREM2PROC(I)
        CALL MPI_IRECV(CRBUF(1,I),MLENGTH,MPI_DOUBLE_PRECISION,&
        SOURCE,TAG,MPI_COMM_WORLD,REQHANDLE(I),IERROR)
    END DO

    DO I=1,NCYCREM
        K = IBFCYCREM(I)
        IC1 = IF2C(K,1)
        TAG = CYCREM2F_PROC(I)
        DEST = CYCREM2PROC(I)
        M = 1
        DO IV=1,4
            CSBUF(M)   = XV(IVCELL(IC1,IV))
            CSBUF(M+1) = YV(IVCELL(IC1,IV))
            M = M + 2
        END DO
        CALL MPI_SEND(CSBUF,MLENGTH,MPI_DOUBLE_PRECISION,&
        DEST,TAG,MPI_COMM_WORLD,IERROR)
    END DO

    CALL MPI_WAITALL(NCYCREM,REQHANDLE,ISTAT,IERROR)

    DO I=1,NCYCREM
        KLOC = IBFCYCREM(I)
        IC2  = IF2C(KLOC,2)

        M = 1
        DO IV=1,4
            XVR = CRBUF(M,I)
            YVR = CRBUF(M+1,I)
            M = M + 2

            DO IV2=1,2
                XVLOCAL = XV(IF2V(KLOC,IV2))
	            YVLOCAL = YV(IF2V(KLOC,IV2))

                epsX = abs( abs(XVLOCAL-XVR) - DXCYCL )
                epsY = abs( abs(YVLOCAL-YVR) - DYCYCL )

                IF ( (epsX.lt.tolCYC.and.abs(YVR-YVLOCAL).lt.tolCYC).OR.& 
                (epsY.lt.tolCYC.and.abs(XVR-XVLOCAL).lt.tolCYC) ) THEN
                    IVCELL(IC2,IV) = IF2V(KLOC,IV2)
                END IF
                
            END DO
        END DO
    END DO

    ! if (rank.eq.3) then
    !     I = 2
    !     K = IBFCYCREM(I)
    !     IC1 = IF2C(K,1)
    !     IC2 = IF2C(K,2)
    !     print *,IVCELL(IC1,:)
    !     print *,IVCELL(IC2,:)
    !     PRINT *,'NCELL=',NCELL,'NCELLTOT=',NCELLTOT
    ! end if

    END SUBROUTINE GETIVCELL_CYC