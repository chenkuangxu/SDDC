    subroutine MAPCYCLOC
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: nf,nf2,K,K2
    double precision :: XFIN,YFIN,XFOU,YFOU,epsX,epsY

    NPAIR = 0
    ALLOCATE(CYC2FPAIR(NCYCLOC/2,2))

    DO nf = 1,NCYCLOC

        K = IBFCYCLOC(nf)
        XFIN = ( XV(IF2V(K,1)) + XV(IF2V(K,2)) ) / 2
        YFIN = ( YV(IF2V(K,1)) + YV(IF2V(K,2)) ) / 2

        DO nf2 = 1,NCYCLOC

            IF (nf2.gt.nf) THEN
                K2 = IBFCYCLOC(nf2)
                XFOU = ( XV(IF2V(K2,1))+XV(IF2V(K2,2)) ) / 2
                YFOU = ( YV(IF2V(K2,1))+YV(IF2V(K2,2)) ) / 2
                epsX = abs(abs(XFIN-XFOU) - DXCYCL)
                epsY = abs(abs(YFIN-YFOU) - DYCYCL)

                IF ( (epsX.lt.tolCYC.and.abs(YFIN-YFOU).lt.tolCYC).OR.& 
                (epsY.lt.tolCYC.and.abs(XFIN-XFOU).lt.tolCYC) ) THEN
                    NPAIR = NPAIR + 1
                    CYC2FPAIR(NPAIR,1) = K
                    CYC2FPAIR(NPAIR,2) = K2
                END IF
            END IF
        END DO
    END DO

    IF (NPAIR/=NCYCLOC/2) THEN
        IF (rank.eq.0) PRINT *,'error in matching CYCLOC faces !!!'
        STOP
    END IF

    end subroutine MAPCYCLOC