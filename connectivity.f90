    subroutine CONNECTIVITY
    use setup2d
    implicit none
    include 'mpif.h'

    INTEGER :: IC,IV,K,I,IC1,IC2,IV1,IV2,IVBASE
    INTEGER :: IC1V1,IC1V2,IC1V3,IC1V4,IFFA,ICC
    INTEGER :: IC2V1,IC2V2,IC2V3,IC2V4,IVV1,IVV2
    INTEGER :: IFA, IFACE, IFACEOLD, IBOUN,NF
    INTEGER :: IPROD11,IPROD12,IPROD21,IPROD22,L
    INTEGER :: KSTABASE, KENDBASE

    INTEGER, DIMENSION (:), ALLOCATABLE :: IFILL
    INTEGER, DIMENSION (:), ALLOCATABLE :: IDUMMY
    INTEGER :: MAXFACE
    INTEGER,ALLOCATABLE :: IF2C_tmp(:,:),IF2V_tmp(:,:)

    INTEGER,DIMENSION(2) :: IVtest
    INTEGER,DIMENSION(2) :: ICtest

    MAXFACE = NCELL*4

    ALLOCATE(IC2F(NCELL,4))
    ALLOCATE(IF2C_tmp(MAXFACE,2))
    ALLOCATE(IF2V_tmp(MAXFACE,2))
    ALLOCATE(IFILL(NVERT),IDUMMY(NVERT))
    ALLOCATE(ICVSTA(NVERT+1))

!---Calculate the cells surrounding each vertex IV
    DO IV=1,NVERT
        IFILL(IV)=0
    END DO

!---Create IFILL for knowing number of cells which the vertex is
!associated with
    DO IC=1,NCELL
        DO K=1,4
            IV = IVCELL(IC,K)
            IFILL(IV)=IFILL(IV)+1
        END DO
    END DO

!--------checking number
    DO IV=1,NVERT
        IF(IFILL(IV).EQ.0) THEN
            write (*,*) 'IFILL =', &
            IFILL(IV), 'at vertex', IV
            STOP
        END IF
    END DO

    K=1
    DO IV=1,NVERT
        ICVSTA(IV)=K
        IDUMMY(IV)=ICVSTA(IV)
        K=K+IFILL(IV)
    END DO

!----------
    ICVSTA(NVERT+1)=ICVSTA(NVERT)+IFILL(NVERT)
    ALLOCATE(ICVERT(ICVSTA(NVERT+1)))

    DO IC=1,NCELL
        DO K=1,4
            IV=IVCELL(IC,K)
            L=IDUMMY(IV)
            ICVERT(L)=IC
            IDUMMY(IV)=IDUMMY(IV)+1
        END DO
    END DO

    !---Variable IVFACE(IFACE,1-4)
    IVFACE(1,1)=1
    IVFACE(1,2)=2
!
    IVFACE(2,1)=2
    IVFACE(2,2)=3
!
    IVFACE(3,1)=3
    IVFACE(3,2)=4
!
    IVFACE(4,1)=4
    IVFACE(4,2)=1

    DO IC=1,NCELL
        DO IFA=1,4
            IC2F(IC,IFA)=0
        END DO
    END DO

    DO IFACE=1,MAXFACE
        IF2C_tmp(IFACE,1) = 0
        IF2C_tmp(IFACE,2) = 0
    END DO

    IFACE = 0

    DO 100 IC=1,NCELL
        DO 110 IFA=1,4
            IFACEOLD=IFACE
            IF (IC2F(IC,IFA).NE.0) GOTO 110
            IVBASE=IVCELL(IC,IVFACE(IFA,1))
            IV1   =IVCELL(IC,IVFACE(IFA,2))
            KSTABASE=ICVSTA(IVBASE)
            KENDBASE=ICVSTA(IVBASE+1)-1 
            DO 120 K=KSTABASE,KENDBASE
            ICC=ICVERT(K)

            IF (ICC.EQ.IC) GOTO 120
            DO 130 IFFA=1,4
            IVV1=IVCELL(ICC,IVFACE(IFFA,1))
            IVV2=IVCELL(ICC,IVFACE(IFFA,2))

            IPROD11=(IVBASE-IVV1)
            IPROD12=(IVBASE-IVV2)

            IPROD21=(IV1-IVV1)
            IPROD22=(IV1-IVV2)

            IF ((IPROD11.EQ.0.AND.IPROD22.EQ.0).OR. &
            (IPROD12.EQ.0.AND.IPROD21.EQ.0)) THEN

                IFACE=IFACE+1
                IC2F(ICC,IFFA)=IFACE

!-- Face to vertex structures
                IF2V_tmp(IFACE,1) = IVBASE
                IF2V_tmp(IFACE,2) = IV1

                IF2C_tmp(IFACE,1) = IC
                IF2C_tmp(IFACE,2) = ICC

                GOTO 115

            END IF

130         CONTINUE
120         CONTINUE

            IF (IF2C_tmp(IFACEOLD+1,1).EQ.0) THEN
            IF2C_tmp(IFACEOLD+1,1)=IC
            IF2C_tmp(IFACEOLD+1,2)=0
            IFACE=IFACEOLD+1
!----------- for near boundary faces-------------
            IF2V_tmp(IFACE,1) = IVBASE
            IF2V_tmp(IFACE,2) = IV1
            END IF

115         IC2F(IC,IFA)=IFACE
110     CONTINUE
100 CONTINUE

    NFACE=IFACE

    DEALLOCATE(IFILL,IDUMMY)
    ALLOCATE(IF2C(NFACE,4))
    ALLOCATE(IF2V(NFACE,2))
    IF2C(1:NFACE,1:2)=IF2C_tmp(1:NFACE,1:2)
    IF2V(1:NFACE,:)=IF2V_tmp(1:NFACE,:)
    DEALLOCATE(IF2C_tmp,IF2V_tmp)

    IF2C(:,3) = 0
    IF2C(:,4) = 0


    DO NF=1,NFACE

        IC1=IF2C(NF,1)
        IC2=IF2C(NF,2)
        IV1= IF2V(NF,1)
        IV2= IF2V(NF,2)
        IC1V1= IVCELL(IC1,1)
        IC1V2= IVCELL(IC1,2)
        IC1V3= IVCELL(IC1,3)
        IC1V4= IVCELL(IC1,4)

        ! Face vertex sequence is given by the left element
        IF (IC1V1.EQ.IV1.AND.IC1V2.EQ.IV2) THEN
            IF2C(NF,3)=1
        ELSE IF (IC1V2.EQ.IV1.AND.IC1V3.EQ.IV2) THEN
            IF2C(NF,3)=2
        ELSE IF (IC1V3.EQ.IV1.AND.IC1V4.EQ.IV2) THEN
            IF2C(NF,3)=3
        ELSE IF (IC1V4.EQ.IV1.AND.IC1V1.EQ.IV2) THEN
            IF2C(NF,3)=4
        ELSE
            PRINT *,'Something wrong in connectivity!'
            PRINT *,'NF=',NF
            PRINT *,IC1V1,IC1V2,IC1V3,IC1V4
            PRINT *,IV1,IV2
            stop
        END IF


        ! Since all element vertices are arranged counter-clockwise
        ! Right element Face2vertex is the opposite
        IF (IC2.NE.0) THEN
            IC2V1= IVCELL(IC2,1)
            IC2V2= IVCELL(IC2,2)
            IC2V3= IVCELL(IC2,3)
            IC2V4= IVCELL(IC2,4)
            IF (IC2V2.EQ.IV1.AND.IC2V1.EQ.IV2) THEN
                IF2C(NF,4)=1
            ELSE IF (IC2V3.EQ.IV1.AND.IC2V2.EQ.IV2) THEN
                IF2C(NF,4)=2
            ELSE IF (IC2V4.EQ.IV1.AND.IC2V3.EQ.IV2) THEN
                IF2C(NF,4)=3
            ELSE IF (IC2V1.EQ.IV1.AND.IC2V4.EQ.IV2) THEN
                IF2C(NF,4)=4
            ELSE
                PRINT *,'Something wrong in connectivity!'
                stop
            END IF
        END IF

    END DO

    ! ! examine the connectivity of a particular face
    ! IFACE = 10
    ! do I = 1,2
    !     IVtest(I) = IF2V(IFACE,I)
    ! end do
    ! write(*,*) 'processor',rank,'the two nodal points associated with this face are',&
    ! IVtest(:),'(local)'
    ! write(*,*) 'processor',rank,'the two nodal points associated with this face are',&
    ! IV2IVG(IVtest(1)),IV2IVG(IVtest(2)),'(global)'
    ! do I = 1,2
    !     ICtest(I) = IF2C(IFACE,I)
    ! end do
    ! write(*,*) 'processor',rank,'the two adjacent cells are ',ICtest(:),'(local)'
    ! write(*,*) 'processor',rank,'the two adjacent cells are ',&
    ! IC2ICG(ICtest(1)),IC2ICG(ICtest(2)),'(global)'

    end subroutine CONNECTIVITY


    subroutine CHECK_MESH
    use setup2d
    implicit none
    include 'mpif.h'
    
    ! Check whether every cell connectivity is counter-clockwise
    ! If not, reverse the local cell connectivity
    
    integer :: IC,IV,IVF(4),K,checkmark,IFA
    double precision :: xxs(2,4),facepro(4),sumpro,x1,x2,y1,y2
    
    !---Variable IVFACE(IFACE,1-4)
    IVFACE(1,1)=1
    IVFACE(1,2)=2
    !
    IVFACE(2,1)=2
    IVFACE(2,2)=3
    !
    IVFACE(3,1)=3
    IVFACE(3,2)=4
    !
    IVFACE(4,1)=4
    IVFACE(4,2)=1
        
    checkmark = 1
    
    DO IC=1,NCELLGLOB
        DO K=1,4
            IV=IVGCELLG(IC,K)
            xxs(1,K) = XVG(IV)
            xxs(2,K) = YVG(IV)
        END DO
            
        DO IFA = 1,4
            x1 = xxs(1,IVFACE(IFA,1))
            x2 = xxs(1,IVFACE(IFA,2))
            y1 = xxs(2,IVFACE(IFA,1))
            y2 = xxs(2,IVFACE(IFA,2))
            facepro(IFA) = (x2 - x1) * (y1 + y2)
        END DO
    
        sumpro = facepro(1)+facepro(2)+facepro(3)+facepro(4)
        
        if (sumpro.gt.0) then
            IVF(1:4) = IVGCELLG(IC,1:4)
            IVGCELLG(IC,1) = IVF(2)
            IVGCELLG(IC,2) = IVF(1)
            IVGCELLG(IC,3) = IVF(4)
            IVGCELLG(IC,4) = IVF(3)
    
            checkmark = 0
        end if
    
    END DO
    
    IF (checkmark==1.AND.RANK.EQ.0) THEN
        PRINT *,'Local cell connectivity is ALL counter-clockwise! Correct!'
    ELSE IF (RANK.EQ.0) THEN
        PRINT *,'Some cell connectivity is reversed to make all cell connectivity counter-clockwise!'
    END IF
    
    end subroutine CHECK_MESH