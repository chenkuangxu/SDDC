    subroutine CONNECT_BDRY
    use setup2d
    implicit none
    include 'mpif.h'

    INTEGER :: IB,IB2,K,IVB1,IVB3,ICSTA1,ICEND1,ICSTA3,ICEND3
    INTEGER :: ICC,ICC1,ICC3,K1,K3,IFA,IV1,IV2,IDIF11,IDIF21,IDIF13,IDIF23
    INTEGER :: IC1,IC2
    INTEGER,ALLOCATABLE :: IBFPROC_tmp(:),IBFCYCLOC_tmp(:),IBFCYCREM_tmp(:),IVCELL_tmp(:,:)
    DOUBLE PRECISION :: XFIN,YFIN,XFOU,YFOU,epsX,epsY

    ALLOCATE(BOUNFACEINDEX(NFACE))
    DO K = 1,NFACE
        BOUNFACEINDEX(K) = 0
    END DO

! First we need to find the NINLET,NOUTLET,NSYMP,NWALL,NCYCL,NFREE
    NINLET=0
    NOUTLET=0
    NSYMP=0
    NWALL=0
    NCYCL=0
    NFREE=0

    DO 10 IB=1,NBOUNDEFINED

        K=0

        IVB1=IVG2IV(IVGBOUN(IB,1))
        IVB3=IVG2IV(IVGBOUN(IB,2)) 

        IF ((IVB1.EQ.0).OR.(IVB3.EQ.0)) THEN
            GOTO 10
        END IF

        ICSTA1=ICVSTA(IVB1)
        ICEND1=ICVSTA(IVB1+1)-1
        ICSTA3=ICVSTA(IVB3)
        ICEND3=ICVSTA(IVB3+1)-1
        ICC=0

! Looping over all cells touching vertex IVB1
        DO K1=ICSTA1,ICEND1
            ICC1=ICVERT(K1)

! Looping over all cells touching vertex IVB3
            DO K3=ICSTA3,ICEND3
                ICC3=ICVERT(K3)
                
                IF (ICC1.EQ.ICC3) ICC=ICC1
! Meaning we found the cell at the boundary, cell ICC
            END DO
        END DO

        IF (ICC.EQ.0) THEN
            WRITE(*,*) 'Does not belong to processor!!'
            GOTO 10
        END IF

        DO IFA=1,4
            IV1=IVCELL(ICC,IVFACE(IFA,1))
            IV2=IVCELL(ICC,IVFACE(IFA,2))

            IDIF11=IV1-IVB1
            IDIF21=IV2-IVB1
            IDIF13=IV1-IVB3
            IDIF23=IV2-IVB3 

            IF ((IDIF11.EQ.0.OR.IDIF21.EQ.0).AND. &
            (IDIF13.EQ.0.OR.IDIF23.EQ.0)) GOTO 41
        END DO

        WRITE(*,*) 'Should not be here'

41      K=IC2F(ICC,IFA)
        BOUNFACEINDEX(K)=1

        IF (K.EQ.0) THEN
            WRITE(*,*) 'Problem with the attachement of boundaries'
            WRITE(*,*) 'IB =',IB,'ICC =',ICC,'rank=',rank
            WRITE(*,*) 'ICC=',ICC,'IFA=',IFA
            WRITE(*,*) 'IVB1=',IVB1,'IVB3=',IVB3
            WRITE(*,*) 'Execution paused !'
! Meaning there is one face defined in NAMEBND
! that was not detected to be a bdy face
        ENDIF

        if (IBOUNTYPE(IB).EQ.'INLE') then
            NINLET=NINLET+1
        end if

        if (IBOUNTYPE(IB).EQ.'OUTL') then
            NOUTLET=NOUTLET+1
        end if

        if (IBOUNTYPE(IB).EQ.'SYMP') then
            NSYMP=NSYMP+1
        end if

        if (IBOUNTYPE(IB).EQ.'WALL') then
            NWALL=NWALL+1
        end if

        if (IBOUNTYPE(IB).EQ.'CYCL') then
            NCYCL=NCYCL+1
        end if

        if (IBOUNTYPE(IB).EQ.'FREE') then
            NFREE=NFREE+1
        end if
                
10  END DO

    ALLOCATE(IBFINL(NINLET),IBFOUT(NOUTLET),IBFSYMP(NSYMP),IBFWAL(NWALL),&
    IBFCYC(NCYCL),IBFREE(NFREE))

! ====================================================================================
    NBOUN = 0
    NINLET=0
    NOUTLET=0
    NSYMP=0
    NWALL=0
    NCYCL=0
    NFREE=0

    DO 20 IB=1,NBOUNDEFINED

        K=0

        IVB1=IVG2IV(IVGBOUN(IB,1))
        IVB3=IVG2IV(IVGBOUN(IB,2)) 

        IF ((IVB1.EQ.0).OR.(IVB3.EQ.0)) THEN
            GOTO 20
        END IF

        ICSTA1=ICVSTA(IVB1)
        ICEND1=ICVSTA(IVB1+1)-1
        ICSTA3=ICVSTA(IVB3)
        ICEND3=ICVSTA(IVB3+1)-1
        ICC=0

! Looping over all cells touching vertex IVB1
        DO K1=ICSTA1,ICEND1
            ICC1=ICVERT(K1)

! Looping over all cells touching vertex IVB3
            DO K3=ICSTA3,ICEND3
                ICC3=ICVERT(K3)
                
                IF (ICC1.EQ.ICC3) ICC=ICC1
! Meaning we found the cell at the boundary, cell ICC
            END DO
        END DO

        IF (ICC.EQ.0) THEN
            WRITE(*,*) 'Does not belong to processor!!'
            GOTO 20
        END IF

        DO IFA=1,4
            IV1=IVCELL(ICC,IVFACE(IFA,1))
            IV2=IVCELL(ICC,IVFACE(IFA,2))

            IDIF11=IV1-IVB1
            IDIF21=IV2-IVB1
            IDIF13=IV1-IVB3
            IDIF23=IV2-IVB3 

            IF ((IDIF11.EQ.0.OR.IDIF21.EQ.0).AND. &
            (IDIF13.EQ.0.OR.IDIF23.EQ.0)) GOTO 21
        END DO

        WRITE(*,*) 'Should not be here'

21      K=IC2F(ICC,IFA)
        BOUNFACEINDEX(K)=1

        IF (K.EQ.0) THEN
            WRITE(*,*) 'Problem with the attachement of boundaries'
            WRITE(*,*) 'IB =',IB,'ICC =',ICC
            WRITE(*,*) 'Execution paused !'
! Meaning there is one face defined in NAMEBND
! that was not detected to be a bdy face
        ENDIF

        if (IBOUNTYPE(IB).EQ.'INLE') then
            NINLET=NINLET+1
            NBOUN=NBOUN+1
            IBFINL(NINLET)=K
        end if

        if (IBOUNTYPE(IB).EQ.'OUTL') then
            NOUTLET=NOUTLET+1
            NBOUN=NBOUN+1
            IBFOUT(NOUTLET)=K
        end if

        if (IBOUNTYPE(IB).EQ.'SYMP') then
            NSYMP=NSYMP+1
            NBOUN=NBOUN+1
            IBFSYMP(NSYMP)=K
        end if

        if (IBOUNTYPE(IB).EQ.'WALL') then
            NWALL=NWALL+1
            NBOUN=NBOUN+1
            IBFWAL(NWALL)=K
        end if

        if (IBOUNTYPE(IB).EQ.'CYCL') then
            NCYCL=NCYCL+1
            NBOUN=NBOUN+1
            IBFCYC(NCYCL)=K
        end if

        if (IBOUNTYPE(IB).EQ.'FREE') then
            NFREE=NFREE+1
            NBOUN=NBOUN+1
            IBFREE(NFREE)=K
        end if
                
20  END DO

    NCELLTOT = NCELL

! ====================================================================================
    NINTER = 0
    DO K=1,NFACE
        IC2=IF2C(K,2)
        IF (IC2.NE.0) THEN
            NINTER = NINTER+1
        END IF
    END DO
    ALLOCATE(IBFINTER(NINTER))

    NINTER = 0
    DO K=1,NFACE
        IC2=IF2C(K,2)
        IF (IC2.NE.0) THEN
            NINTER = NINTER+1
            IBFINTER(NINTER) = K
        END IF
    END DO

! ====================================================================================
    
    ALLOCATE(IBFPROC_tmp(NFACE))
    ALLOCATE(IF2IBPROC(NFACE))

    IF2IBPROC(:) = 0

    NPROCINT = 0

    DO K=1,NFACE
        IC1=IF2C(K,1)
        IC2=IF2C(K,2)

        IF (IC2.EQ.0.AND.BOUNFACEINDEX(K).EQ.0) THEN
            NPROCINT=NPROCINT+1
            IBFPROC_tmp(NPROCINT)=K
            IF2IBPROC(K)=NPROCINT

            NCELLTOT=NCELLTOT + 1
            IF2C(K,2)=NCELLTOT 
        END IF

    END DO

    allocate(IBFPROC(NPROCINT))
    IBFPROC(1:NPROCINT) = IBFPROC_tmp(1:NPROCINT)
    deallocate(IBFPROC_tmp)

! ====================================================================================
    ALLOCATE(IBFCYCLOC_tmp(NFACE),IBFCYCREM_tmp(NFACE))
    ALLOCATE(IF2IBCYCREM(NFACE))
    IF2IBCYCREM(:) = 0
    NCYCLOC = 0
    NCYCREM = 0

    DO 160 IB=1,NCYCL
        K = IBFCYC(IB)
        XFIN = (XV(IF2V(K,1)) + XV(IF2V(K,2))) / 2
        YFIN = (YV(IF2V(K,1)) + YV(IF2V(K,2))) / 2
        
        DO IB2=1,NCYCL
            K = IBFCYC(IB2)
            XFOU = (XV(IF2V(K,1)) + XV(IF2V(K,2))) / 2
            YFOU = (YV(IF2V(K,1)) + YV(IF2V(K,2))) / 2
            epsX = abs(abs(XFIN-XFOU) - DXCYCL)
            epsY = abs(abs(YFIN-YFOU) - DYCYCL)

            IF ( (abs(XFIN-XFOU).lt.tolCYC.AND.epsY.lt.tolCYC).OR.&
            (abs(YFIN-YFOU).lt.tolCYC.AND.epsX.lt.tolCYC) ) THEN
                NCYCLOC=NCYCLOC+1
                IBFCYCLOC_tmp(NCYCLOC)=IBFCYC(IB)
                GOTO 160
            END IF
        END DO

        NCYCREM = NCYCREM + 1
        IBFCYCREM_tmp(NCYCREM) = IBFCYC(IB)
        IF2IBCYCREM(IBFCYC(IB)) = IB

        NCELLTOT = NCELLTOT + 1
        IF2C(IBFCYC(IB),2) = NCELLTOT

160 END DO

    allocate(IBFCYCLOC(NCYCLOC),IBFCYCREM(NCYCREM))
    IBFCYCLOC(1:NCYCLOC) = IBFCYCLOC_tmp(1:NCYCLOC)
    deallocate(IBFCYCLOC_tmp)
    IBFCYCREM(1:NCYCREM) = IBFCYCREM_tmp(1:NCYCREM)
    deallocate(IBFCYCREM_tmp)

    allocate(IVCELL_tmp(NCELL,4))
    IVCELL_tmp = IVCELL
    deallocate(IVCELL)
    allocate(IVCELL(NCELLTOT,4))
    IVCELL = 0
    IVCELL(1:NCELL,:) = IVCELL_tmp(1:NCELL,:)
    deallocate(IVCELL_tmp)

    WRITE(*,*) '--------------------------------'
    WRITE(*,*) 'Info from connectivity'
    WRITE(*,*) '--------------------------------'
    WRITE(*,*) 'RANK =',rank
    WRITE(*,*) 'Number of cells (NCELL) =',NCELL
    WRITE(*,*) 'Total number of cells (NCELLTOT) &
                  (internal+procint+ncycrem)=',NCELLTOT
    WRITE(*,*) 'Total number of faces: ',NFACE
    WRITE(*,*) 'Number of Interfaces=',NINTER
    WRITE(*,*) 'Number of total boundary faces=',NBOUN
    WRITE(*,*) 'Boundary faces attached'
    WRITE(*,*) 'Inlet faces attached=',NINLET
    WRITE(*,*) 'Outlet faces attached=',NOUTLET
    WRITE(*,*) 'Symmetry faces attached=',NSYMP
    WRITE(*,*) 'Wall faces attached=',NWALL
    WRITE(*,*) 'Local Cyclic faces=',NCYCLOC
	WRITE(*,*) 'Remote Cyclic faces=',NCYCREM
    WRITE(*,*) 'Processor interfaces attached=',NPROCINT
    WRITE(*,*) 'Free boundary faces=',NFREE
    WRITE(*,*) '--------------------------------'


    DEALLOCATE(IBOUNTYPE,BOUNFACEINDEX)

    ! IF (RANK.EQ.1) THEN
    !     PRINT *,'------'
    !     PRINT *,IBFCYCLOC(:)
    !     PRINT *,IV2IVG(IF2V(39,1)),IV2IVG(IF2V(39,2))
    ! END IF

    end subroutine CONNECT_BDRY