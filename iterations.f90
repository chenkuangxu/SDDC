    SUBROUTINE iterations
    
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: checkstat,ic,is,js,k,i
    double precision :: residnorm,rbuf(NPROC),residglob,divB,dbuf(NPROC),divBglob

    if (rank.eq.0) then
        open(11,file='resid.out',status='unknown')
        write(11,*) 'ITER    ctime     RESIDNORM   divB' 
        close(11)
        write(*,*) 'Start of iterations'
    end if

    checkstat = 0

    do while (checkstat/=1)
        iter = iter + 1

        if(k_stage==5) then
            CALL FivestageRK
        else if (k_stage==3) then
            CALL RK_S5P3
        else
            if (rank.eq.0) then
                print *,'k_stage can be only 3 or 5 now! Cannot run other values!'
            end if
            stop
        end if

        ! Calculate the residual norm
        if (mod(iter,10)==1) then 
            residnorm = 0.d0
            divB = 0.d0
            do ic=1,NCELL
                do js=1,N
                do is=1,N
                do k=1,9
                    residnorm = residnorm + (resid(k,is,js,ic))**2
                end do
                divB = divB + (nablaQs(6,1,is,js,ic)+nablaQs(7,2,is,js,ic))**2
                end do
                end do 
            end do
            ! Collect the residual norm from each processor
            CALL MPI_GATHER(residnorm, 1, MPI_DOUBLE_PRECISION, rbuf, 1, &
            MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)

            CALL MPI_GATHER(divB, 1, MPI_DOUBLE_PRECISION, dbuf, 1, &
            MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)

            if (rank.eq.0) then
                residglob = 0.d0
                divBglob = 0.d0
                do i=1,NPROC
                    residglob = residglob + rbuf(i)
                    divBglob = divBglob + dbuf(i)
                enddo
                residglob= sqrt(residglob/(NCELLGLOB*N*N*9))
                divBglob = sqrt(divBglob/(NCELLGLOB*N*N))
                divBglob_report = divBglob
                open(11,file='resid.out',position='append')
                write(11,*) iter, ctime, residglob, divBglob
                close(11)
                write (*,*) 'Rank:',rank,'Iteration:',  &
                iter, 'Residual: ', residglob,'divB: ',divBglob,&
                'eigvmax=',eigvmax_g,'mu_max=',epsilon0
            end if
          
        end if

        if(iter>=MAXITER) checkstat=1

        if (mod(iter,nwrite) == 1.or.iter.ge.MAXITER) then
            call tecplotter2dsetup
            call tecplotter2d
            call tecpost_parallel
            call WRITE_ALL_DATA_PARALLEL
        end if

    end do ! End of iteration loop
    
    END SUBROUTINE iterations
    
    !**********************************************
    subroutine FiveStageRk

    use setup2d

    IMPLICIT NONE
    include 'mpif.h'

    double precision,pointer,dimension(:,:,:,:):: qj0,qj2,qj3,resid0
    double precision :: old_time,rdt
    double precision :: a11,a21,a22,a23,a31,a32,a33,a41,&
    a42,a43,a51,a52,a53,a54,a55,a56

    allocate(qj0(9,N,N,NCELL),qj2(9,N,N,NCELL), &
    qj3(9,N,N,NCELL),resid0(9,N,N,NCELL))

    old_time = ctime
    qj0=Q

!   stage 1
    CALL compresid
    a11=0.39175222700392d0
    rdt = dt*a11
    Q =qj0+rdt*resid
    ctime = old_time+rdt

!   stage 2
    CALL compresid
    a21=0.44437049406734d0
    a22= 0.55562950593266d0
    a23=0.36841059262959d0
    Q= a21*qj0 +a22*Q  +a23*resid*dt
    qj2=Q
    rdt = rdt*a22+a23*dt
    ctime = old_time+rdt

!   stage 3
    CALL compresid
    a31=0.6201018513854d0
    a32=0.3798981486146d0
    a33=0.25189177424738d0
    Q=a31*qj0+a32*Q+a33*dt*resid
    qj3=Q
    rdt = rdt*a32+a33*dt
    ctime = old_time+rdt

!   stage 4
    CALL compresid
    a41=0.17807995410773d0
    a42=0.82192004589227d0
    a43=0.54497475021237d0
    Q=a41*qj0+a42*Q +a43*dt*resid
    resid0=resid
    rdt = rdt*a42+a43*dt
    ctime = old_time+rdt

!   stage 5
    CALL compresid
    a51=0.00683325884039d0
    a52=0.51723167208978d0
    a53=0.12759831133288d0
    a54=0.34833675773694d0
    a55=0.08460416338212d0
    a56= 0.22600748319395d0
    Q=a51*qj0+a52*qj2+a53*qj3 +a54*Q+a55*dt*resid0+a56*dt*resid
    ctime = old_time+dt

    deallocate(qj0,qj2,qj3,resid0)
    return
    end subroutine FiveStageRk

    subroutine RK_S5P3

    use setup2d
    implicit none
          
    double precision :: ctime0, dt1, dt2, dt3, dt4
    double precision, parameter :: a30=0.355909775063327d0
    double precision, parameter :: a32=0.644090224936674d0
    double precision, parameter :: a40=0.367933791638137d0
    double precision, parameter :: a43=0.632066208361863d0
    double precision, parameter :: a52=0.237593836598569d0
    double precision, parameter :: a54=0.762406163401431d0
    double precision, parameter :: b10=0.377268915331368d0
    double precision, parameter :: b21=0.377268915331368d0
    double precision, parameter :: b32=0.242995220537396d0
    double precision, parameter :: b43=0.238458932846290d0
    double precision, parameter :: b54=0.287632146308408d0
    double precision,allocatable,dimension(:,:,:,:) :: Qk0,Qka
    
    allocate(Qk0(9,N,N,NCELL),Qka(9,N,N,NCELL))
          
    Qk0 = Q
    ctime0 = ctime
          
    !stage 1
    call compresid
    Q = Q + b10*dt*resid
    dt1 = b10*dt
    ctime = ctime0 + dt1
          
    !stage 2
    call compresid
    Q = Q + b21*dt*resid
    dt2 = dt1 + b21*dt
    ctime = ctime0 + dt2
          
    Qka = Q
          
    !stage 3
    call compresid
    Q = a30*Qk0 + a32*Q + b32*dt*resid
    dt3 = a32*dt2 + b32*dt
    ctime = ctime0 + dt3
          
    !stage 4
    call compresid
    Q = a40*Qk0 + a43*Q + b43*dt*resid
    dt4 = a43*dt3 + b43*dt
    ctime = ctime0 + dt4
          
    !stage 5
    call compresid
    Q = a52*Qka + a54*Q + b54*dt*resid
    ctime = ctime0 + dt
          
    end subroutine RK_S5P3