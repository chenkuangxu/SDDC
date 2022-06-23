    program SD2D_MHD
    use setup2d
    implicit none
    include 'mpif.h'

    double precision :: t_start,t_finish,t_resolution

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    NPROC = size

    if (rank.eq.0) t_start = mpi_wtime()

    call readinput

    call READ_CELL_DATA

    call READ_VRT_DATA

    call CHECK_MESH

    call READ_METIS

    call GLOB2LOCAL

    if (CURVE_WALL==1) call READ_VRT8_DATA

    call READ_BND_DATA

    call CONNECTIVITY

    call CONNECT_BDRY

    call MATCHPROC

    call MATCHCYCREM

    call GETIVCELL_PROC

    call GETIVCELL_CYC

    call MAPCYCLOC

    call INIT_SETUP

    ! call compresid

    ! call tecplotter2dsetup
    ! call tecplotter2d
    ! call tecpost_parallel

    ! call check_divergence_B

    ! call debug

    call iterations

    if (rank.eq.0) then
        t_finish     = mpi_wtime()
        t_resolution = mpi_wtick()
        print *, "Elapsed time ", t_finish - t_start, " seconds, &
        resolution ", t_resolution
    end if

    CALL MPI_FINALIZE(IERROR)


    end program SD2D_MHD