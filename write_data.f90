    !----------------------------------------
    !Parallel Writing --- By Bin Zhang
    !              Jul.2013
    !----------------------------------------
    SUBROUTINE WRITE_ALL_DATA_PARALLEL

    use setup2d
    IMPLICIT NONE
    include 'mpif.h'

    !basic variables
    integer :: I, is, js, k
    integer, allocatable :: ICG(:)
    character(64) :: filename

    !MPI related variables
    integer :: fh, intsize, dbsize
    integer(kind=MPI_OFFSET_KIND) :: disp 
    integer :: status(MPI_STATUS_SIZE)
    integer :: etype, filetype, tmptype
    integer(kind=MPI_ADDRESS_KIND) :: lb, extent


    !get the size of int and real ...
    call MPI_TYPE_SIZE(MPI_INTEGER,intsize,ierror)
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,dbsize,ierror)


    !now get the displacement by cells
    if(rank .eq. 0) then
        disp = 0
    else
        disp = sum(PROC_NCELL(1:rank))
    end if


    !get the global cell no. 
    allocate(ICG(NCELL))
    do I=1,NCELL
        ICG(I) = IC2ICG(I)
    end do

    !open the file ...
    write(filename,'(A,I8.8,A)') 'restart.', iter, '.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierror)

    !write the header by root process ...
    if(rank.eq.0) then
        call MPI_FILE_WRITE(fh,iter,1,MPI_INTEGER,status,ierror)
        call MPI_FILE_WRITE(fh,ctime,1,MPI_DOUBLE_PRECISION,status,ierror)
    end if

    !write the global cell number first ...
    disp = intsize + dbsize + (intsize+dbsize*N*N*numv)*disp

    call MPI_TYPE_CONTIGUOUS(1,MPI_INTEGER,tmptype,ierror)
    lb=0; extent=intsize+dbsize*N*N*numv; etype=MPI_INTEGER
    call MPI_TYPE_CREATE_RESIZED(tmptype,lb,extent,filetype,ierror)
    call MPI_TYPE_COMMIT(filetype,ierror)

    call MPI_FILE_SET_VIEW(fh,disp,etype,filetype,"native",MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,ICG,NCELL,MPI_INTEGER,status,ierror)


    !then write the solution ...
    disp=disp+intsize

    call MPI_TYPE_CONTIGUOUS(N*N*numv,MPI_DOUBLE_PRECISION,tmptype,ierror)
    lb=0; extent=intsize+dbsize*N*N*numv; etype=MPI_DOUBLE_PRECISION
    call MPI_TYPE_CREATE_RESIZED(tmptype,lb,extent,filetype,ierror)
    call MPI_TYPE_COMMIT(filetype,ierror)

    call MPI_FILE_SET_VIEW(fh,disp,etype,filetype,"native",MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,Q,NCELL*N*N*numv,MPI_DOUBLE_PRECISION,status,ierror)


    !close the file ...
    call MPI_FILE_CLOSE(fh,ierror)

    !deallocate memory ...
    deallocate(ICG)

    END SUBROUTINE WRITE_ALL_DATA_PARALLEL
