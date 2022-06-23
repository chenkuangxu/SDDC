    SUBROUTINE READ_ALL_DATA_BINARY_MPI

    use setup2d
    IMPLICIT NONE
    INCLUDE 'mpif.h'
  
    integer :: I,is,js,k,ICG,IC,fh
    integer :: status(MPI_STATUS_SIZE)
    double precision,allocatable :: tmp(:,:,:)
  
    open(newunit=fh,file=NAMERESTART,access='stream')
  
    read(fh) iter
    read(fh) ctime
  
    !now read cell no. and solutions...
    allocate(tmp(numv,N,N))
  
    do I=1, NCELLGLOB
     
        read(fh) ICG
        read(fh) tmp(1:numv,1:N,1:N)
  
        IC=ICG2IC(ICG)
  
        if(IC.ne.0) Q(1:numv,1:N,1:N,IC)=tmp(1:numv,1:N,1:N)
    end do
  
    close(fh)
  
    deallocate(tmp)
  
    END SUBROUTINE READ_ALL_DATA_BINARY_MPI

    SUBROUTINE READ_ALL_DATA_BINARY

    use setup2d
    IMPLICIT NONE
    INCLUDE 'mpif.h'
      
    integer :: I,is,js,k,ICG,IC,fh
    integer :: status(MPI_STATUS_SIZE)
    double precision,allocatable :: tmp(:,:,:)
      
    open(newunit=fh,file=NAMERESTART,access='stream')
      
    read(fh) iter
    read(fh) ctime
      
    ! now read cell no. and solutions...
    allocate(tmp(numv,N,N))
      
    do I=1,NCELLGLOB
      
        read(fh) ICG
        read(fh) tmp(1:numv,1:N,1:N)
      
        IC=ICG2IC(ICG)
      
        if(IC.ne.0) Q(1:numv,1:N,1:N,IC)=tmp(1:numv,1:N,1:N)
    end do
      
    close(fh)
      
    deallocate(tmp)
      
    END SUBROUTINE READ_ALL_DATA_BINARY