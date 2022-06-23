    subroutine readinput
    use setup2d
    implicit none
    include 'mpif.h'

    open(30,FILE='QUAD.INP')

    read(30,*)
    read(30,*)
    read(30,*) NAMECEL
    read(30,*) NAMEVRT
    read(30,*) NAMEVRT8
    read(30,*) NAMEBND
    read(30,*) METIS_CELL
    read(30,*) NAMERESTART
    read(30,*)
    read(30,*) N
    read(30,*) CURVE_WALL
    read(30,*) vismode
    read(30,*) mu
    read(30,*) eta0
    read(30,*)
    read(30,*) k_stage
    read(30,*) dt
    read(30,*) MAXITER
    read(30,*) 
    read(30,*) DXCYCL
    read(30,*) DYCYCL
    read(30,*)
    read(30,*) ARTIV
    read(30,*)
    read(30,*) NRE
    read(30,*) nwrite
    read(30,*)
    read(30,*) restart
    read(30,*) restart_ord
    read(30,*)
    read(30,*) rinf
    read(30,*) pinf
    read(30,*) uinf
    read(30,*) vinf
    read(30,*) Bxinf
    read(30,*) Byinf
    
    close(30)


    end subroutine readinput