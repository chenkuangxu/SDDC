    subroutine tecplotter2dsetup
    use setup2d
    implicit none
    include 'mpif.h'

    double precision :: tolplot,xsi,eta,rx(2)
    integer :: icell,k,IC2V(4),Np,gnid,j,enode,iface,ii,jj,&
    existnode,checkmatch,kmax
    double precision,allocatable :: xxs(:,:)
    double precision,allocatable :: plotX_tmp(:,:)

    Np = NRE
    call allocate_plot_memory

    plotnodes = 0
    plotcells = 0
    tolplot = 1d-5
    gnode_factor = 0
    isfacenod = 0

    if (CURVE_WALL==0) then
        allocate(xxs(2,4))
        kmax = 4
    else if (CURVE_WALL==1) then
        allocate(xxs(2,8))
        kmax = 8
    end if

    ! First set up the 4 vertices of each cell
    ! Basically creating new global nodes
    do icell = 1,NCELL
        do k = 1,4
            IC2V(k) = IVCELL(icell,k)
        end do

        ! 4 vertices, lets set global node number for corr i,j
        gnumnode(1,1,icell) = IC2V(1)
        gnumnode(Np+1,1,icell) = IC2V(2)
        gnumnode(Np+1,Np+1,icell) = IC2V(3)
        gnumnode(1,Np+1,icell) = IC2V(4)

        do k = 1,4
            gnid = IC2V(k)
            if (gnode_factor(gnid)==0) then
                plotnodes = plotnodes+1
                plotX(1,gnid) = XV(gnid)
                plotX(2,gnid) = YV(gnid)
            end if
            gnode_factor(gnid) = gnode_factor(gnid) + 1
        end do

    end do ! end of loop over cells

    if (plotnodes/=NVERT) then
        write(*,*) 'something wrong in tecplotter2dsetup, my friend'
    end if

    do icell = 1,NCELL
        if (CURVE_WALL==0) then
            do j = 1,4
                xxs(1,j) = XV(IVCELL(icell,j))
                xxs(2,j) = YV(IVCELL(icell,j))
            end do
        else if (CURVE_WALL==1) then
            do j = 1,8
                xxs(1:2,j) = xxf_cell(1:2,j,icell)
            end do
        end if

        ! Let's loop over faces to check if 'nodalized'
        do k = 1,4
            iface = IC2F(icell,k)
            do enode = 1,Np-1
                select case(k)
                case(1)
                    ii = enode+1
                    jj = 1
                case(2)
                    ii = Np+1
                    jj = enode+1
                case(3)
                    ii = Np+1-enode
                    jj = Np+1
                case(4)
                    ii = 1
                    jj = Np+1-enode
                end select

                xsi = dble(ii-1)/dble(Np)
                eta = dble(jj-1)/dble(Np)

                call xyCoor_atGps(kmax,xxs,xsi,eta,rx)

                if (isfacenod(iface)==1) then
                    checkmatch = 0
                    do existnode=1,Np-1
                        gnid = gfacenode(existnode,iface)
                        if (distcalc(rx,gnid) <= tolplot) then
                            gnumnode(ii,jj,icell) = gnid
                            gnode_factor(gnid) = gnode_factor(gnid)+1
                            checkmatch = 1
                        end if
                    end do

                    if (checkmatch==0) then
                        write(*,*) 'something wrong in tecplotter2dsetup, face!'
                        print *,'rx',rx
                        do existnode=1,Np-1
                            gnid = gfacenode(existnode,iface)
                            print *,gnid,XV(gnid),YV(gnid)
                        end do
                    end if

                else
                    plotnodes = plotnodes+1
                    plotX(1,plotnodes) = rx(1)
                    plotX(2,plotnodes) = rx(2)
                    gnumnode(ii,jj,icell) = plotnodes
                    gfacenode(enode,iface) = plotnodes
                    gnode_factor(plotnodes) = gnode_factor(plotnodes)+1
                end if

            end do ! end of loop over nodes on face

            if (isfacenod(iface)==0) isfacenod(iface) = 1

        end do ! end of loop over nodes on face

        do ii = 2,Np
        do jj = 2,Np
            xsi = dble(ii-1)/dble(Np)
            eta = dble(jj-1)/dble(Np)
            call xyCoor_atGps(kmax,xxs,xsi,eta,rx)
            plotnodes = plotnodes + 1
            plotX(1,plotnodes) = rx(1)
            plotX(2,plotnodes) = rx(2)
            gnumnode(ii,jj,icell) = plotnodes
            gnode_factor(plotnodes) = gnode_factor(plotnodes)+1
        end do
        end do

        do ii = 1,Np
        do jj = 1,Np
            plotcells = plotcells + 1
            gnumcell(ii,jj,icell) = plotcells
            connec_c2n(1,plotcells) = gnumnode(ii,jj,icell)
            connec_c2n(2,plotcells) = gnumnode(ii+1,jj,icell)
            connec_c2n(3,plotcells) = gnumnode(ii+1,jj+1,icell)
            connec_c2n(4,plotcells) = gnumnode(ii,jj+1,icell)
        end do
        end do

    end do ! loop over the cells

    allocate(plotX_tmp(2,plotnodes))
    plotX_tmp(1:2,1:plotnodes) = plotX(1:2,1:plotnodes)
    deallocate(plotX)
    allocate(plotX(2,plotnodes))
    plotX(1:2,1:plotnodes) = plotX_tmp(1:2,1:plotnodes)

    ! write(*,*) 'rank=',rank,'plotnodes,plotcells:',plotnodes,plotcells
    
    contains

    double precision function distcalc(xx,gnid)
    use setup2d
    implicit none

    integer,intent(in) :: gnid
    double precision,intent(in) :: xx(2)

    distcalc = sqrt( (xx(1)-plotX(1,gnid))**2 + &
    (xx(2)-plotX(2,gnid))**2 )

    distcalc = abs(distcalc)

    end function distcalc


    end subroutine tecplotter2dsetup


    subroutine allocate_plot_memory
    use setup2d
    implicit none

    allocate(gnumcell(NRE,NRE,NCELL),connec_c2n(4,NRE**2*NCELL),plotX(2,(NRE+1)**2*NCELL))
    allocate(gnode_factor((NRE+1)**2*NCELL),gnumnode(NRE+1,NRE+1,NCELL))
    allocate(isfacenod(NFACE),gfacenode(NRE-1,NFACE))

    allocate(PROC_plotnodes(NPROC))
    allocate(PROC_plotcells(NPROC))

    end subroutine allocate_plot_memory