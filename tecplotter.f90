    subroutine tecplotter2d
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: icell,ii,jj,Np,gnid,is,js,nd
    double precision :: xsi,eta,NORM_U,NORM_B
    double precision,allocatable :: Qs(:)
    
    Np = NRE
    allocate(Qs(numv))
    
    allocate(plotQ(9,plotnodes))
    allocate(plotAV(2,plotnodes))
    ! nine variables for plotting
    ! rho u v w pr Bx By Bz Az

    plotQ(:,:) = 0.d0
    plotAV(:,:) = 0.d0

    do icell = 1,NCELL
        do jj = 1,Np+1
        do ii = 1,Np+1
            xsi = dble(ii-1)/dble(Np)
            eta = dble(jj-1)/dble(Np)
            gnid = gnumnode(ii,jj,icell)
            Qs = 0.d0

            do js = 1,N
            do is = 1,N
                Qs(1:numv) = Qs(1:numv) + Q(1:numv,is,js,icell) * &
                hval(N,is,xsi) * hval(N,js,eta)
            end do
            end do
            plotQ(1:numv,gnid) = plotQ(1:numv,gnid) + Qs(1:numv)
            
            plotAV(1,gnid) = plotAV(1,gnid) + muAV(1,icell)
            plotAV(2,gnid) = plotAV(2,gnid) + muAV(2,icell)
            
        end do
        end do
    end do

    do nd = 1,plotnodes
        Qs(1:numv) = plotQ(1:numv,nd)/dble(gnode_factor(nd))
        plotQ(1,nd) = Qs(1)
        plotQ(2,nd) = Qs(2)/Qs(1)
        plotQ(3,nd) = Qs(3)/Qs(1)
        plotQ(4,nd) = Qs(4)/Qs(1)
        plotQ(6:9,nd) = Qs(6:9)
        NORM_U = plotQ(2,nd)**2+plotQ(3,nd)**2+plotQ(4,nd)**2
        NORM_B = plotQ(6,nd)**2+plotQ(7,nd)**2+plotQ(8,nd)**2
        plotQ(5,nd) = (Qs(5) - 0.5d0*Qs(1)*NORM_U - 0.5d0*NORM_B)*(gam-1)
        plotAV(1:2,nd) = plotAV(1:2,nd)/dble(gnode_factor(nd))
    end do

    contains

    !***************************************
	DOUBLE PRECISION FUNCTION hval(np,i,xval)
	
	use	setup2d
	IMPLICIT NONE
	
	integer, intent(in)	:: np, i
	double precision, intent(in)	:: xval
	double precision	:: hvaln, hvald
	integer	:: s
		
	hvaln = 1.d0
	hvald = 1.d0
	do s=1,np
		if( s/=i) then
			hvaln = hvaln * ( xval - Xs(s) ) 
			hvald = hvald * ( Xs(i) - Xs(s) )
		end if		
	end do
	
	hval = hvaln/hvald	
			
	END FUNCTION hval		


    end subroutine tecplotter2d

!================================================================
!          PARALLEL WRITING OF TECPLOT BINARY FILES
!                   Bin Zhang, Jul.2013
!================================================================
    SUBROUTINE tecpost_parallel
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: plotnodes_g, plotcells_g
    integer :: plotoffset_nodes, plotoffset_cells
    integer :: status(MPI_STATUS_SIZE)
    integer (kind=MPI_OFFSET_KIND) :: disp
    integer :: fh, INTSIZE, DBSIZE, REALSIZE

    integer, parameter :: nvar=13
    integer, parameter :: max_string_len=64
    integer :: itmp(max_string_len)
    character(max_string_len) :: filename, string
    character*6 :: decimal
    character*6  :: PATH = './TEC/'
    real(8) :: minmax(2,nvar), minmaxtmp(2,nvar)
    real(4),parameter :: ZONEMARKER=299.0, EOHMARKER=357.0
    integer :: pt

    disp = 0

    !gather each proc's number of plot nodes and cells
    call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    call MPI_ALLGATHER(plotnodes,1,MPI_INTEGER,PROC_plotnodes,1,MPI_INTEGER,MPI_COMM_WORLD,IERROR)
    call MPI_ALLGATHER(plotcells,1,MPI_INTEGER,PROC_plotcells,1,MPI_INTEGER,MPI_COMM_WORLD,IERROR)
    
    !get the total number of plot nodes and cells
    call MPI_ALLREDUCE(plotnodes,plotnodes_g,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
    call MPI_ALLREDUCE(plotcells,plotcells_g,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)

    !get the node and cell offset for each processor
    if (rank .eq. 0) then
        plotoffset_nodes = 0
        plotoffset_cells = 0
    else
        plotoffset_nodes = sum(PROC_plotnodes(1:rank))
        plotoffset_cells = sum(PROC_plotcells(1:rank))
    end if

    !update the vertice numbers
    connec_c2n(:,:) = connec_c2n(:,:) + plotoffset_nodes

    !get the integer size and double precision size
    call MPI_TYPE_SIZE(MPI_INTEGER,INTSIZE,IERROR)
    call MPI_TYPE_SIZE(MPI_REAL,REALSIZE,IERROR)
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,DBSIZE,IERROR)

    !get the min and max value of each varialbe
    !will be used in tecplot file header
    minmaxtmp(1, 1) = minval(plotX(1,:))
    minmaxtmp(1, 2) = minval(plotX(2,:))
    minmaxtmp(1, 3) = minval(plotQ(1,:))
    minmaxtmp(1, 4) = minval(plotQ(2,:))
    minmaxtmp(1, 5) = minval(plotQ(3,:))
    minmaxtmp(1, 6) = minval(plotQ(4,:))
    minmaxtmp(1, 7) = minval(plotQ(5,:))
    minmaxtmp(1, 8) = minval(plotQ(6,:))
    minmaxtmp(1, 9) = minval(plotQ(7,:))
    minmaxtmp(1,10) = minval(plotQ(8,:))
    minmaxtmp(1,11) = minval(plotQ(9,:))
    minmaxtmp(1,12) = minval(plotAV(1,:))
    minmaxtmp(1,13) = minval(plotAV(2,:))
    
    minmaxtmp(2, 1) = maxval(plotX(1,:))
    minmaxtmp(2, 2) = maxval(plotX(2,:))
    minmaxtmp(2, 3) = maxval(plotQ(1,:))
    minmaxtmp(2, 4) = maxval(plotQ(2,:))
    minmaxtmp(2, 5) = maxval(plotQ(3,:))
    minmaxtmp(2, 6) = maxval(plotQ(4,:))
    minmaxtmp(2, 7) = maxval(plotQ(5,:))
    minmaxtmp(2, 8) = maxval(plotQ(6,:))
    minmaxtmp(2, 9) = maxval(plotQ(7,:))
    minmaxtmp(2,10) = maxval(plotQ(8,:))
    minmaxtmp(2,11) = maxval(plotQ(9,:))
    minmaxtmp(2,12) = maxval(plotAV(1,:))
    minmaxtmp(2,13) = maxval(plotAV(2,:))

    call MPI_REDUCE(minmaxtmp(1,:),minmax(1,:),nvar,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,IERROR)
    call MPI_REDUCE(minmaxtmp(2,:),minmax(2,:),nvar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,IERROR)

    call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    !-----------------------------------------------------
    !now open a new tecplot binary file...
    !-----------------------------------------------------
    write(decimal,'(f6.3)') ctime-int(ctime)
    write(filename,'(a,i7.7,a4,a)')'TEC',int(ctime), decimal(3:6), '.plt'
    call MPI_FILE_OPEN(MPI_COMM_WORLD,PATH//filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,IERROR)

    !------------------------------------------------------
    !let the root proc write the FILE HEADER...
    !------------------------------------------------------
    if (rank .eq. 0) then
        disp = 0  !displacement from file beginning

        string = '#!TDV111'  !Magic Number (version)
        call MPI_FILE_WRITE(fh,string,8,MPI_CHARACTER,status,IERROR)
        disp = disp + 8

        itmp(1) = 1     !integer value of 1
        itmp(2) = 0     !file type (0=full,1=grid,2=solution)
        call MPI_FILE_WRITE(fh,itmp,2,MPI_INTEGER,status,IERROR)
        disp = disp + INTSIZE*2

        string = 'step data'    !title of data set
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        call MPI_FILE_WRITE(fh,nvar,1,MPI_INTEGER,status,IERROR)
        disp = disp + INTSIZE   !number of variables

        string = 'X'   !write variables names
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'Y'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'RHO'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'U'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'V'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'W'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'P'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'Bx'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'By'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'Bz'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'Az'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'muAV'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        string = 'etaAV'
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        call MPI_FILE_WRITE(fh,ZONEMARKER,1,MPI_REAL,status,IERROR)
        disp = disp + REALSIZE   !the Zone Marker

        string = 'ZONE 001'   !Zone Name
        call write_tecplot_string(string)
        disp = disp + INTSIZE*(len_trim(string)+1)

        itmp(1) = -1    !Parent Zone
        itmp(2) =  0    !Strand ID
        call MPI_FILE_WRITE(fh,itmp,2,MPI_INTEGER,status,IERROR)
        disp = disp + INTSIZE*2

        call MPI_FILE_WRITE(fh,ctime,1,MPI_DOUBLE_PRECISION,status,IERROR)
        disp = disp + DBSIZE    !Solution Time

        itmp(1) = -1    !Zone Color
        itmp(2) =  3    !Zone Type (3=FEQUADRILATERAL)
        itmp(3) =  0    !Data Packing (0=block, 1=point)
        itmp(4) =  0    !Var Location (0=nodes, 1=specify)
        itmp(5) =  0    !Face neighbour connectivity supplied?
        itmp(6) =  0    !Any user define face neighbour?
        itmp(7) =  plotnodes_g    !number of points
        itmp(8) =  plotcells_g    !number of elements
        itmp(9) =  0    !ICellDim (set to 0)
        itmp(10)=  0    !JCellDim (set to 0)
        itmp(11)=  0    !KCellDim (set to 0)
        itmp(12)=  0    !Any auxiliary data?
        call MPI_FILE_WRITE(fh,itmp,12,MPI_INTEGER,status,IERROR)
        disp = disp + INTSIZE*12

        call MPI_FILE_WRITE(fh,EOHMARKER,1,MPI_REAL,status,IERROR)
        call MPI_FILE_WRITE(fh,ZONEMARKER,1,MPI_REAL,status,IERROR)
        disp = disp + REALSIZE*2

        itmp(1:nvar) = 1  !Variable format: 1=float,2=double...
        itmp(nvar+1) = 0  !Has passive variables?
        itmp(nvar+2) = 0  !Has variable sharing?
        itmp(nvar+3) =-1  !Zone to share connectivity(-1=no)
        call MPI_FILE_WRITE(fh,itmp,nvar+3,MPI_INTEGER,status,IERROR)
        disp = disp + INTSIZE*(nvar+3)

        call MPI_FILE_WRITE(fh,minmax,2*nvar,MPI_DOUBLE_PRECISION,status,IERROR)
        disp = disp + DBSIZE*2*nvar     !min and max of each variable

    end if

    call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    !let all procs know the header displacement
    call MPI_BCAST(disp,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)

    !---- write 'X'
    disp = disp + REALSIZE*plotoffset_nodes
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotX(1,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'Y'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotX(2,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'RHO'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotQ(1,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'U'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotQ(2,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'V'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotQ(3,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'W'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotQ(4,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'P'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotQ(5,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'Bx'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotQ(6,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'By'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotQ(7,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'Bz'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotQ(8,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'Az'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotQ(9,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'muAV'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotAV(1,:)),plotnodes,MPI_REAL,status,IERROR)

    !---- write 'etaAV'
    disp = disp + REALSIZE*plotnodes_g
    call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,MPI_REAL,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,real(plotAV(2,:)),plotnodes,MPI_REAL,status,IERROR)

    !------------------------------------------------------
    !start writing connectivity...
    !------------------------------------------------------
    disp = disp + REALSIZE*plotnodes_g - REALSIZE*plotoffset_nodes + 4*INTSIZE*plotoffset_cells
    connec_c2n(:,:)=connec_c2n(:,:)-1   !NOTE: node no. starts from 0 in binary file!!!!!!!
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER,MPI_INTEGER,'native',MPI_INFO_NULL,IERROR)
    call MPI_FILE_WRITE_ALL(fh,connec_c2n,4*plotcells,MPI_INTEGER,status,IERROR)


    !------------------------------------------------------
    !DONE. Close the file...
    !------------------------------------------------------
    call MPI_FILE_CLOSE(fh,IERROR)

    DEALLOCATE(plotQ)
    deallocate(plotAV)

    deallocate(gnumcell,connec_c2n,plotX)
    deallocate(gnode_factor,gnumnode)
    deallocate(isfacenod,gfacenode)
	
	deallocate(PROC_plotnodes)
    deallocate(PROC_plotcells)

    contains
    !---------------------------------------------------
    !Convert string to tecplot format (integers plus 0)
    !---------------------------------------------------
    subroutine write_tecplot_string(string)

    implicit none

    character(len=*) :: string
    integer :: i, n

    n = len_trim(string)
    do i=1,n
        itmp(i) = ichar(string(i:i))
    end do
    itmp(n+1) = 0

    call MPI_FILE_WRITE(fh,itmp,n+1,MPI_INTEGER,status,IERROR)

    end subroutine write_tecplot_string

    END SUBROUTINE tecpost_parallel