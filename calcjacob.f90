    subroutine calcjacob
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: is,js,k, ic,ifp,jfp
    integer, dimension(4) :: IC2V
    double precision :: xpsi,ypsi,xeta,yeta
    double precision,allocatable :: xxs(:,:)
    integer :: kmax

    if (CURVE_WALL==1) then
        kmax = 8
    else if (CURVE_WALL==0) then
        kmax = 4
    end if

    allocate(xxs(2,kmax))
    
    do ic=1,NCELL

        if (CURVE_WALL==0) then
            do k =1,kmax
                IC2V(k)=IVCELL(ic,k)
            end do
        else if (CURVE_WALL==1) then
            do k = 1,kmax
                xxs(1:2,k) = xxf_cell(1:2,k,ic)
            end do
        end if 
   
        ! loop to find the |J| at the solution points
        ! also compute the volume of cell
        do js=1,N
            do is=1,N
                xpsi = 0.d0
                xeta = 0.d0
                ypsi = 0.d0
                yeta = 0.d0
                if (CURVE_WALL==0) then
                    do k=1,kmax
                    xpsi = xpsi + dmdxs(k,1,is,js) * XV(IC2V(k))
                    ypsi = ypsi + dmdxs(k,1,is,js) * YV(IC2V(k))
                    xeta = xeta + dmdxs(k,2,is,js) * XV(IC2V(k))
                    yeta = yeta + dmdxs(k,2,is,js) * YV(IC2V(k))
                end do
                else if (CURVE_WALL==1) then
                    do k = 1,kmax
                    xpsi = xpsi + dmdxs(k,1,is,js) * xxs(1,k)
                    ypsi = ypsi + dmdxs(k,1,is,js) * xxs(2,k)
                    xeta = xeta + dmdxs(k,2,is,js) * xxs(1,k)
                    yeta = yeta + dmdxs(k,2,is,js) * xxs(2,k)
                    end do
                end if
                Jac(is,js,ic) = xpsi*yeta - xeta*ypsi
            end do
        end do
   
        ! loop to find the S vector at the flux points (family 1)
   
        do jfp=1,N
            do ifp=1,N+1
                xpsi = 0.d0
                xeta = 0.d0
                ypsi = 0.d0
                yeta = 0.d0

                if (CURVE_WALL==0) then
                    do k=1,kmax
                    xpsi = xpsi + dmdxf1(k,1,ifp,jfp) * XV(IC2V(k))
                    ypsi = ypsi + dmdxf1(k,1,ifp,jfp) * YV(IC2V(k))
                    xeta = xeta + dmdxf1(k,2,ifp,jfp) * XV(IC2V(k))
                    yeta = yeta + dmdxf1(k,2,ifp,jfp) * YV(IC2V(k))
                end do
                else if (CURVE_WALL==1) then
                    do k = 1,kmax
                    xpsi = xpsi + dmdxf1(k,1,ifp,jfp) * xxs(1,k)
                    ypsi = ypsi + dmdxf1(k,1,ifp,jfp) * xxs(2,k)
                    xeta = xeta + dmdxf1(k,2,ifp,jfp) * xxs(1,k)
                    yeta = yeta + dmdxf1(k,2,ifp,jfp) * xxs(2,k)
                    end do
                end if 
                
                ! S are the elements of the inverse jacobian
                S1(1,1,ifp,jfp,ic) = yeta
                S1(1,2,ifp,jfp,ic) = -xeta
   
                S1(2,1,ifp,jfp,ic) = -ypsi
                S1(2,2,ifp,jfp,ic) = xpsi
            end do
        end do
   
        ! loop to find the S vector at the flux points (family 2)
   
        do jfp=1,N+1
            do ifp=1,N
   
                xpsi = 0.d0
                xeta = 0.d0
                ypsi = 0.d0
                yeta = 0.d0
                
                if (CURVE_WALL==0) then
                    do k=1,kmax
                    xpsi = xpsi + dmdxf2(k,1,ifp,jfp) * XV(IC2V(k))
                    ypsi = ypsi + dmdxf2(k,1,ifp,jfp) * YV(IC2V(k))
                    xeta = xeta + dmdxf2(k,2,ifp,jfp) * XV(IC2V(k))
                    yeta = yeta + dmdxf2(k,2,ifp,jfp) * YV(IC2V(k))
                    end do
                else if (CURVE_WALL==1) then
                    do k = 1,kmax
                    xpsi = xpsi + dmdxf2(k,1,ifp,jfp) * xxs(1,k)
                    ypsi = ypsi + dmdxf2(k,1,ifp,jfp) * xxs(2,k)
                    xeta = xeta + dmdxf2(k,2,ifp,jfp) * xxs(1,k)
                    yeta = yeta + dmdxf2(k,2,ifp,jfp) * xxs(2,k)
                    end do
                end if
   
                ! S are the elements of the inverse jacobian
   
                S2(1,1,ifp,jfp,ic) = yeta
                S2(1,2,ifp,jfp,ic) = -xeta
   
                S2(2,1,ifp,jfp,ic) = -ypsi
                S2(2,2,ifp,jfp,ic) = xpsi
   
            end do
        end do
    end do
   

    end subroutine calcjacob