    SUBROUTINE INIT_SETUP
    use setup2d
    implicit none
    include 'mpif.h'
    
    integer :: ic,is,js,iv
    double precision :: Xis,Xjs,xx(2)
    double precision,allocatable :: xxs(:,:)

    numv = 9

    allocate(Xs(N),Xf(N+1))
    allocate(Lmat(N+1,N),Mmat(N+1,N))

    CALL setsandpoints(N,Xs,Xf)
    CALL computeLandM
    
    allocate(XXsolu(2,N,N,NCELL))
    allocate(XXfluxi(2,N+1,N,NCELL))
    allocate(XXfluxj(2,N,N+1,NCELL))
    
    allocate(Jac(N,N,NCELL))
    allocate(S1(2,2,N+1,N,NCELL))
    allocate(S2(2,2,N,N+1,NCELL))

    allocate(iface2fp(N,4))
    allocate(jface2fp(N,4))

    allocate(Q(numv,N,N,NCELL))
    allocate(Qvfi(numv,N+1,N,NCELL))
    allocate(Qvfj(numv,N,N+1,NCELL))
    allocate(F1(numv,N+1,N,NCELL))
    allocate(G2(numv,N,N+1,NCELL))
    allocate(resid(numv,N,N,NCELL))

    if (vismode==1) then
        allocate(Fv1(numv,N+1,N,NCELL))
        allocate(Gv2(numv,N,N+1,NCELL))
        if (ARTIV==1) then
            allocate(SM_IND(2,NCELL))
            allocate(muAV(2,NCELL))
            allocate(muAVfi(2,2,NCELL),muAVfj(2,2,NCELL))
            allocate(Xs_low(N-1),Xf_low(N))
            call setsandpoints(N-1,Xs_low,Xf_low)
        end if
    end if

    allocate(nablaQs(numv,2,N,N,NCELL))
    allocate(nablaQvfi(numv,2,N+1,N,NCELL))
    allocate(nablaQvfj(numv,2,N,N+1,NCELL))
    
    ! for processor interfaces
    allocate(Qfl_p(numv,N,NPROCINT))
    allocate(Qfr_p(numv,N,NPROCINT))
    ! for processor interfaces VISCOUS FLUX
    allocate(Qfl_p2(numv,N,NPROCINT))
    allocate(Qfr_p2(numv,N,NPROCINT))
    
    ! for remote cyclic faces
    allocate(Qfl_c(numv,N,NCYCREM))
    allocate(Qfr_c(numv,N,NCYCREM))
    ! for remote cyclic faces VISCOUS FLUX
    allocate(Qfl_c2(numv,N,NCYCREM))
    allocate(Qfr_c2(numv,N,NCYCREM))

    if (CURVE_WALL==1) then
        allocate(dmdxs(8,2,N,N))
        allocate(dmdxf1(8,2,N+1,N))
        allocate(dmdxf2(8,2,N,N+1))
        call MAPDER_8_NODE
        allocate(xxs(2,8))
    else if (CURVE_WALL==0) then
        allocate(dmdxs(4,2,N,N))
        allocate(dmdxf1(4,2,N+1,N))
        allocate(dmdxf2(4,2,N,N+1))
        call MAPDER
        allocate(xxs(2,4))
    end if

    call calcjacob

    do ic = 1,NCELL
        if (CURVE_WALL==0) then
            do iv = 1,4
            xxs(1,iv) = XV(IVCELL(ic,iv))
            xxs(2,iv) = YV(IVCELL(ic,iv))
            end do
        else if (CURVE_WALL==1) then
        end if
        
        do js = 1,N
        do is = 1,N
  
            Xis = Xs(is)
            Xjs = Xs(js)
  
            if (CURVE_WALL==0) then
                call xyCoor_atGps(4,xxs,Xis,Xjs,xx)
            else if (CURVE_WALL==1) then
                call xyCoor_atGps(8,xxs,Xis,Xjs,xx)
            end if
  
            XXsolu(1,is,js,ic) = xx(1)
            XXsolu(2,is,js,ic) = xx(2)
        end do
        end do
  
        ! i-direction
        do js = 1,N
        do is = 1,N+1
  
            Xis = Xf(is)
            Xjs = Xs(js)
  
            if (CURVE_WALL==0) then
                call xyCoor_atGps(4,xxs,Xis,Xjs,xx)
            else if (CURVE_WALL==1) then
                call xyCoor_atGps(8,xxs,Xis,Xjs,xx)
            end if
  
            XXfluxi(1,is,js,ic) = xx(1)
            XXfluxi(2,is,js,ic) = xx(2)
        end do
        end do
        
        ! j-direction
        do js = 1,N+1
        do is = 1,N
  
            Xis = Xs(is)
            Xjs = Xf(js)
  
            if (CURVE_WALL==0) then
                call xyCoor_atGps(4,xxs,Xis,Xjs,xx)
            else if (CURVE_WALL==1) then
                call xyCoor_atGps(8,xxs,Xis,Xjs,xx)
            end if
  
            XXfluxj(1,is,js,ic) = xx(1)
            XXfluxj(2,is,js,ic) = xx(2)
            
        end do
        end do

    end do

    call face2fpconnec

    if (NWALL.gt.0) then
        allocate(dedxi(N,N,NCELL))
        allocate(dedeta(N,N,NCELL))
    end if

    if (restart==0) then
        call SETINITIALCOND
    else if (restart==1) then
        if (restart_ord==N) then
            call READ_ALL_DATA_BINARY
        else
            if (rank.eq.0) print *,'restart_ord not equal to N, not available now'
        end if
    end if

    !! divergence cleaning
    alpha = 2.d0
    eigvmax = 0.d0
    umax = 0.d0

    END SUBROUTINE INIT_SETUP

    SUBROUTINE setsandpoints(np, Xs_tmp, Xf_tmp)
    use setup2d
    IMPLICIT NONE
    INTEGER :: i
    INTEGER, INTENT(IN) :: np
    DOUBLE PRECISION, INTENT(OUT) :: Xs_tmp(np), Xf_tmp(np + 1)

    DO i = 1, np
        Xs_tmp(i) = 0.5D0 * (1.0D0 - cos(dble(2 * i-1) * pi/(2 * dble(np))))
    END DO
    
    if (np == 1)then
        Xf_tmp(1) = 0.d0
        Xf_tmp(2) = 1.d0
    else if (np == 2)then
        Xf_tmp(1) = 0.d0
        Xf_tmp(2) = 0.5d0
        Xf_tmp(3) = 1.d0
    else if (np == 3)then
        Xf_tmp(1) = 0.d0
        Xf_tmp(2) = (1.d0 - dsqrt(3.d0)/3.d0)/2.d0
        Xf_tmp(3) = (1.d0 + dsqrt(3.d0)/3.d0)/2.d0
        Xf_tmp(4) = 1.d0
    else if (np == 4)then
        Xf_tmp(1) = 0.d0
        Xf_tmp(2) = (1.d0 - dsqrt(15.d0)/5.d0)/2.d0
        Xf_tmp(3) = 0.5d0
        Xf_tmp(4) = (1.d0 + dsqrt(15.d0)/5.d0)/2.d0
        Xf_tmp(5) = 1.d0
    else if (np == 5)then
        Xf_tmp(1) = 0.d0
        Xf_tmp(2) = (1.d0 - dsqrt(525.d0 + 70.d0 * dsqrt(30.d0))/35.d0)/2.d0
        Xf_tmp(3) = (1.d0 - dsqrt(525.d0-70.d0 * dsqrt(30.d0))/35.d0)/2.d0
        Xf_tmp(4) = (1.d0 + dsqrt(525.d0-70.d0 * dsqrt(30.d0))/35.d0)/2.d0
        Xf_tmp(5) = (1.d0 + dsqrt(525.d0 + 70.d0 * dsqrt(30.d0))/35.d0)/2.d0
        Xf_tmp(6) = 1.d0
    else if (np == 6)then
        Xf_tmp(1) = 0.d0
        Xf_tmp(2) = (1. - dsqrt(245.d0 + 14.d0 * dsqrt(70.d0))/21.d0)/2.d0
        Xf_tmp(3) = (1. - dsqrt(245.d0-14.d0 * dsqrt(70.d0))/21.d0)/2.d0
        Xf_tmp(4) = 0.5d0
        Xf_tmp(5) = (1. + dsqrt(245.d0-14.d0 * dsqrt(70.d0))/21.d0)/2.d0
        Xf_tmp(6) = (1. + dsqrt(245.d0 + 14.d0 * dsqrt(70.d0))/21.d0)/2.d0
        Xf_tmp(7) = 1.d0
    else if (np == 7)then
        Xf_tmp(1) = 0.d0
        Xf_tmp(2) = (1.d0 - 0.932469514203d0)/2.d0
        Xf_tmp(3) = (1.d0 - 0.661209386466d0)/2.d0
        Xf_tmp(4) = (1.d0 - 0.238619186083d0)/2.d0
        Xf_tmp(5) = (1.d0 + 0.238619186083d0)/2.d0
        Xf_tmp(6) = (1.d0 + 0.661209386466d0)/2.d0
        Xf_tmp(7) = (1.d0 + 0.932469514203d0)/2.d0
        Xf_tmp(8) = 1.d0
    else if (np == 8)then
        Xf_tmp(1) = 0.d0
        Xf_tmp(2) = (1.d0 - 0.949107912343d0)/2.d0
        Xf_tmp(3) = (1.d0 - 0.741531185599d0)/2.d0
        Xf_tmp(4) = (1.d0 - 0.405845151377d0)/2.d0
        Xf_tmp(5) = 0.5d0
        Xf_tmp(6) = (1.d0 + 0.405845151377d0)/2.d0
        Xf_tmp(7) = (1.d0 + 0.741531185599d0)/2.d0
        Xf_tmp(8) = (1.d0 + 0.949107912343d0)/2.d0
        Xf_tmp(9) = 1.d0
    else if (np == 9)then
        Xf_tmp(1) = 0.d0
        Xf_tmp(2) = (1.d0 - 0.960289856498d0)/2.d0
        Xf_tmp(3) = (1.d0 - 0.796666477414d0)/2.d0
        Xf_tmp(4) = (1.d0 - 0.525532409916d0)/2.d0
        Xf_tmp(5) = (1.d0 - 0.183434642496d0)/2.d0
        Xf_tmp(6) = (1.d0 + 0.183434642496d0)/2.d0
        Xf_tmp(7) = (1.d0 + 0.525532409916d0)/2.d0
        Xf_tmp(8) = (1.d0 + 0.796666477414d0)/2.d0
        Xf_tmp(9) = (1.d0 + 0.960289856498d0)/2.d0
        Xf_tmp(10) = 1.d0
    else if (np == 10)then
        Xf_tmp(1) = 0.d0
        Xf_tmp(2) = (1.d0 - 0.968160239508d0)/2.d0
        Xf_tmp(3) = (1.d0 - 0.836031107327d0)/2.d0
        Xf_tmp(4) = (1.d0 - 0.613371432701d0)/2.d0
        Xf_tmp(5) = (1.d0 - 0.324253423404d0)/2.d0
        Xf_tmp(6) = 0.5d0
        Xf_tmp(7) = (1.d0 + 0.324253423404d0)/2.d0
        Xf_tmp(8) = (1.d0 + 0.613371432701d0)/2.d0
        Xf_tmp(9) = (1.d0 + 0.836031107327d0)/2.d0
        Xf_tmp(10) = (1.d0 + 0.968160239508d0)/2.d0
        Xf_tmp(11) = 1.d0
    end if
    
    END SUBROUTINE setsandpoints
    
    
    SUBROUTINE computeLandM
    use setup2d
    IMPLICIT NONE

    integer :: np
    integer :: is, ifp, s, k, js
    double precision :: hhval
    double precision :: num, den, xval, sum
    
    np = N
    
    do ifp = 1, np + 1
        do is = 1, np
            Lmat(ifp, is) = hhval(is, Xf(ifp))
        end do
    end do
    
    do ifp = 1, np + 1
        do is = 1, np
       
            xval = Xs(is)
            sum = 0.0d0
            do k = 1, np + 1
                if (k/= ifp)then
                    num = 1.0d0
                    den = 1.0d0
                    do s = 1, np + 1
                    if (s/= ifp .and. s/= k)then
                        num = num * (xval - Xf(s))
                    end if
                    if (s/= ifp)then
                        den = den * (Xf(ifp) - Xf(s))
                    end if
       
                    end do
       
                    sum = sum + num/den
                end if
            end do ! end of loop over k
       
            Mmat(ifp, is) = sum
       
        end do
    end do 

    ! write(*,*) "Lmat"
    ! do ifp=1,np+1
    !   write(*,*) (Lmat(ifp,is), is=1,np)
    ! end do

    ! write(*,*) "Mmat"
    ! do ifp=1,np+1
    !   write(*,*) (Mmat(ifp,is), is=1,np)
    ! end do

    ! write (*,*) 'N=', np, 'Lmat and Mmat is generated'

    END SUBROUTINE computeLandM
    
    
    FUNCTION hhval(i, xval)
    use setup2d
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: i
    INTEGER :: np
    double precision, intent(in) :: xval
    double precision :: hhval, hvaln, hvald
    integer :: s
    
    hvaln = 1.0d0
    hvald = 1.0d0
    np = n
    do s = 1, np
        if (s/= i)then
            hvaln = hvaln * (xval - Xs(s))
            hvald = hvald * (Xs(i) - Xs(s))
        end if
    end do
    
    hhval = hvaln/hvald
    
    END FUNCTION hhval

    ! SUBROUTINE TO CONNECTS THE FLUX POINT INDICES TO THE LOCAL FACE INDEX
    SUBROUTINE face2fpconnec

    use setup2d
  
    IMPLICIT NONE
  
    integer	:: nfp,iface
  
    ! for face 1
    iface = 1
    do nfp=1,N
        iface2fp(nfp,iface) = nfp
        jface2fp(nfp,iface) = 1
    end do
  
    ! for face 2
    iface = 2
    do nfp=1,N
        iface2fp(nfp,iface) = N+1
        jface2fp(nfp,iface) = nfp
    end do
  
    ! for face 3
    iface = 3
    do nfp=1,N
        iface2fp(nfp,iface) = N - nfp + 1
        jface2fp(nfp,iface) = N+1
    end do
  
    ! for face 4
    iface = 4
    do nfp=1,N
        iface2fp(nfp,iface) = 1
        jface2fp(nfp,iface) = N - nfp + 1
    end do
  
    END SUBROUTINE face2fpconnec
  