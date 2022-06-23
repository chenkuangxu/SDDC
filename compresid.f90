    subroutine compresid
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: is,js,ic,rfp
    double precision :: dFpsi(9),dGeta(9),dFvpsi(9),dGveta(9),Qs(9),divB,&
    Bx,By,Bz,ux,uy,uz,psi

    call compflux

    do ic=1,NCELL

        ! need to find the flux derivatives at each solution pt
   
        do js=1,N
        do is=1,N

            Qs(1:9) = Q(1:9,is,js,ic)
   
            ! compute the dFdpsi derivative
            dFpsi(1:9) = 0.d0
            dFvpsi(1:9) = 0.d0
            do rfp=1,N+1
                dFpsi(1:9) = dFpsi(1:9) + F1(1:9,rfp,js,ic)*Mmat(rfp,is)
                if (vismode==1) then
                dFvpsi(1:9) = dFvpsi(1:9) + Fv1(1:9,rfp,js,ic)*Mmat(rfp,is)
                end if
            end do
   
            ! compute the dGdeta derivative
            dGeta(1:9) = 0.d0
            dGveta(1:9) = 0.d0
            do rfp=1,N+1
                dGeta(1:9) = dGeta(1:9) + G2(1:9,is,rfp,ic)*Mmat(rfp,js)
                if (vismode==1) then
                dGveta(1:9) = dGveta(1:9) + Gv2(1:9,is,rfp,ic)*Mmat(rfp,js)
                end if
            end do
   
            resid(1:9,is,js,ic) = -( dFpsi(1:9) + dGeta(1:9) )
   
            if (vismode==1) then
                resid(1:9,is,js,ic) = resid(1:9,is,js,ic) + &
                (dFvpsi(1:9)+dGveta(1:9))
            end if
   
        end do
        end do
   
    end do
    
    do ic=1,NCELL
        do js=1,N
        do is=1,N
            resid(1:9,is,js,ic)=resid(1:9,is,js,ic)/Jac(is,js,ic)
            !! source term
            divB = nablaQs(6,1,is,js,ic) + nablaQs(7,2,is,js,ic)
            Bx = Q(6,is,js,ic)
            By = Q(7,is,js,ic)
            Bz = Q(8,is,js,ic)
            ux = Q(2,is,js,ic) / Q(1,is,js,ic)
            uy = Q(3,is,js,ic) / Q(1,is,js,ic)
            uz = Q(4,is,js,ic) / Q(1,is,js,ic)
            psi= Q(9,is,js,ic)

            resid(2,is,js,ic) = resid(2,is,js,ic) - divB * Bx
            resid(3,is,js,ic) = resid(3,is,js,ic) - divB * By
            resid(4,is,js,ic) = resid(4,is,js,ic) - divB * Bz
            resid(5,is,js,ic) = resid(5,is,js,ic) - divB * (ux*Bx + uy*By + uz*Bz)
            resid(6,is,js,ic) = resid(6,is,js,ic) - divB * ux
            resid(7,is,js,ic) = resid(7,is,js,ic) - divB * uy
            resid(8,is,js,ic) = resid(8,is,js,ic) - divB * uz
            resid(9,is,js,ic) = resid(9,is,js,ic) - alpha * psi
            !!
        end do
        end do
    end do

    end subroutine compresid
    