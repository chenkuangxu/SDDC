    SUBROUTINE SETINITIALCOND
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: ic,is,js
    double precision :: Qvi(numv),x,y,u,v,w,rho,Bx,By,Bz,&
    NORM_U,NORM_B,psi,pr,xx(2)

    do ic=1,NCELL
        do js=1,N
        do is=1,N
            xx(1:2) = XXsolu(1:2,is,js,ic)
            x = xx(1)
            y = xx(2)
            rho = rinf
            u = uinf + 0.02d0 * ran1()
            v = vinf + 0.02d0 * ran1()
            w = 0.d0
            pr = pinf
            Bx = Bxinf
            By = Byinf
            Bz = 0.d0
            psi= 0.d0
            NORM_U = u**2+v**2+w**2
            NORM_B = Bx**2+By**2+Bz**2
            Qvi(1) = rho
            Qvi(2) = rho*u
            Qvi(3) = rho*v
            Qvi(4) = rho*w
            Qvi(5) = pr/(gam-1) + 0.5d0*rho*NORM_U + 0.5d0*NORM_B&
            + 0.5d0*psi**2.d0
            Qvi(6) = Bx
            Qvi(7) = By
            Qvi(8) = Bz
            Qvi(9) = psi
            Q(1:numv,is,js,ic) = Qvi(1:numv)
        end do
        end do
    end do

    ctime = 1d-6
    iter = 0

    contains
	
	double precision function ran1()  !RETURNS RANDOM NUMBER BETWEEN 0 - 1
		
	double precision :: x
		
	!BUILT IN FORTRAN 90 RANDOM NUMBER FUNCTION
	call random_number(x)
		
	ran1 = x
	
	end function ran1

    END SUBROUTINE SETINITIALCOND