    SUBROUTINE mapder

    use setup2d

    IMPLICIT NONE

    integer :: is,js,ifp,jfp
    double precision :: psi, eta

    ! compute the derivatives at the solution points

    do js=1,N
    do is=1,N
        
        psi = Xs(is)
        eta = Xs(js)
        !
        dmdxs(1,1,is,js) = eta - 1.d0
        dmdxs(1,2,is,js) = psi - 1.d0
        !
        dmdxs(2,1,is,js) = 1.d0 - eta
        dmdxs(2,2,is,js) = -psi
        !
        dmdxs(3,1,is,js) = eta
        dmdxs(3,2,is,js) = psi
        !
        dmdxs(4,1,is,js) = -eta
        dmdxs(4,2,is,js) = 1.d0 - psi

    end do
    end do

    ! compute the derivatives at the flux points (family 1)
    do jfp=1,N
    do ifp=1,N+1
        !
        psi = Xf(ifp)
        eta = Xs(jfp)
        !
        dmdxf1(1,1,ifp,jfp) = eta - 1.d0
        dmdxf1(1,2,ifp,jfp) = psi - 1.d0
        !
        dmdxf1(2,1,ifp,jfp) = 1.d0 - eta
        dmdxf1(2,2,ifp,jfp) = -psi
        !
        dmdxf1(3,1,ifp,jfp) = eta
        dmdxf1(3,2,ifp,jfp) = psi
        !
        dmdxf1(4,1,ifp,jfp) = -eta
        dmdxf1(4,2,ifp,jfp) = 1.d0 - psi
        !
    end do
    end do
    !
    ! compute the derivatives at the flux points (family 2)

    do jfp=1,N+1
    do ifp=1,N
        psi = Xs(ifp)
        eta = Xf(jfp)
        dmdxf2(1,1,ifp,jfp) = eta - 1.d0
        dmdxf2(1,2,ifp,jfp) = psi - 1.d0
        !
        dmdxf2(2,1,ifp,jfp) = 1.d0 - eta
        dmdxf2(2,2,ifp,jfp) = -psi
        !
        dmdxf2(3,1,ifp,jfp) = eta
        dmdxf2(3,2,ifp,jfp) = psi
        !
        dmdxf2(4,1,ifp,jfp) = -eta
        dmdxf2(4,2,ifp,jfp) = 1.d0 - psi
        !
    end do
    end do

    END SUBROUTINE mapder

!****************************************************************************
    SUBROUTINE MapBaseFunc(func,xi,eta)
    IMPLICIT NONE
    DOUBLE PRECISION :: func(4),xi,eta
    func(1)= (1.0-xi)*(1.0-eta)
    func(2)= (1.0-eta)*xi
    func(3)= xi*eta
    func(4)= (1.0-xi)*eta
    RETURN
    END SUBROUTINE MapBaseFunc
!****************************************************************************

!****************************************************************************
    SUBROUTINE MapBaseFunc8(func8,xi,eta)
     IMPLICIT NONE
     DOUBLE PRECISION :: func8(8),xi,eta
     func8(1)= (1.d0-xi)*(1.d0-eta)*(-2.d0*xi-2.d0*eta+1.d0)
     func8(2)= xi*(1.d0-eta)*(2.d0*xi-2.d0*eta-1.d0)
     func8(3)= xi*eta*(2.d0*xi+2.d0*eta-3.d0)
     func8(4)= (1.d0-xi)*eta*(-2.d0*xi+2.d0*eta-1.d0)
     func8(5)= 4.d0*(1.d0-eta)*(1.d0-xi)*xi
     func8(6)= 4.d0*eta*(1.d0-eta)*xi
     func8(7)= 4.d0*(1.d0-xi)*xi*eta
     func8(8)= 4.d0*(1.d0-eta)*(1.d0-xi)*eta
     RETURN
     END SUBROUTINE MapBaseFunc8
 !***************************************************************************
 
 
 !!!!!!!!!!!!orientation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !#     4---------7--------3
 !      |                  |
 !      |                  |
 !      8                  6
 !      |                  |
 !      |                  |
 !      |                  |
 !      1---------5--------2
 !****************************************************************
     SUBROUTINE MAPDER_8_NODE
     USE setup2d
     IMPLICIT NONE
   
 
     INTEGER	:: IS,JS,IFP,JFP
     double precision :: psi, eta
     ! compute the derivatives at the solution points
     do js=1,N
     do is=1,N
         psi = Xs(is)
         eta = Xs(js)
         !
         dmdxs(1,1,is,js) = &
              (1.d0-eta)*(4.d0*psi+2.d0*eta-3.d0)
         dmdxs(1,2,is,js) = &
              (1.d0-psi) * (4.d0*eta+2.d0*psi-3.d0)
 
         dmdxs(2,1,is,js) = &
              (1.d0-eta) * (4.d0*psi - 2.d0*eta -1.d0)
         dmdxs(2,2,is,js) =  &
              psi * (4.d0*eta - 2.d0*psi -1.d0)
         !
         dmdxs(3,1,is,js) =  &
              eta * (4.d0*psi + 2.d0*eta -3.d0)
         dmdxs(3,2,is,js) =  &
              psi * (4.d0*eta + 2.d0*psi -3.d0)
         !
         dmdxs(4,1,is,js) = &
              eta * (4.d0*psi -2.d0*eta -1.d0)
         dmdxs(4,2,is,js) = &
              (1.d0-psi)*(4.d0*eta - 2.d0*psi -1.d0)
         !
         dmdxs(5,1,is,js) = 4.d0*(1.d0-eta)*(1.d0-2.d0*psi)
         dmdxs(5,2,is,js) = 4.d0*psi*(psi-1.d0)
         !
         dmdxs(6,1,is,js) =  4.d0*eta*(1.d0-eta)
         dmdxs(6,2,is,js) =  4.d0*psi*(1.d0-2.d0*eta)
         !
         dmdxs(7,1,is,js) =  4.d0*eta*(1.d0-2.d0*psi)
         dmdxs(7,2,is,js) =  4.d0*psi*(1.d0-psi)
         !
         dmdxs(8,1,is,js) =  4.d0*eta*(eta-1.d0)
         dmdxs(8,2,is,js) =  4.d0*(1.d0-psi)*(1.d0-2.d0*eta)
     end do
     end do
 
     !	% compute the derivatives at the flux points (family 1)
     do jfp=1,N
     do ifp=1,N+1
         !
         psi = Xf(ifp)
         eta = Xs(jfp)
         !
         dmdxf1(1,1,ifp,jfp) = (1.d0-eta)*(4.d0*psi+2.d0*eta-3.d0)
         dmdxf1(1,2,ifp,jfp) = (1.d0-psi) * (4.d0*eta+2.d0*psi-3.d0)
 
         dmdxf1(2,1,ifp,jfp) = &
              (1.d0-eta) * (4.d0*psi - 2.d0*eta -1.d0)
         dmdxf1(2,2,ifp,jfp) = &
              psi * (4.d0*eta - 2.d0*psi -1.d0)
         !
         dmdxf1(3,1,ifp,jfp) = &
              eta * (4.d0*psi + 2.d0*eta -3.d0)
         dmdxf1(3,2,ifp,jfp) = &
              psi * (4.d0*eta + 2.d0*psi -3.d0)
         !
         dmdxf1(4,1,ifp,jfp) = &
              eta * (4.d0*psi -2.d0*eta -1.d0)
         dmdxf1(4,2,ifp,jfp) = (1.d0-psi)*(4.d0*eta - 2.d0*psi -1.d0)
 
         dmdxf1(5,1,ifp,jfp) =  4.d0*(1.d0-eta)*(1.d0-2.d0*psi)
         dmdxf1(5,2,ifp,jfp) = 4.d0*psi*(psi-1.d0)
 
         dmdxf1(6,1,ifp,jfp) = &
              4.d0*eta*(1.d0-eta)
         dmdxf1(6,2,ifp,jfp) = &
              4.d0*psi*(1.d0-2.d0*eta)
 
         dmdxf1(7,1,ifp,jfp) = &
              4.d0*eta*(1.d0-2.d0*psi)
         dmdxf1(7,2,ifp,jfp) = &
              4.d0*psi*(1.d0-psi)
 
         dmdxf1(8,1,ifp,jfp) = &
              4.d0*eta*(eta-1.d0)
         dmdxf1(8,2,ifp,jfp) = &
              4.d0*(1.d0-psi)*(1.d0-2.d0*eta)
     !
     end do
     end do
     !
     ! compute the derivatives at the flux points (family 2)
 
     do jfp=1,N+1
     do ifp=1,N
         psi = Xs(ifp)
         eta = Xf(jfp)
         dmdxf2(1,1,ifp,jfp) = &
              (1.d0-eta)*(4.d0*psi+2.d0*eta-3.d0)
         dmdxf2(1,2,ifp,jfp) = &
              (1.d0-psi) * (4.d0*eta+2.d0*psi-3.d0)
         !
         dmdxf2(2,1,ifp,jfp) = &
              (1.d0-eta) * (4.d0*psi - 2.d0*eta -1.d0)
         dmdxf2(2,2,ifp,jfp) = &
              psi * (4.d0*eta - 2.d0*psi -1.d0)
         !
         dmdxf2(3,1,ifp,jfp) = &
              eta * (4.d0*psi + 2.d0*eta -3.d0)
         dmdxf2(3,2,ifp,jfp) = &
              psi * (4.d0*eta + 2.d0*psi -3.d0)
         !
         dmdxf2(4,1,ifp,jfp) = &
              eta * (4.d0*psi -2.d0*eta -1.d0)
         dmdxf2(4,2,ifp,jfp) =  (1.d0-psi)*(4.d0*eta - 2.d0*psi -1.d0)
 
         dmdxf2(5,1,ifp,jfp) =  4.d0*(1.d0-eta)*(1.d0-2.d0*psi)
         dmdxf2(5,2,ifp,jfp) =  4.d0*psi*(psi-1.d0)
 
         dmdxf2(6,1,ifp,jfp) =  4.d0*eta*(1.d0-eta)
         dmdxf2(6,2,ifp,jfp) =  4.d0*psi*(1.d0-2.d0*eta)
 
         dmdxf2(7,1,ifp,jfp) = &
              4.d0*eta*(1.d0-2.d0*psi)
         dmdxf2(7,2,ifp,jfp) = &
              4.d0*psi*(1.d0-psi)
 
         dmdxf2(8,1,ifp,jfp) = &
              4.d0*eta*(eta-1.d0)
         dmdxf2(8,2,ifp,jfp) = &
              4.d0*(1.d0-psi)*(1.d0-2.d0*eta)
         !
     end do
     end do
 
     END SUBROUTINE MAPDER_8_NODE

    SUBROUTINE xyCoor_atGps(ns,xxs,xi,eta,xx)
    IMPLICIT NONE
    integer,intent(in)::ns
    integer :: i,nb
    double precision :: xxs(2,ns),xx(2),xi,eta
    double precision,pointer,dimension(:)::func
    allocate(func(ns))
    if(ns==4) call MapBaseFunc(func,xi,eta)
    if(ns==8) call MapBaseFunc8(func,xi,eta)

    do i=1,2
     xx(i)= 0.d0
          do nb = 1,ns
          xx(i) = xx(i)+xxs(i,nb)*func(nb)
          end do
     end do
     deallocate(func)
     return
     END SUBROUTINE xyCoor_atGps