    subroutine getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)
    use setup2d,only: gam,eigvmax,c_h,umax
    implicit none
    
    integer,intent(in) :: sign_l,sign_r
    double precision :: Qfl(9),Qfr(9),Fnl(9),Fnr(9),eigv,normf(2)
    double precision :: ul,vl,wl,ur,vr,wr,Vn_l,Vn_r,pl,pr,& 
    Bxl,Byl,Bzl,Bxr,Byr,Bzr,NORM_Bl,NORM_Br,Bx,By,Bz,Bn_l,Bn_r,psil,psir
    double precision :: rhoi,c_a,c_d,c_n,c_tmp,c_f
    double precision :: unit_x,unit_y,magnorm,umaxtmp

    ! get the magnitude of normal vector |J| div(\xi) or |J| div(\eta) from left cell
    magnorm = dsqrt(normf(1)**2+normf(2)**2)

    ! unit normal vector pointing from left cell to right cell
    unit_x = normf(1)/magnorm*sign_l
    unit_y = normf(2)/magnorm*sign_l 
    
    ! compute the Rosanov flux on physical domain along vector (unit_x,unit_y)
    
    !for the left side
    rhoi  = 1.d0/Qfl(1)
    ul    = rhoi*Qfl(2)
    vl    = rhoi*Qfl(3)
    wl    = rhoi*Qfl(4)
    Bxl   = Qfl(6)
    Byl   = Qfl(7)
    Bzl   = Qfl(8)
    psil  = Qfl(9)
    
    NORM_Bl = Bxl**2+Byl**2+Bzl**2
    
    pl    = (gam-1.d0)*(Qfl(5) - 0.5d0*Qfl(1)*(ul**2+vl**2+wl**2) &
    - NORM_Bl*0.5d0 - 0.5d0*psil**2)
    pl    = pl + 0.5d0*NORM_Bl  ! Total pressure = pressure + magnetic pressure
    Vn_l  = ul*unit_x+vl*unit_y
    Bn_l  = Bxl*unit_x + Byl*unit_y

    !for the right side
    rhoi  = 1.d0/Qfr(1)
    ur    = rhoi*Qfr(2)
    vr    = rhoi*Qfr(3)
    wr    = rhoi*Qfr(4)
    Bxr   = Qfr(6)
    Byr   = Qfr(7)
    Bzr   = Qfr(8)
    psir  = Qfr(9)
    
    NORM_Br = Bxr**2 + Byr**2 + Bzr**2
    
    pr    = (gam-1.d0)*(Qfr(5) - 0.5d0*Qfr(1)*(ur**2+vr**2+wr**2) &
    - NORM_Br*0.5d0 - 0.5d0*psir**2)
    pr    = pr + NORM_Br*0.5d0  ! Total pressure = pressure + magnetic pressure
    Vn_r  = ur*unit_x+vr*unit_y
    Bn_r  = Bxr*unit_x + Byr*unit_y
    
    ! calculate the normal flux component at left side
    Fnl(1) = Qfl(1)*Vn_l
    Fnl(2) = Qfl(2)*Vn_l + pl*unit_x - Bxl*Bn_l
    Fnl(3) = Qfl(3)*Vn_l + pl*unit_y - Byl*Bn_l
    Fnl(4) = Qfl(4)*Vn_l -Bzl*Bn_l
    Fnl(5) = (Qfl(5)-0.5d0*psil**2+pl)*Vn_l - Bn_l*(ul*Bxl+vl*Byl+wl*Bzl)
    Fnl(6) = (vl*Bxl-ul*Byl)*unit_y
    Fnl(7) = (ul*Byl-vl*Bxl)*unit_x
    Fnl(8) = Bzl*Vn_l - wl*Bn_l  
     
    ! calculate the normal flux component at right side
    Fnr(1) = Qfr(1)*Vn_r
    Fnr(2) = Qfr(2)*Vn_r + pr*unit_x - Bxr*Bn_r
    Fnr(3) = Qfr(3)*Vn_r + pr*unit_y - Byr*Bn_r
    Fnr(4) = Qfr(4)*Vn_r - Bzr*Bn_r
    Fnr(5) = (Qfr(5)-0.5d0*psir**2+pr)*Vn_r - Bn_r*(ur*Bxr+vr*Byr+wr*Bzr)
    Fnr(6) = (vr*Bxr-ur*Byr)*unit_y
    Fnr(7) = (ur*Byr-vr*Bxr)*unit_x
    Fnr(8) = Bzr*Vn_r - wr*Bn_r  

    !!!!!!!!! add divergence cleaning part
    Fnl(5) = Fnl(5) + c_h * psil * Bn_l
    Fnr(5) = Fnr(5) + c_h * psir * Bn_r

    Fnl(6) = Fnl(6) + c_h * psil * unit_x
    Fnr(6) = Fnr(6) + c_h * psir * unit_x

    Fnl(7) = Fnl(7) + c_h * psil * unit_y
    Fnr(7) = Fnr(7) + c_h * psir * unit_y

    Fnl(9) = c_h * Bn_l
    Fnr(9) = c_h * Bn_r
    !!!!!!!!!
    
    ! calculate the eignvalue (fast magneto-acoustic wave speed)
    pl = pl - 0.5d0*NORM_Bl
    pr = pr - 0.5d0*NORM_Br
    Bx = 0.5d0*(Bxl+Bxr)
    By = 0.5d0*(Byl+Byr)
    Bz = 0.5d0*(Bzl+Bzr)
    
    c_a  = gam*(pl+pr)/(Qfl(1)+Qfr(1))
    c_d  = (Bx**2+By**2+Bz**2)/(0.5d0*(Qfl(1)+Qfr(1)))    
    c_n  = (Bx*unit_x+By*unit_y)**2/(0.5d0*(Qfl(1)+Qfr(1)))
    c_tmp= dsqrt((c_a+c_d)**2-4.0d0*c_a*c_n)
    
    c_f  = dsqrt(0.5d0*(c_a+c_d+c_tmp))
    
    eigv = 0.5d0*dabs(Vn_l+Vn_r) + c_f
        
    ! Rusanov flux for left cell on physical domain
    Fnl(1:9) = 0.5d0*(Fnl(1:9)+Fnr(1:9)-eigv*(Qfr(1:9)-Qfl(1:9)))

    ! Transfer to computational domain for left cell
    Fnl(1:9) = magnorm*sign_l*Fnl(1:9)

    ! Transfer to computational domain for right cell
    Fnr(1:9) = sign_l*sign_r*Fnl(1:9)

    if (eigv.gt.eigvmax) then
        eigvmax = eigv
    end if

    umaxtmp = max(0.5d0*dabs(ul+ur),0.5d0*dabs(vl+vr),0.5d0*dabs(wl+wr))
    if (umaxtmp.gt.umax) then
        umax = umaxtmp
    end if

    end subroutine getrusanovflux

    subroutine getIfluxvectors(Qi,Fi,Gi)
    use setup2d,only: gam,c_h,umax
    implicit none

    double precision,intent(in) :: Qi(9)
    double precision,intent(out) :: Fi(9),Gi(9)
    double precision :: rho,rhou,rhov,rhow,rhoE,pr,Bx,By,Bz,NORM_B,psi,u,v,w,umaxtmp

    rho  = Qi(1)
    rhou = Qi(2)
    rhov = Qi(3)
    rhow = Qi(4)
    rhoE = Qi(5)
    Bx   = Qi(6)
    By   = Qi(7)
    Bz   = Qi(8)
    psi  = Qi(9)

    NORM_B = Bx**2+By**2+Bz**2

    !compute the pressure
    pr = (gam-1.d0)*(rhoE - 0.5d0*(rhou**2+rhov**2+rhow**2)/rho&
    -0.5d0*NORM_B-0.5d0*psi**2)

    !compute flux vector F
    Fi(1) = rhou
    Fi(2) = (rhou**2)/rho + pr + (0.5d0*NORM_B-Bx**2)
    Fi(3) = rhou*rhov/rho - Bx*By
    Fi(4) = rhou*rhow/rho - Bx*Bz
    Fi(5) = (rhou/rho)*(rhoE-0.5d0*psi**2+pr+0.5d0*NORM_B) &
    - (rhou/rho*Bx+rhov/rho*By+rhow/rho*Bz)*Bx
    Fi(6) = 0.d0
    Fi(7) = rhou/rho*By-rhov/rho*Bx
    Fi(8) = rhou/rho*Bz-rhow/rho*Bx 

    !compute flux vector G
    Gi(1) = rhov
    Gi(2) = rhou*rhov / rho - Bx*By
    Gi(3) = (rhov**2)/rho + pr + (0.5d0*NORM_B-By**2)
    Gi(4) = rhow*rhov/rho - By*Bz
    Gi(5) = (rhov/rho)*(rhoE-0.5d0*psi**2+pr+0.5d0*NORM_B) &
    - (rhou/rho*Bx+rhov/rho*By+rhow/rho*Bz)*By
    Gi(6) = rhov/rho*Bx - rhou/rho*By 
    Gi(7) = 0.d0
    Gi(8) = rhov/rho*Bz - rhow/rho*By  

    !!!!!!!!! add divergence cleaning part
    Fi(5) = Fi(5) + c_h * psi * Bx
    Gi(5) = Gi(5) + c_h * psi * By

    Fi(6) = Fi(6) + c_h * psi
    Gi(7) = Gi(7) + c_h * psi

    Fi(9) = c_h * Bx
    Gi(9) = c_h * By
    !!!!!!!!!

    u = rhou/rho
    v = rhov/rho
    w = rhow/rho
  
    umaxtmp = max(dabs(u),dabs(v),dabs(w))
    if (umaxtmp.gt.umax) then
        umax = umaxtmp
    end if

    end subroutine getIfluxvectors 

    !subroutine to calculate viscous flux
    !Qv and nablaQ are all smoothed ones for the interface flux points
    SUBROUTINE getVfluxvectors(Qv,nablaQ,Fv,Gv)

    use setup2d
    IMPLICIT NONE

    double precision,intent(in) :: Qv(numv),nablaQ(numv,2)
    double precision,intent(out) :: Fv(numv),Gv(numv)
    double precision :: u, v, w, NORM_U, Bx, By, Bz, Jx, Jy, Jz, NORM_B, inte
    double precision :: rho, rhou, rhov, rhow, rhoE, psi
    double precision :: ux, vx, uy, vy, wx, wy, ex, ey, rhox, rhoy, Bxy, Byx, Bzy, Bzx
    double precision :: du(2),dv(2),dw(2),drho(2),dBx(2),dBy(2),dBz(2),&
    de(2),dKE(2),dpre_B(2),dpsi(2),dpsi2(2)
  
    rho  = Qv(1)
    rhou = Qv(2)
    rhov = Qv(3)
    rhow = Qv(4)
    rhoE = Qv(5)
    Bx   = Qv(6)
    By   = Qv(7)
    Bz   = Qv(8)
    psi  = Qv(9)
  
    u = rhou/rho
    v = rhov/rho
    w = rhow/rho
  
    NORM_U = u**2  + v**2  + w**2
    NORM_B = Bx**2 + By**2 + Bz**2
  
    inte = (rhoE - 0.5d0*NORM_U*rho - 0.5d0*NORM_B - 0.5d0*psi**2)/rho
  
    drho(:) =  nablaQ(1,:)
    du(:)   = (nablaQ(2,:)-drho(:)*u)/rho
    dv(:)   = (nablaQ(3,:)-drho(:)*v)/rho
    dw(:)   = (nablaQ(4,:)-drho(:)*w)/rho
    dBx(:)  =  nablaQ(6,:)
    dBy(:)  =  nablaQ(7,:)
    dBz(:)  =  nablaQ(8,:)
    dpsi(:) =  nablaQ(9,:)
  
    dKE(:)    = 0.5d0*NORM_U*drho(:) + rho*(u*du(:)+v*dv(:)+w*dw(:))
    dpre_B(:) = Bx*dBx(:) + By*dBy(:) + Bz*dBz(:)

    dpsi2(:) = psi*dpsi(:)
  
    de(:)  = (nablaQ(5,:) - dKE(:) - dpre_B(:) - drho(:)*inte - dpsi2(:))/rho
  
    rhox = drho(1)
    rhoy = drho(2)
  
    ux = du(1)
    uy = du(2)
    vx = dv(1)
    vy = dv(2)
    wx = dw(1)
    wy = dw(2)
    ex = de(1)
    ey = de(2)
  
    Bxy = dBx(2)
    Byx = dBy(1)
    Bzx = dBz(1)
    Bzy = dBz(2)
  
    ! Viscous term from hydrodynamics
    Fv(1) = 0.d0
    Fv(2) = mu*(2.d0*ux+lambda*(ux+vy))
    Fv(3) = mu*(vx+uy)
    Fv(4) = mu*wx
    Fv(5) = u*Fv(2)+v*Fv(3)+w*Fv(4)+mu*gam*ex/prandt
  
    Gv(1) = 0.d0
    Gv(2) = mu*(vx+uy)
    Gv(3) = mu*(2.d0*vy+lambda*(ux+vy))
    Gv(4) = mu*wy
    Gv(5) = Gv(2)*u+Gv(3)*v+Gv(4)*w+mu*gam*ey/prandt
  
    ! Resistive term from magnetic
  
    Jx = Bzy   - 0.0d0
    Jy = 0.0d0 - Bzx
    Jz = Byx   - Bxy
  
    Fv(6) = 0.0d0
    Fv(7) = eta0*Jz
    Fv(8) = -eta0*Jy
    Fv(5) = Fv(5) + Fv(7)*By + Fv(8)*Bz
    
    Gv(6) = -eta0*Jz
    Gv(7) = 0.0d0
    Gv(8) = eta0*Jx
    Gv(5) = Gv(5) + Gv(6)*Bx + Gv(8)*Bz
  
    Fv(9) = 0.d0
    Gv(9) = 0.d0
  
    END SUBROUTINE getVfluxvectors

    SUBROUTINE getVfluxvectors_AV(Qv,nablaQ,Fv,Gv,fam,ifp,jfp,ic)
    use setup2d, only: numv,muAVfi,muAVfj,Xs,Xf,prandt,gam,lambda

    double precision :: Qv(numv),nablaQ(numv,2),Fv(numv),Gv(numv)
    integer,intent(in) :: ifp,jfp,ic,fam
    double precision :: epsilon_loc(numv),xsi,eta,mu_art,eta_art
    double precision :: u, v, w, NORM_U, Bx, By, Bz, NORM_B, inte
    double precision :: rho, rhou, rhov, rhow, rhoE, psi
    double precision :: ux, vx, uy, vy, wx, wy, ex, ey, rhox, rhoy, Bxy, Byx, Bzy, Bzx
    double precision :: Jx, Jy, Jz
    double precision :: du(2),dv(2),dw(2),drho(2),dBx(2),dBy(2),dBz(2),&
    de(2),dKE(2),dpre_B(2),dpsi(2),dpsi2(2)

    if (fam==1) then
        xsi = Xf(ifp)
        eta = Xs(jfp)
    else 
        xsi = Xs(ifp)
        eta = Xf(jfp)
    end if

    if (fam==1) then
        epsilon_loc(1:2) = muAVfi(1:2,1,ic) * (1.d0-xsi) + &
        muAVfi(1:2,2,ic) * xsi
    else
        epsilon_loc(1:2) = muAVfj(1:2,1,ic) * (1.d0-eta) + &
        muAVfj(1:2,2,ic) * eta
    end if

    ! mu_art = epsilon_loc(1)
    mu_art = 0.d0
    ! eta_art = epsilon_loc(1)*0.4d0
    eta_art = 0.d0

    rho  = Qv(1)
    rhou = Qv(2)
    rhov = Qv(3)
    rhow = Qv(4)
    rhoE = Qv(5)
    Bx   = Qv(6)
    By   = Qv(7)
    Bz   = Qv(8)
    psi  = Qv(9)
  
    u = rhou/rho
    v = rhov/rho
    w = rhow/rho
  
    NORM_U = u**2  + v**2  + w**2
    NORM_B = Bx**2 + By**2 + Bz**2
  
    inte = (rhoE - 0.5d0*NORM_U*rho - 0.5d0*NORM_B - 0.5d0*psi**2)/rho
  
    drho(:) =  nablaQ(1,:)
    du(:)   = (nablaQ(2,:)-drho(:)*u)/rho
    dv(:)   = (nablaQ(3,:)-drho(:)*v)/rho
    dw(:)   = (nablaQ(4,:)-drho(:)*w)/rho
    dBx(:)  =  nablaQ(6,:)
    dBy(:)  =  nablaQ(7,:)
    dBz(:)  =  nablaQ(8,:)
    dpsi(:) =  nablaQ(9,:)
  
    dKE(:)    = 0.5d0*NORM_U*drho(:) + rho*(u*du(:)+v*dv(:)+w*dw(:))
    dpre_B(:) = Bx*dBx(:) + By*dBy(:) + Bz*dBz(:)
  
    dpsi2(:) = psi*dpsi(:)
  
    de(:)  = (nablaQ(5,:) - dKE(:) - dpre_B(:) - drho(:)*inte - dpsi2(:))/rho
  
    rhox = drho(1)
    rhoy = drho(2)
  
    ux = du(1)
    uy = du(2)
    vx = dv(1)
    vy = dv(2)
    wx = dw(1)
    wy = dw(2)
    ex = de(1)
    ey = de(2)
  
    Bxy = dBx(2)
    Byx = dBy(1)
    Bzx = dBz(1)
    Bzy = dBz(2)
  
    ! Viscous term from hydrodynamics
    Fv(1) = 0.d0
    Fv(2) = mu_art*(2.d0*ux+lambda*(ux+vy))
    Fv(3) = mu_art*(vx+uy)
    Fv(4) = mu_art*wx
    Fv(5) = u*Fv(2)+v*Fv(3)+w*Fv(4)+mu_art*gam*ex/prandt
  
    Gv(1) = 0.d0
    Gv(2) = mu_art*(vx+uy)
    Gv(3) = mu_art*(2.d0*vy+lambda*(ux+vy))
    Gv(4) = mu_art*wy
    Gv(5) = Gv(2)*u+Gv(3)*v+Gv(4)*w+mu_art*gam*ey/prandt

    ! Resistive term from magnetic
    Jx = Bzy   - 0.0d0
    Jy = 0.0d0 - Bzx
    Jz = Byx   - Bxy
  
    Fv(6) = 0.0d0
    Fv(7) = eta_art*Jz
    Fv(8) = -eta_art*Jy
    Fv(5) = Fv(5) + Fv(7)*By + Fv(8)*Bz
    
    Gv(6) = -eta_art*Jz
    Gv(7) = 0.0d0
    Gv(8) = eta_art*Jx
    Gv(5) = Gv(5) + Gv(6)*Bx + Gv(8)*Bz

    Fv(numv) = 0.d0
    Gv(numv) = 0.d0

    END SUBROUTINE getVfluxvectors_AV