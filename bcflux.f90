    SUBROUTINE BCFLUX
    use setup2d
    implicit none

    integer :: ifn,iface,icleft,ifacelc,faml,ifpl,jfpl,nfp,ic,sp
    integer :: sign_l,sign_r
    double precision,dimension(9) :: Qfl,Qfr,Qs
    double precision,dimension(9) :: Fnl,Fnr
    double precision,dimension(2) :: normf
    double precision :: eigv
    double precision :: rhor,ur,vr,wr,V2,pr,Bxr,Byr,Bzr,&
    rhow,uw,vw,ww,magnorm,VN,BN

    sign_r = 1
    ! FIRSTLY THE INLET BC'S

    do ifn=1,NINLET
        
        iface = IBFINL(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1

        do nfp=1,N

            ! lets get the left solution vectors at flux point

            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft

            sign_l = 1

            Qfl(1:9) = 0.d0
            if(faml==1) then
                do sp=1,N
                    Qs(1:9) = Q(1:9,sp,jfpl,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(ifpl,sp)
                end do
                if(ifpl==1) sign_l=-1
                normf(1:2) = S1(1,1:2,ifpl,jfpl,ic)
            else
                do sp=1,N
                    Qs(1:9) = Q(1:9,ifpl,sp,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(jfpl,sp)
                end do
                if(jfpl==1) sign_l=-1
                normf(1:2) = S2(2,1:2,ifpl,jfpl,ic)
            end if

            rhor = rinf
            ur   = uinf
            vr   = vinf
            wr   = 0.d0
            V2   = ur*ur+vr*vr+wr*wr
            Bxr  = Bxinf
            Byr  = Byinf
            Bzr  = 0.d0
            
            pr = (gam-1.d0)*(Qfl(5) - 0.5d0*(Qfl(2)**2+Qfl(3)**2+Qfl(4)**2)&
            /Qfl(1) - 0.5d0*(Qfl(6)**2+Qfl(7)**2+Qfl(8)**2) - 0.5d0*Qfl(9)**2)

            Qfr(1) = rhor
            Qfr(2) = rhor*ur
            Qfr(3) = rhor*vr
            Qfr(4) = rhor*wr
            Qfr(5) = pinf/(gam-1)+0.5d0*rhor*(ur**2+vr**2+wr**2) &
            + 0.5d0*(Bxr**2+Byr**2+Bzr**2) + 0.5d0*Qfl(9)**2.d0
            Qfr(6) = Bxr
            Qfr(7) = Byr
            Qfr(8) = Bzr
            Qfr(9) = Qfl(9)

            CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)

            if(faml==1) then
                F1(1:9,ifpl,jfpl,icleft)  = Fnl(1:9)
            else
                G2(1:9,ifpl,jfpl,icleft)  = Fnl(1:9)
            end if

            if (vismode==1) then
                if(faml==1) then
                    Qvfi(1:9,ifpl,jfpl,icleft) = 0.5d0*Qfr(1:9)+0.5d0*Qfl(1:9)
                else
                    Qvfj(1:9,ifpl,jfpl,icleft) = 0.5d0*Qfr(1:9)+0.5d0*Qfl(1:9)
                end if
            end if 

        end do
    end do

    do ifn=1,NOUTLET
        iface = IBFOUT(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1
        do nfp=1,N
            ! lets get the left solution vectors at flux point
            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft

            sign_l = 1
            Qfl(1:9) = 0.d0
            if(faml==1) then
                do sp=1,N
                    Qs(1:9)  = Q(1:9,sp,jfpl,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(ifpl,sp)
                end do
                if(ifpl==1) sign_l=-1
                normf(1:2) = S1(1,1:2,ifpl,jfpl,ic)
            else
                do sp=1,N
                    Qs(1:9) = Q(1:9,ifpl,sp,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(jfpl,sp)
                end do
                if(jfpl==1) sign_l=-1
                normf(1:2) = S2(2,1:2,ifpl,jfpl,ic)
            end if

            rhor = Qfl(1)
            ur   = Qfl(2)/Qfl(1)
            vr   = Qfl(3)/Qfl(1)
            wr   = 0.d0
            V2   = ur*ur+vr*vr+wr*wr
            Bxr  = Bxinf
            Byr  = Byinf
            Bzr  = 0.d0
            
            pr = (gam-1.d0)*(Qfl(5) - 0.5d0*(Qfl(2)**2+Qfl(3)**2+Qfl(4)**2)&
            /Qfl(1) - 0.5d0*(Qfl(6)**2+Qfl(7)**2+Qfl(8)**2) - 0.5d0*Qfl(9)**2)

            Qfr(1) = rhor
            Qfr(2) = rhor*ur
            Qfr(3) = rhor*vr
            Qfr(4) = rhor*wr
            Qfr(5) = pinf/(gam-1)+0.5d0*rhor*(ur**2+vr**2+wr**2) &
            + 0.5d0*(Bxr**2+Byr**2+Bzr**2) + 0.5d0*Qfl(9)**2.d0
            Qfr(6) = Bxr
            Qfr(7) = Byr
            Qfr(8) = Bzr
            Qfr(9) = 0.d0
        
            if (vismode==1) then
                if(faml==1) then
                    Qvfi(1:9,ifpl,jfpl,icleft)  = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                else
                    Qvfj(1:9,ifpl,jfpl,icleft)  = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                end if
            end if

            ! now we have the Q on either side of face
            CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)

            if(faml==1) then
                F1(1:9,ifpl,jfpl,icleft)  = Fnl(1:9)
            else
                G2(1:9,ifpl,jfpl,icleft)  = Fnl(1:9)
            end if
        end do ! loop over flux point on exit boundary face

    end do ! loop over exit boundary faces

    DO ifn=1,NWALL

        iface = IBFWAL(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1

        do nfp=1,N
            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft

            sign_l = 1
            Qfl(1:9) = 0.d0

            if(faml==1) then
                do sp=1,N
                    Qs(1:9) = Q(1:9,sp,jfpl,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(ifpl,sp)
                end do
                if(ifpl==1) sign_l=-1
                normf(1:2) = S1(1,1:2,ifpl,jfpl,ic)
            else
                do sp=1,N
                    Qs(1:9) = Q(1:9,ifpl,sp,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(jfpl,sp)
                end do
                if(jfpl==1) sign_l=-1
                normf(1:2) = S2(2,1:2,ifpl,jfpl,ic)
            end if

            rhow = Qfl(1)
            uw = Qfl(2)/rhow
            vw = Qfl(3)/rhow
            ww = Qfl(4)/rhow
            Qfr(1) = rhow
            Qfr(5) = Qfl(5)
            Qfr(6) = -Qfl(6)
            Qfr(7) = -Qfl(7)
            Qfr(8) = -Qfl(8)
            Qfr(9) = Qfl(9)

            if (vismode==0) then
            
                magnorm = sqrt(normf(1)**2+normf(2)**2)
                VN = uw*normf(1)/magnorm+vw*normf(2)/magnorm

                Qfr(2) = rhow * (uw-2.d0*VN*normf(1)/magnorm)
                Qfr(3) = rhow * (vw-2.d0*VN*normf(2)/magnorm)
                Qfr(4) = -rhow * ww

            else if (vismode==1) then

                Qfr(2) = -rhow * uw
                Qfr(3) = -rhow * vw
                Qfr(4) = -rhow * ww

                if(faml==1) then
                    Qvfi(1:9,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                else
                    Qvfj(1:9,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                end if

            end if

            ! now we have the Q on either side of face
            CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)

            if(faml==1) then
                F1(1:9,ifpl,jfpl,icleft)  = Fnl(1:9)
            else
                G2(1:9,ifpl,jfpl,icleft)  = Fnl(1:9)
            end if
        end do

    END DO  !

    DO ifn=1,NSYMP

        iface = IBFSYMP(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1

        do nfp=1,N
            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft

            sign_l = 1
            Qfl(1:9) = 0.d0
            if(faml==1) then
                do sp=1,N
                    Qs(1:9) = Q(1:9,sp,jfpl,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(ifpl,sp)
                end do
                if(ifpl==1) sign_l=-1
                normf(1:2) = S1(1,1:2,ifpl,jfpl,ic)
            else
                do sp=1,N
                    Qs(1:9) = Q(1:9,ifpl,sp,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(jfpl,sp)
                end do
                if(jfpl==1) sign_l=-1
                normf(1:2) = S2(2,1:2,ifpl,jfpl,ic)
            end if

            rhow = Qfl(1)
	        uw   = Qfl(2) / rhow
	        vw   = Qfl(3) / rhow
            ww   = Qfl(4) / rhow

            Qfr(1) = rinf * 2.d0 - Qfl(1)

            magnorm = sqrt(normf(1)**2+normf(2)**2)
            VN = uw*normf(1)/magnorm+vw*normf(2)/magnorm

            Qfr(2) = Qfr(1) * (uw-2.d0*VN*normf(1)/magnorm)
	        Qfr(3) = Qfr(1) * (vw-2.d0*VN*normf(2)/magnorm)
            Qfr(4) = -Qfr(1) * ww

            pr = (gam-1.d0)*(Qfl(5) - 0.5d0*(Qfl(2)**2+Qfl(3)**2+Qfl(4)**2)&
            /Qfl(1) - 0.5d0*(Qfl(6)**2+Qfl(7)**2+Qfl(8)**2) - 0.5d0*Qfl(9)**2)
            
            BN = Qfl(6)*normf(1)/magnorm+Qfl(7)*normf(2)/magnorm
            ! Qfr(6) = Qfl(6)-2.d0*BN*normf(1)/magnorm
            ! Qfr(7) = Qfl(7)-2.d0*BN*normf(1)/magnorm
            Qfr(6) = Bxinf
            Qfr(7) = Byinf
            Qfr(8) = -Qfl(8)
            Qfr(9) = Qfl(9)

            Qfr(5) = pinf/(gam-1.d0) + &
            0.5d0/Qfr(1)*(Qfr(2)**2 + Qfr(3)**2 + Qfr(4)**2) + &
            0.5d0 * (Qfr(6)**2+Qfr(7)**2+Qfr(8)**2+Qfr(9)**2)

            if (vismode==1) then
                if(faml==1) then
                    Qvfi(1:9,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                else
                    Qvfj(1:9,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                end if
            end if
  
            CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)
  
            if(faml==1) then
                F1(1:9,ifpl,jfpl,icleft)  = Fnl(1:9)
            else
                G2(1:9,ifpl,jfpl,icleft)  = Fnl(1:9)
            end if
       
        end do

    END DO  !

    do ifn=1,NFREE

        iface = IBFREE(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1

        do nfp=1,N

            ! lets get the left solution vectors at flux point
            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft
            sign_l = 1

            if(faml==1) then
                Qfl(1:9) = 0.d0
                do sp=1,N
                    Qs(1:9)  = Q(1:9,sp,jfpl,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(ifpl,sp)
                end do
                if(ifpl==1) sign_l=-1
                normf(1:2) = S1(1,1:2,ifpl,jfpl,ic)
           else
                Qfl(1:9) = 0.d0
                do sp=1,N
                    Qs(1:9) = Q(1:9,ifpl,sp,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(jfpl,sp)
                end do
                if(jfpl==1) sign_l=-1
                normf(1:2) = S2(2,1:2,ifpl,jfpl,ic)
            end if
            
            ! Set the BC on the right sol vector
            Qfr(1:9) = Qfl(1:9)
   
            if (vismode==1) then
                if(faml==1) then
                    Qvfi(1:9,ifpl,jfpl,icleft)  = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                else
                    Qvfj(1:9,ifpl,jfpl,icleft)  = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                end if
            end if
   
            ! now we have the Q on either side of face
            CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)
   
            if(faml==1) then
                F1(1:9,ifpl,jfpl,icleft)  = Fnl(1:9)
            else
                G2(1:9,ifpl,jfpl,icleft)  = Fnl(1:9)
            end if
        end do ! loop over flux point on exit boundary face
   
    end do ! loop over exit boundary faces

    END SUBROUTINE BCFLUX


    SUBROUTINE BCVISFLUX
    use setup2d
    implicit none

    integer :: ifn,iface,icleft,ifacelc,faml,ifpl,jfpl,nfp,ic,sp
    double precision,dimension(9) :: Qfl,Qfl2,Qfr,Qfr2,Qfo,Qs
    double precision :: ur,vr,wr,inte,dur,dvr,dwr,dKE,dBp,dedxifp,dedetafp

    DO ifn=1,NINLET

        iface = IBFINL(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1
        
        do nfp=1,N
            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft
            Qfl(1:9) = 0.d0
            Qfl2(1:9) = 0.d0
            if(faml==1) then
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,sp,jfpl,ic)*Lmat(ifpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,sp,jfpl,ic)*Lmat(ifpl,sp)
                end do
            else
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,ifpl,sp,ic)*Lmat(jfpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,ifpl,sp,ic)*Lmat(jfpl,sp)
                end do
            end if
   
            Qfr(1:9) = Qfl(1:9)
            Qfr2(1:9) = Qfl2(1:9)
           
            if(faml==1) then
                nablaQvfi(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfi(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            else
                nablaQvfj(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfj(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            end if
        
        end do
   
    END DO


    DO ifn=1,NOUTLET

        iface = IBFOUT(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1
        
        do nfp=1,N

            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft
            Qfl(1:9) = 0.d0
            Qfl2(1:9) = 0.d0
           if(faml==1) then
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,sp,jfpl,ic)*Lmat(ifpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,sp,jfpl,ic)*Lmat(ifpl,sp)
                end do
            else
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,ifpl,sp,ic)*Lmat(jfpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,ifpl,sp,ic)*Lmat(jfpl,sp)
                end do
            end if
   
           Qfr(1:9) = Qfl(1:9)
           Qfr2(1:9) = Qfl2(1:9)
   
            if(faml==1) then
                nablaQvfi(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfi(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            else
                nablaQvfj(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfj(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            end if
        end do
   
    END DO  

    DO ifn=1,NWALL
      
        iface = IBFWAL(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1
        
        do nfp=1,N

            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft
            Qfl(1:9) = 0.d0
            Qfl2(1:9) = 0.d0

            Qfo(1:9) = 0.d0

            dedxifp =0.d0
            dedetafp =0.d0

            if(faml==1) then
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,sp,jfpl,ic)*Lmat(ifpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,sp,jfpl,ic)*Lmat(ifpl,sp)

                    Qs(1:9) = Q(1:9,sp,jfpl,ic)
                    Qfo(1:9) = Qfo(1:9) + Qs(1:9)*Lmat(ifpl,sp)

                    dedxifp = dedxifp + dedxi(sp,jfpl,ic)*Lmat(ifpl,sp)
                    dedetafp = dedetafp + dedeta(sp,jfpl,ic)*Lmat(ifpl,sp)

                end do

            else
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,ifpl,sp,ic)*Lmat(jfpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,ifpl,sp,ic)*Lmat(jfpl,sp)

                    Qs(1:9) = Q(1:9,ifpl,sp,ic)
                    Qfo(1:9) = Qfo(1:9) + Qs(1:9)*Lmat(jfpl,sp)
                        
                    dedxifp = dedxifp + dedxi(ifpl,sp,ic)*Lmat(jfpl,sp)
                    dedetafp = dedetafp + dedeta(ifpl,sp,ic)*Lmat(jfpl,sp)
                end do
            end if
            
            ! implement adiabatic wall boundary condition
            ur = Qfo(2)/Qfo(1)
            vr = Qfo(3)/Qfo(1)
            wr = Qfo(4)/Qfo(1)
            inte = Qfo(5)/Qfo(1) - 0.5d0*(ur**2+vr**2+wr**2) - &
            0.5d0*(Qfo(6)**2+Qfo(7)**2+Qfo(8)**2+Qfo(9)**2)/Qfo(1)
    
            Qfr(1) = Qfl(1)
            Qfr(2) = Qfl(2)
            Qfr(3) = Qfl(3)
            Qfr(4) = Qfl(4)
            Qfr(6) = Qfl(6)
            Qfr(7) = Qfl(7)
            Qfr(8) = Qfl(8)
            Qfr(9) = Qfl(9)
            dur = (Qfl(2) - ur*Qfl(1))/Qfo(1)
            dvr = (Qfl(3) - vr*Qfl(1))/Qfo(1)
            dwr = (Qfl(4) - wr*Qfl(1))/Qfo(1)
            !
            dKE = 0.5d0*(ur**2+vr**2+wr**2)*Qfr(1) + &
            Qfo(1)*(ur*dur+vr*dvr+wr*dwr)
            !
            dBp = Qfo(6)*Qfl(6)+Qfo(7)*Qfl(7)+Qfo(8)*Qfl(8)+Qfo(9)*Qfl(9)
    
            if (faml == 1) then
                Qfr(5) = -Qfl(5) +2.0*inte*Qfr(1)+2.0*dKE+dBp + &
                2.0 *dedetafp * S1(2,1,ifpl,jfpl,ic)*Qfo(1)
            else
                Qfr(5) = -Qfl(5) +2.0*inte*Qfr(1)+2.0*dKE+dBp + &
                2.0 *dedxifp * S2(1,1,ifpl,jfpl,ic)*Qfo(1)
            end if
            !
            Qfr2(1) = Qfl2(1)
            Qfr2(2) = Qfl2(2)
            Qfr2(3) = Qfl2(3)
            Qfr2(4) = Qfl2(4)
            Qfr2(6) = Qfl2(6)
            Qfr2(7) = Qfl2(7)
            Qfr2(8) = Qfl2(8)
            Qfr2(9) = Qfl2(9)
            !
            dur = (Qfl2(2) - ur*Qfl2(1))/Qfo(1)
            dvr = (Qfl2(3) - vr*Qfl2(1))/Qfo(1)
            dwr = (Qfl2(4) - wr*Qfl2(1))/Qfo(1)
            !
            dKE = 0.5d0*(ur**2+vr**2+wr**2)*Qfr2(1) + &
            Qfo(1)*(ur*dur+vr*dvr+wr*dwr)
            !
            dBp = Qfo(6)*Qfl(6)+Qfo(7)*Qfl(7)+Qfo(8)*Qfl(8)+Qfo(9)*Qfl(9)

            if (faml == 1) then
                Qfr2(5) = -Qfl2(5) +2.0*inte*Qfr2(1) +2.0*dKE + dBp +&
                2.0*dedetafp * S1(2,2,ifpl,jfpl,ic)*Qfo(1)
            else
                Qfr2(5) = -Qfl2(5) +2.0*inte*Qfr2(1) +2.0*dKE + dBp +&
                2.0*dedxifp * S2(1,2,ifpl,jfpl,ic)*Qfo(1)
            end if

    
            if(faml==1) then
                nablaQvfi(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfi(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            else
                nablaQvfj(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfj(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            end if
         end do

    end do

    DO ifn=1,NSYMP

        iface = IBFSYMP(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1

        do nfp=1,N
            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft
            Qfl(1:9) = 0.d0
            Qfl2(1:9) = 0.d0
            if(faml==1) then
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,sp,jfpl,ic)*Lmat(ifpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,sp,jfpl,ic)*Lmat(ifpl,sp)
                end do
            else
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,ifpl,sp,ic)*Lmat(jfpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,ifpl,sp,ic)*Lmat(jfpl,sp)
                end do
            end if
   
            Qfr(1:9) = Qfl(1:9)
            Qfr2(1:9) = Qfl2(1:9)
           
            if(faml==1) then
                nablaQvfi(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfi(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            else
                nablaQvfj(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfj(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            end if
        end do
   
    END DO 


    DO ifn=1,NFREE

        iface = IBFREE(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1

        do nfp=1,N

            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft
            Qfl(1:9) = 0.d0
            Qfl2(1:9) = 0.d0

            if(faml==1) then
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,sp,jfpl,ic)*Lmat(ifpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,sp,jfpl,ic)*Lmat(ifpl,sp)
                end do
            else
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,ifpl,sp,ic)*Lmat(jfpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,ifpl,sp,ic)*Lmat(jfpl,sp)
                end do
            end if
   
            Qfr(1:9) = Qfl(1:9)
            Qfr2(1:9) = Qfl2(1:9)
           
            if(faml==1) then
                nablaQvfi(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfi(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            else
                nablaQvfj(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfj(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            end if
        end do
   
    END DO  

    END SUBROUTINE BCVISFLUX