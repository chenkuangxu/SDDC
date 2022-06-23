    SUBROUTINE CYCLOCFLUX
    use setup2d
    implicit none
    include 'mpif.h'

    integer		::	icleft,icright,ifacelc,ifacerc,faml,famr
    integer		::	nfp,ifpl,jfpl,ifpr,jfpr,ic,K1,K2,sp,ipair,sign_l,sign_r,k
      
    double precision,dimension(9) :: Qs,Qfl,Qfr
    double precision,dimension(9) :: Fnl,Fnr
    double precision, dimension(2) :: normf
    double precision :: eigv,Az
      
    ! calculate the interface fluxes at all interior faces (non-boundary faces)
    do ipair=1,NPAIR
      
        K1 = CYC2FPAIR(ipair,1)
        K2 = CYC2FPAIR(ipair,2)
      
        icleft  = IF2C(K1,1)
        icright = IF2C(K2,1)
      
        ifacelc = IF2C(K1,3)
        ifacerc = IF2C(K2,3)
      
        faml = mod(ifacelc,2) + 1
        famr = mod(ifacerc,2) + 1
      
        ! we now need to compute the rusanov flux at N flux points on the face
      
        do nfp=1,N
      
            ! lets get the left solution vectors at flux point
      
            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft
      
            sign_l = 1
            Qfl(1:9) = 0.d0
            if(faml==1) then
                do sp=1,N
                    Qs(1:9)	= Q(1:9,sp,jfpl,ic)
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
      
            ! lets get the right solution vector at flux point
      
            ifpr = iface2fp(N-nfp+1,ifacerc)
            jfpr = jface2fp(N-nfp+1,ifacerc)
            ic = icright
      
            sign_r = 1
            Qfr(1:9) = 0.d0
            if(famr==1) then
                do sp=1,N
                    Qs(1:9)	= Q(1:9,sp,jfpr,ic)
                    Qfr(1:9) = Qfr(1:9) + Qs(1:9)*Lmat(ifpr,sp)
                end do
                if(ifpr==N+1) sign_r=-1
            else
                do sp=1,N
                    Qs(1:9) = Q(1:9,ifpr,sp,ic)
                    Qfr(1:9) = Qfr(1:9) + Qs(1:9)*Lmat(jfpr,sp)
                end do
                if(jfpr==N+1) sign_r=-1
            end if
      
            IF (vismode==1) THEN
                if(faml==1) then
                    Qvfi(1:9,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                else
                    Qvfj(1:9,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                end if
      
                if(famr==1) then
                    Qvfi(1:9,ifpr,jfpr,icright) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                else
                    Qvfj(1:9,ifpr,jfpr,icright) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                end if
            END IF
              
            ! now we have the Q on either side of face
            CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,normf,sign_l,sign_r,eigv)
      
            if(faml==1) then
                F1(1:9,ifpl,jfpl,icleft) = Fnl(1:9)
            else
                G2(1:9,ifpl,jfpl,icleft) = Fnl(1:9)
            end if
      
            if(famr==1) then
                F1(1:9,ifpr,jfpr,icright) = Fnr(1:9)
            else
                G2(1:9,ifpr,jfpr,icright) = Fnr(1:9)
            end if
      
        end do	! do loop over points on interface
      
    end do	! do loop over interior faces

    END SUBROUTINE CYCLOCFLUX

    SUBROUTINE CYCLOCVISFLUX
    use setup2d
    implicit none
    include 'mpif.h'
        
    integer		::	icleft,icright,ifacelc,ifacerc,faml,famr
    integer		::	nfp,ifpl,jfpl,ifpr,jfpr,ic,K1,K2,sp,ipair,k
              
     double precision, dimension(9)	:: Qfl,Qfr,Qfl2,Qfr2
              
    ! calculate the interface fluxes at all interior faces (non-boundary faces)
    do ipair=1,NPAIR
              
        K1 = CYC2FPAIR(ipair,1)
        K2 = CYC2FPAIR(ipair,2)
              
        icleft  = IF2C(K1,1)
        icright = IF2C(K2,1)
              
        ifacelc = IF2C(K1,3)
        ifacerc = IF2C(K2,3)
              
        faml = mod(ifacelc,2) + 1
        famr = mod(ifacerc,2) + 1
              
        ! we now need to compute the rusanov flux at N flux points on the face
        do nfp=1,N
              
            ! lets get the left solution vectors at flux point
              
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
              
            ! lets get the right solution vector at flux point
              
            ifpr = iface2fp(N-nfp+1,ifacerc)
            jfpr = jface2fp(N-nfp+1,ifacerc)
            ic = icright
              
            Qfr(1:9) = 0.d0
            Qfr2(1:9) = 0.d0
    
            if(famr==1) then
                do sp=1,N
                    Qfr(1:9) = Qfr(1:9) + nablaQs(1:9,1,sp,jfpr,ic)*Lmat(ifpr,sp)
                    Qfr2(1:9) = Qfr2(1:9) + nablaQs(1:9,2,sp,jfpr,ic)*Lmat(ifpr,sp)
                end do
            else
                do sp=1,N
                    Qfr(1:9) = Qfr(1:9) + nablaQs(1:9,1,ifpr,sp,ic)*Lmat(jfpr,sp)
                    Qfr2(1:9) = Qfr2(1:9) + nablaQs(1:9,2,ifpr,sp,ic)*Lmat(jfpr,sp)
                end do
            end if
              
            if(faml==1) then
                nablaQvfi(1:9,1,ifpl,jfpl,icleft) = 0.5*Qfl(1:9)+0.5*Qfr(1:9)
                nablaQvfi(1:9,2,ifpl,jfpl,icleft) = 0.5*Qfl2(1:9)+0.5*Qfr2(1:9)
            else
                nablaQvfj(1:9,1,ifpl,jfpl,icleft) = 0.5*Qfl(1:9)+0.5*Qfr(1:9)
                nablaQvfj(1:9,2,ifpl,jfpl,icleft) = 0.5*Qfl2(1:9)+0.5*Qfr2(1:9)
            end if
    
            if(famr==1) then
                nablaQvfi(1:9,1,ifpr,jfpr,icright) = 0.5*Qfl(1:9)+0.5*Qfr(1:9)
                nablaQvfi(1:9,2,ifpr,jfpr,icright) = 0.5*Qfl2(1:9)+0.5*Qfr2(1:9)
            else
                nablaQvfj(1:9,1,ifpr,jfpr,icright) = 0.5*Qfl(1:9)+0.5*Qfr(1:9)
                nablaQvfj(1:9,2,ifpr,jfpr,icright) = 0.5*Qfl2(1:9)+0.5*Qfr2(1:9)
            end if
              
        end do	! do loop over points on interface
              
    end do	! do loop over interior faces
        
    END SUBROUTINE CYCLOCVISFLUX