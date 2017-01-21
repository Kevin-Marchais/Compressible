module numerique

  use donnees
  use ordre2

  implicit none

contains

  subroutine Euler_explicite
    real(PR), dimension(4,2:Nx,2:Ny) :: fluF, fluG
    integer :: i,j

    call reconstruction

    do i=2,Nx
       do j=2,Ny
          !if(iter==1 .and. i==2 .and. j==202) print*,iter,i,j,U_d(:,i,j)
          fluF(:,i,j) = flux_HLLx(U_g(:,i,j), U_d(:,i,j))
          fluG(:,i,j) = flux_HLLy(U_b(:,i,j), U_h(:,i,j))
       end do
    end do
    
    !F(i) = F(i-1/2)

    do i = 2,Nx-1
       do j = 2,Ny-1
           if(iter==0 .and. i==2 .and. j==202) then
              !print*,iter,i,j,dtsdy!fluG(:,i,j+1)!U0(:,i,j)
           end if
          U(:,i,j) = U0(:,i,j) - dtsdx*(fluF(:,i+1,j)-fluF(:,i,j)) - dtsdy*(fluG(:,i,j+1)-fluG(:,i,j))
         
       end do
    end do

    call CLim(U)
    
  end subroutine Euler_explicite

  subroutine CLim(CL)
    integer  :: i,j
    real(PR) :: x,y
    real(PR),dimension(4,Nx,Ny) :: CL
    
    do i=1,Nx
       x = (i-1)*dx !x(i-1/2)
       if(x>0.5_PR*Lx) then
          CL(1,i,1)  = rhoBD
          CL(2,i,1)  = rhoBD*uxBD
          CL(3,i,1)  = rhoBD*uyBD
          CL(4,i,1)  = pBD/(gamma-1) + 0.5_PR*rhoBD*(uxBD**2+uyBD**2)

          CL(1,i,Ny) = rhoHD
          CL(2,i,Ny) = rhoHD*uxHD
          CL(3,i,Ny) = rhoHD*uyHD
          CL(4,i,Ny) = pHD/(gamma-1) + 0.5_PR*rhoHD*(uxHD**2+uyHD**2)
       else
          CL(1,i,1)  = rhoBG
          CL(2,i,1)  = rhoBG*uxBG
          CL(3,i,1)  = rhoBG*uyBG
          CL(4,i,1)  = pBG/(gamma-1) + 0.5_PR*rhoBG*(uxBG**2+uyBG**2)

          CL(1,i,Ny) = rhoHG
          CL(2,i,Ny) = rhoHG*uxHG
          CL(3,i,Ny) = rhoHG*uyHG
          CL(4,i,Ny) = pHG/(gamma-1) + 0.5_PR*rhoHG*(uxHG**2+uyHG**2)
       end if
    end do

    do j=1,Ny
       y = (j-1)*dy !y(j-1/2)
       if(y>0.5_PR*Ly) then
          CL(1,1,j) = rhoHG
          CL(2,1,j) = rhoHG*uxHG
          CL(3,1,j) = rhoHG*uyHG
          CL(4,1,j) = pHG/(gamma-1) + 0.5_PR*rhoHG*(uxHG**2+uyHG**2)

          CL(1,Nx,j) = rhoHD
          CL(2,Nx,j) = rhoHD*uxHD
          CL(3,Nx,j) = rhoHD*uyHD
          CL(4,Nx,j) = pHD/(gamma-1) + 0.5_PR*rhoHD*(uxHD**2+uyHD**2)
       else
          CL(1,1,j) = rhoBG
          CL(2,1,j) = rhoBG*uxBG
          CL(3,1,j) = rhoBG*uyBG
          CL(4,1,j) = pBG/(gamma-1) + 0.5_PR*rhoBG*(uxBG**2+uyBG**2)

          CL(1,Nx,j) = rhoBD
          CL(2,Nx,j) = rhoBD*uxBD
          CL(3,Nx,j) = rhoBD*uyBD
          CL(4,Nx,j) = pBD/(gamma-1) + 0.5_PR*rhoBD*(uxBD**2+uyBD**2)
       end if
    end do
    
  end subroutine CLim

  function F(X)
    real(PR), dimension(4) :: F
    real(PR), dimension(4) :: X
    real(PR)   :: pression
    integer :: i,j

    pression = (gamma-1._PR)*(X(4)-0.5_PR*X(1)*(X(2)**2+X(3)**2))

    
    F(1) = X(2)
    F(2) = X(2)**2/X(1) + pression
    F(3) = X(2)*X(3)/X(1)
    F(4) = (X(4)+pression)*X(2)/X(1)

  end function F

  function G(X)
    real(PR), dimension(4) :: G
    real(PR), dimension(4) :: X
    real(PR)  :: pression

    pression = (gamma-1._PR)*(X(4)-0.5_PR*X(1)*(X(2)**2+X(3)**2))    

    G(1) = X(3)
    G(2) = X(2)*X(3)/X(1)
    G(3) = X(3)**2/X(1) + pression
    G(4) = (X(4)+pression)*X(3)/X(1)

  end function G
  
  function flux_HLLx(Ug, Ud)
    real(PR), dimension(size(U,1)), intent(in) :: Ug,Ud
    real(PR), dimension(size(U,1))     :: flux_HLLx
    real(PR)                           :: bm,bp,cd,cg
    integer :: i

    cg = sqrt(gamma*(gamma-1._PR)*(Ug(4)-0.5_PR*Ug(1)*(Ug(2)**2+Ug(3)**2))/Ug(1))
    !print*,iter, Ud
    cd = sqrt(gamma*(gamma-1._PR)*(Ud(4)-0.5_PR*Ud(1)*(Ud(2)**2+Ud(3)**2))/Ud(1))
    
    bp = max(Ug(2)/Ug(1) + cg, Ud(2)/Ud(1) + cd)
    bm = min(Ug(2)/Ug(1) - cg, Ud(2)/Ud(1) - cd)

    if (bp<=0._PR .and. bm>=0._PR) then
       flux_HLLx = (bp*F(Ug)-bm*F(Ud)+bp*bm*(Ud-Ug))/(bp-bm)
    else if (bm>0) then
       flux_HLLx = F(Ug)
    else
       flux_HLLx = F(Ud)
    end if

  end function flux_HLLx

  function flux_HLLy(Ub, Uh)
    real(PR), dimension(size(U,1)), intent(in) :: Uh, Ub
    real(PR), dimension(size(U,1))    :: flux_HLLy
    real(PR)                           :: bm,bp,cb,ch
    integer :: i

    cb = sqrt(gamma*(gamma-1._PR)*(Ub(4)-0.5_PR*Ub(1)*(Ub(2)**2+Ub(3)**2))/Ub(1))
    ch = sqrt(gamma*(gamma-1._PR)*(Uh(4)-0.5_PR*Uh(1)*(Uh(2)**2+Uh(3)**2))/Uh(1))
    
    bp = max(Ub(2)/Ub(1) + cb, Uh(2)/Uh(1) + ch)
    bm = min(Ub(2)/Ub(1) - cb, Uh(2)/Uh(1) - ch)

    if (bp*bm<=0) then
       flux_HLLy = (bp*G(Ub)-bm*G(Uh)+bp*bm*(Uh-Ub))/(bp-bm)
    else if (bm>0) then
       flux_HLLy = F(Ub)
    else
       flux_HLLy = F(Uh)
    end if

  end function flux_HLLy
  
  subroutine calcul_dt
    real(PR),dimension(Nx,Ny) :: pression,son

    
    
    dt = 0.5_PR*cfl*min(dx/max(maxval(bg),maxval(bd)),dy/max(maxval(bb),maxval(bh)))
    !dt = 0.0000000002_PR
    dtsdx = dt/dx
    dtsdy = dt/dy
  end subroutine calcul_dt
  
end module numerique
