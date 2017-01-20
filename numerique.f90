module numerique

  use donnees

  implicit none

contains

  subroutine Euler_explicite
    real(PR), dimension(4,2:Nx,2:Ny) :: fluF, fluG
    integer :: i,j

    !print*, U0(:,228,395)
    
    do i=2,Nx
       do j=2,Ny
          fluF(:,i,j) = flux_num(U0(:,i-1,j), U0(:,i,j), F(:,i-1,j), F(:,i,j), bg(i), bd(i))
          fluG(:,i,j) = flux_num(U0(:,i,j-1), U0(:,i,j), G(:,i,j-1), G(:,i,j), bb(j), bh(j))
       end do
    end do
    
    !F(i) = F(i-1/2)

    do i = 2,Nx-1
       do j = 2,Ny-1
          U(:,i,j) = U0(:,i,j) - dtsdx*(fluF(:,i+1,j)-fluF(:,i,j)) - dtsdy*(fluG(:,i,j+1)-fluG(:,i,j))
       end do
    end do

    call CLim(U)
    
  end subroutine Euler_explicite

  subroutine RK2(alpha)
    real(PR), dimension(4,2:Nx,2:Ny)     :: fluF, fluG, fluF1, fluG1
    real(PR), dimension(4,Nx,Ny)         :: U1,F1,G1
    real(PR), dimension(4,2:Nx-1,2:Ny-1) :: k1,k2
    real(PR)                             :: alpha
    integer :: i,j

    print*,U0(:,228,395)

    !k1
    do i=2,Nx
       do j=2,Ny
          fluF(:,i,j) = flux_num(U0(:,i-1,j), U0(:,i,j), F(:,i-1,j), F(:,i,j), bg(i), bd(i))
          fluG(:,i,j) = flux_num(U0(:,i,j-1), U0(:,i,j), G(:,i,j-1), G(:,i,j), bb(j), bh(j))
       end do
    end do
    
    do i=2,Nx-1
       do j=2,Ny-1
          k1(:,i,j) = - dtsdx*(fluF(:,i+1,j)-fluF(:,i,j)) - dtsdy*(fluG(:,i,j+1)-fluG(:,i,j))
          U1(:,i,j) = U0(:,i,j) + dt*k1(:,i,j)/(2._PR*alpha)
       end do
    end do
    
    ! k2
    call CLim(U1)
    !if(iter==66)
    !print*,"U0 =",U0(3,228,395)
    !print*,"k1 =",k1(3,228,395)
    F1 = calculF(U1)
    G1 = calculG(U1)
    
    do i=2,Nx
       do j=2,Ny
          fluF1(:,i,j) = flux_num(U1(:,i-1,j), U1(:,i,j), F1(:,i-1,j), F1(:,i,j), bg(i), bd(i))
          fluG1(:,i,j) = flux_num(U1(:,i,j-1), U1(:,i,j), G1(:,i,j-1), G1(:,i,j), bb(j), bh(j))
       end do
    end do
    
    do i=2,Nx-1
       do j=2,Ny-1
          k2(:,i,j) = - dtsdx*(fluF1(:,i+1,j)-fluF1(:,i,j)) - dtsdy*(fluG1(:,i,j+1)-fluG1(:,i,j))
          U(:,i,j) = U0(:,i,j) + (1-alpha)*dt*k1(:,i,j) + alpha*dt*k2(:,i,j)
       end do
    end do
    
    call CLim(U)
    
  end subroutine RK2

  subroutine RK4
    real(PR), dimension(4,2:Nx,2:Ny) :: fluF, fluG, fluF1, fluG1, fluF2,fluG2, fluF3, fluG3
    real(PR), dimension(4,Nx,Ny)     :: U1,U2,U3,F1,F2,F3,G1,G2,G3,k1,k2,k3,k4
    real(PR)                         :: alpha
    integer :: i,j

    !k1
    do i=2,Nx
       do j=2,Ny
          fluF(:,i,j) = flux_num(U0(:,i-1,j), U0(:,i,j), F(:,i-1,j), F(:,i,j), bg(i), bd(i))
          fluG(:,i,j) = flux_num(U0(:,i,j-1), U0(:,i,j), G(:,i,j-1), G(:,i,j), bb(j), bh(j))
       end do
    end do
    
    do i=2,Nx-1
       do j=2,Ny-1
          k1(:,i,j) = - dtsdx*(fluF(:,i+1,j)-fluF(:,i,j)) - dtsdy*(fluG(:,i,j+1)-fluG(:,i,j))
          U1(:,i,j) = U0(:,i,j) + dt*k1(:,i,j)/2._PR
       end do
    end do
    
    ! k2
    
    call CLim(U1)
    F1 = calculF(U1)
    G1 = calculG(U1)
    
    do i=2,Nx
       do j=2,Ny
          fluF1(:,i,j) = flux_num(U1(:,i-1,j), U1(:,i,j), F1(:,i-1,j), F1(:,i,j), bg(i), bd(i))
          fluG1(:,i,j) = flux_num(U1(:,i,j-1), U1(:,i,j), G1(:,i,j-1), G1(:,i,j), bb(j), bh(j))
       end do
    end do
    
    do i=2,Nx-1
       do j=2,Ny-1
          k2(:,i,j) = - dtsdx*(fluF1(:,i+1,j)-fluF1(:,i,j)) - dtsdy*(fluG1(:,i,j+1)-fluG1(:,i,j))
          U2(:,i,j) = U0(:,i,j) + dt*k2(:,i,j)/2._PR
       end do
    end do
    
    ! k3
    call CLim(U2)
    F2 = calculF(U2)
    G2 = calculG(U2)

    do i=2,Nx
       do j=2,Ny
          fluF2(:,i,j) = flux_num(U2(:,i-1,j), U2(:,i,j), F2(:,i-1,j), F2(:,i,j), bg(i), bd(i))
          fluG2(:,i,j) = flux_num(U2(:,i,j-1), U2(:,i,j), G2(:,i,j-1), G2(:,i,j), bb(j), bh(j))
       end do
    end do
    
    do i=2,Nx-1
       do j=2,Ny-1
          k3(:,i,j) = - dtsdx*(fluF2(:,i+1,j)-fluF2(:,i,j)) - dtsdy*(fluG2(:,i,j+1)-fluG2(:,i,j))
          U3(:,i,j) = U0(:,i,j) + dt*k3(:,i,j)
       end do
    end do

    !k4
    call CLim(U3)
    F3 = calculF(U3)
    G3 = calculG(U3)
    
    do i=2,Nx
       do j=2,Ny
          fluF3(:,i,j) = flux_num(U3(:,i-1,j), U3(:,i,j), F3(:,i-1,j), F3(:,i,j), bg(i), bd(i))
          fluG3(:,i,j) = flux_num(U3(:,i,j-1), U3(:,i,j), G3(:,i,j-1), G3(:,i,j), bb(j), bh(j))
       end do
    end do

    do i=2,Nx-1
       do j=2,Ny-1
          k4(:,i,j) = - dtsdx*(fluF3(:,i+1,j)-fluF3(:,i,j)) - dtsdy*(fluG3(:,i,j+1)-fluG3(:,i,j))
          U(:,i,j) = U0(:,i,j) + dt/6._PR * (k1(:,i,j)+2*k2(:,i,j)+2*k3(:,i,j)+k4(:,i,j))
       end do
    end do

    ! RK4
    
    call CLim(U)
    
  end subroutine RK4

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

  function calculF(X)
    real(PR), dimension(4,Nx,Ny) :: calculF
    real(PR), dimension(4,Nx,Ny) :: X
    real(PR), dimension(Nx,Ny)   :: pression
    integer :: i,j

    do i=1,Nx
       do j=1,Ny
          !if(iter>65)print*,iter,i,j,"ok1",X(4,i,j),X(1,i,j),X(2,i,j),X(3,i,j)
          pression(i,j) = (gamma-1._PR)*(X(4,i,j)-0.5_PR*X(1,i,j)*(X(2,i,j)**2+X(3,i,j)**2))
          !if(iter>65)print*,iter,i,j,"ok2"
       end do
    end do
    
    calculF(1,:,:) = X(2,:,:)
    calculF(2,:,:) = X(2,:,:)**2/X(1,:,:) + pression
    calculF(3,:,:) = X(2,:,:)*X(3,:,:)/X(1,:,:)
    calculF(4,:,:) = (X(4,:,:)+pression)*X(2,:,:)/X(1,:,:)

  end function calculF

  function calculG(X)
    real(PR), dimension(4,Nx,Ny) :: calculG
    real(PR), dimension(4,Nx,Ny) :: X
    real(PR), dimension(Nx,Ny)   :: pression

    pression = (gamma-1._PR)*(X(4,:,:)-0.5_PR*X(1,:,:)*(X(2,:,:)**2+X(3,:,:)**2))    

    calculG(1,:,:) = X(3,:,:)
    calculG(2,:,:) = X(2,:,:)*X(3,:,:)/X(1,:,:)
    calculG(3,:,:) = X(3,:,:)**2/X(1,:,:) + pression
    calculG(4,:,:) = (X(4,:,:)+pression)*X(3,:,:)/X(1,:,:)

  end function calculG

  subroutine physique
    integer :: i,j
    
    rho = U0(1,:,:)
    ux  = U0(2,:,:)/U0(1,:,:)
    uy  = U0(3,:,:)/U0(1,:,:)
    E   = U0(4,:,:)
    p   = (gamma-1)*(U0(4,:,:)-0.5_PR*U0(1,:,:)*(U0(2,:,:)**2+U0(3,:,:)**2))
    c   = sqrt(gamma*p/rho)
    
    ! Pas optimisé calcul de p avant, puis dans calculF et dans calculG
    F = calculF(U0) 
    G = calculG(U0)

    ! vitesses
    !bg(i) = min(minval(ux(i,:)-c(i,:)),minval(ux(i-1,:)-c(i-1,:)))
    do i=2,Nx
       bg(i) = min(ux(i,1)-c(i,1),ux(i-1,1)-c(i-1,1))
       bd(i) = max(ux(i,1)+c(i,1),ux(i-1,1)+c(i-1,1))
       do j=2,Ny
          if(bg(i) > min(ux(i,j)-c(i,j),ux(i-1,j)-c(i-1,j))) bg(i)=min(ux(i,j)-c(i,j),ux(i-1,j)-c(i-1,j))
          if(bd(i) < max(ux(i,j)+c(i,j),ux(i-1,j)+c(i-1,j))) bd(i)=max(ux(i,j)+c(i,j),ux(i-1,j)+c(i-1,j))
       end do
    end do
    do j=2,Nx
       bb(j) = min(uy(1,j)-c(1,j),uy(1,j-1)-c(1,j-1))
       bh(j) = max(uy(1,j)+c(1,j),uy(1,j-1)+c(1,j-1))
       do i=2,Nx
          if(bb(j) > min(uy(i,j)-c(i,j),uy(i,j-1)-c(i,j-1))) bb(j)=min(uy(i,j)-c(i,j),uy(i,j-1)-c(i,j-1))
          if(bh(j) < max(uy(i,j)+c(i,j),uy(i,j-1)+c(i,j-1))) bh(j)=max(uy(i,j)+c(i,j),uy(i,j-1)+c(i,j-1))
       end do
    end do
    
  end subroutine physique

  function flux_num(Ug, Ud, Fg, Fd, bm, bp)
    real(PR), dimension(:), intent(in) :: Ug,Ud,Fg,Fd
    real(PR), dimension(size(Ug,1))    :: flux_num
    real(PR)                           :: bm,bp
    real(PR),dimension(size(Ug,1)) :: num,denom,div
    integer :: i

    if (bp*bm<=0) then
       flux_num = (bp*Fg-bm*Fd+bp*bm*(Ud-Ug))/(bp-bm)
    else if (bm>0) then
       flux_num = Fg
    else
       flux_num = Fd
    end if

  end function flux_num
  
  subroutine calcul_dt
    dt = 0.5_PR*cfl*min(dx/max(maxval(bg),maxval(bd)),dy/max(maxval(bb),maxval(bh)))
    dtsdx = dt/dx
    dtsdy = dt/dy
  end subroutine calcul_dt
  
end module numerique
