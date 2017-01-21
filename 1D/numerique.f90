module numerique

  use donnees

  implicit none

contains

  subroutine calcul_U1
    real(PR), allocatable, dimension(:,:) :: flu
    integer :: i,k
    
    allocate(flu(size(U,1),2:size(U,2)))
    
    do i=2,size(U,2)
       flu(:,i) = flux_num(UG(:,i), UD(:,i-1), F(UG(:,i)), F(UD(:,i-1)), bg(i), bd(i))
    end do
    
    !F(i) = F(i-1/2)
    
    do i = 2,size(U,2)-1
       U1(:,i) = U(:,i) - dtsdx*(flu(:,i+1)-flu(:,i))
    end do

    call CLim

    deallocate(flu)
    
  end subroutine calcul_U1

  subroutine CLim
    
    U1(1,Nx) = rhoD
    U1(2,Nx) = rhoD*uxD
    U1(3,Nx) = pD/(gamma-1) + 0.5_PR*rhoD*uxD**2

    U1(1,1) = rhoG
    U1(2,1) = rhoG*uxG
    U1(3,1) = pG/(gamma-1) + 0.5_PR*rhoG*uxG**2
    
  end subroutine CLim

  subroutine vitesses(XG,XD,bm,bp)
    real(PR),dimension(size(U,1),size(U,2)),intent(in) :: XG,XD
    real(PR),dimension(2:Nx),intent(out) :: bp,bm
    integer  :: i
    real(PR),dimension(size(U,2)) :: pressionG, pressionD, c_sonG, c_sonD
    if(iter==10)print*,XD(1,102),XD(2,102),XD(3,102)
    do i=1,Nx
       pressionG(i) = (gamma-1._PR)*(XG(3,i)-0.5_PR*XG(1,i)*XG(2,i)**2)
       c_sonG(i) = sqrt(gamma*pressionG(i)/XG(1,i))

       !if(iter==1303)print*,iter,i,XD(1,i),XD(2,i),XD(3,i)
       pressionD(i) = (gamma-1._PR)*(XD(3,i)-0.5_PR*XD(1,i)*XD(2,i)**2)
       c_sonD(i) = sqrt(gamma*pressionD(i)/XD(1,i))
    end do
    
    ! vitesses
    do i=2,Nx
       bm(i) = min(XG(2,i)/XG(1,i)-c_sonG(i),XD(2,i-1)/XD(1,i-1)-c_sonD(i-1),0._PR)
       bp(i) = max(XG(2,i)/XG(1,i)+c_sonG(i),XD(2,i-1)/XD(1,i-1)+c_sonD(i-1),0._PR)
    end do
    
  end subroutine vitesses

  function F(X)
    real(PR),dimension(size(U,1)) :: F,X
    real(PR) :: pression
    
    pression = (gamma-1._PR)*(X(3)-0.5_PR*X(1)*X(2)**2)
    
    F(1) = X(2)
    F(2) = X(2)**2/X(1) + pression
    F(3) = (X(3)+pression)*X(2)/X(1)
    
  end function F

  function flux_num(Ug,Ud,Fg,Fd,bm,bp)
    real(PR), dimension(:), intent(in) :: Ug,Ud,Fg,Fd
    real(PR), dimension(size(Ug,1))    :: flux_num
    real(PR),intent(in)                :: bp,bm
    integer                            :: i

    if (bp*bm<=0) then
       do i=1,size(Ug)
          flux_num(i) = (bp*Fg(i)-bm*Fd(i)+bp*bm*(Ud(i)-Ug(i)))/(bp-bm)
       end do
    else if (bm>0) then
       do i=1,size(Ug)
          flux_num(i) = Fg(i)
       end do
    else
       do i=1,size(Ug)
          flux_num(i) = Fd(i)
       end do
    end if

  end function flux_num

  subroutine reconstruction
    integer :: i
    
    if (xorder==1)then
       UG(:,1) = U(:,1)
       UD(:,1) = U(:,1)
       do i=2,Nx
          UG(:,i) = U(:,i-1)
          UD(:,i) = U(:,i)
       end do
    else if (xorder==2) then
       do i=1,Nx
          UG(:,i) = U(:,i) - 0.5_PR*dx*dudx(i)
          UD(:,i) = U(:,i) + 0.5_PR*dx*dudx(i)
       end do
    end if
    
  end subroutine reconstruction

  function dudx(i)
    real(PR),dimension(size(U,1)) :: dudx
    integer :: i

    if(i/=1 .and. i/=Nx)then
       dudx = minmod(theta*(U(:,i+1)-U(:,i))/dx,(U(:,i+1)-U(:,i-1))/(2*dx),theta*(U(:,i)-U(:,i-1))/dx)
    else if(i==1)then
       dudx = theta*(U(:,i+1)-U(:,i))/dx
    else
       dudx = theta*(U(:,i)-U(:,i-1))/dx
    end if
    
  end function dudx

  function minmod(a1,a2,a3)
    real(PR),dimension(size(U,1)) :: minmod
    real(PR),dimension(size(U,1)) :: a1,a2,a3
    integer :: i

    !do i=1,size(U,1)
    minmod = min(abs(a1),abs(a2),abs(a3))
    !end do
    
  end function minmod
  
  subroutine calcul_dt
    dt    = 0.5_PR*cfl*dx/(max(maxval(abs(bd)),maxval(abs(bg))))
    dtsdx = dt/dx
  end subroutine calcul_dt

  subroutine config1

    rhoG = 1._PR
    uxG  = 0._PR
    pG   = 1._PR

    rhoD = 0.125_PR
    uxD  = 0._PR
    pD   = 0.1_PR

  end subroutine config1
  
end module numerique
