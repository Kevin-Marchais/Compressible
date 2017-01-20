module numerique

  use donnees

  implicit none

contains

  subroutine calcul_U1
    real(PR), allocatable, dimension(:,:) :: flu
    integer :: i,k
    
    allocate(flu(size(U,1),2:size(U,2)))
    
    do i=2,size(U,2)
       flu(:,i) = flux_num(U(:,i-1), U(:,i), F(:,i-1), F(:,i), bg(i), bd(i))
    end do
    
    !F(i) = F(i-1/2)
    
    do i = 2,size(U,2)-1
       do k=1,size(U,1)
          U1(k,i) = U(k,i) - dtsdx*(flu(k,i+1)-flu(k,i))
       end do
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

  subroutine physique
    integer  :: i
    
    do i=1,Nx
       rho(i) = U(1,i)
       ux(i)  = U(2,i)/U(1,i)
       E(i)   = U(3,i)
       p(i)   = (gamma-1)*(U(3,i)-0.5_PR*U(1,i)*U(2,i)**2)
       c(i)   = sqrt(gamma*p(i)/rho(i))

       if(i/=1) then
          F(1,i) = U(2,i)
          F(2,i) = U(2,i)**2/U(1,i) + p(i)
          F(3,i) = (U(3,i)+p(i))*U(2,i)/U(1,i)
       end if
    end do

    ! vitesses
    do i=2,Nx
       bg(i) = min(ux(i)-c(i),ux(i-1)-c(i-1))
       bd(i) = max(ux(i)+c(i),ux(i-1)+c(i-1))
    end do
    
  end subroutine physique

  function flux_num(Ug,Ud,Fg,Fd,bm,bp)
    real(PR), dimension(:), intent(in) :: Ug,Ud,Fg,Fd
    real(PR), dimension(size(Ug,1))    :: flux_num
    real(PR)                           :: bp,bm
    integer                            :: i

    if (bp*bm<=0) then
       do i=1,size(Ug)
          !print*,i
          !print*,bp,Fg(i),i !Fg(3) est débile
          !print*,bm*Fd(i)
          !print*,bp*bm*(Ud(i)-Ug(i))
          flux_num(i) = (bp*Fg(i)-bm*Fd(i)+bp*bm*(Ud(i)-Ug(i)))/(bp-bm)
          !print*,i
          !write(*,*)
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
  
  subroutine calcul_dt
    dt    = 0.5_PR*cfl*dx/(max(maxval(abs(bd)),maxval(abs(bg))))
    dtsdx = dt/dx
  end subroutine calcul_dt
  
end module numerique
