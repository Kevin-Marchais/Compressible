module ordre2
  
  use donnees

  implicit none

contains

  subroutine reconstruction!(Ug,Ud,Uh,Ub)
    real(PR),dimension(4,2:Nx,Ny) :: thetaX!,Ug,Ud
    real(PR),dimension(4,Nx,2:Ny) :: thetaY!,Uh,Ub
    integer :: i,j

    if(xorder==1)then
       do i=1,Nx
          do j=1,Ny
             if(i/=1)then
                U_d(:,i,j) = U0(:,i,j)
                U_g(:,i,j) = U0(:,i-1,j)
             end if
             if (j/=1)then
                U_h(:,i,j) = U0(:,i,j)
                U_b(:,i,j) = U0(:,i,j-1)
             end if
          end do
       end do
    else if(xorder==2)then

       do i=1,Nx
          do j=1,Ny
             if(i/=1) then
                thetaX(:,i,j) = (U0(:,i+1,j)-U0(:,i,j))/(U0(:,i,j)-U0(:,i-1,j))
             end if
             if(j/=1)then
                thetaY(:,i,j) = (U0(:,i,j+1)-U0(:,i,j))/(U0(:,i,j)-U0(:,i,j-1))
             end if
          end do
       end do
       
       
       do i=1,Nx
          do j=1,Ny
             if(i/=1)then
             U_d(:,i,j) = U0(:,i,j)   - 0.5_PR*dx*dudx(i-1,j)*phi(thetaX(:,i,j))
             U_g(:,i,j) = U0(:,i-1,j) + 0.5_PR*dx*dudx(i-1,j)*phi(thetaX(:,i,j))
          end if
          if(j/=1)then
             U_h(:,i,j) = U0(:,i,j)   - 0.5_PR*dy*dudy(i,j-1)*phi(thetaY(:,i,j))
             U_b(:,i,j) = U0(:,i,j-1) + 0.5_PR*dy*dudy(i,j-1)*phi(thetaY(:,i,j))
          end if
          end do
       end do
       
    end if

  end subroutine reconstruction

  function phi(a)
    real(PR),dimension(4) :: a,phi
    integer :: i

    do i=1,4
       if(a(i)<0._PR) then
          phi(i)=0._PR
       else
          phi(i) = min(1._PR,a(i))
       end if
    end do
    
  end function phi

  function dudx(i,j)
    integer :: i,j
    real(PR),dimension(size(U0,1)) :: dudx

    if (i/=1 .and. i/=Nx) then
       dudx = (U0(:,i+1,j)-U0(:,i-1,j))/(2*dx)
    else if(i==1)then
       dudx = (U0(:,i+1,j)-U0(:,i,j))/dx
    else
       dudx = (U0(:,i,j)-U0(:,i-1,j))/dx
    end if
    
  end function dudx

  function dudy(i,j)
    integer :: i,j
    real(PR),dimension(size(U0,1)) :: dudy

    if (j/=1 .and. j/=Ny) then
       dudy = (U0(:,i,j+1)-U0(:,i,j-1))/(2*dy)
    else if(j==1)then
       dudy = (U0(:,i,j+1)-U0(:,i,j))/dy
    else
       dudy = (U0(:,i,j)-U0(:,i,j-1))/dy
    end if
    
  end function dudy
  
end module ordre2
