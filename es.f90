module es
  use donnees
  
  implicit none

contains

  subroutine lit_entree
    
    open(1,file='param.dat')
    read(1,*) Lx, Ly
    read(1,*) Nx, Ny
    read(1,*) cfl
    read(1,*) Tmax
    read(1,*) pHD, rhoHD, uxHD, uyHD !1
    print*,"HD",pHD,rhoHD,uxHD,uyHD
    read(1,*) pHG, rhoHG, uxHG, uyHG !2
    print*,"HG",pHG, rhoHG, uxHG, uyHG
    read(1,*) pBG, rhoBG, uxBG, uyBG !3
    print*,"BG",pBG, rhoBG, uxBG, uyBG
    read(1,*) pBD, rhoBD, uxBD, uyBD !4
    print*,"BD",pBD, rhoBD, uxBD, uyBD
    close(1)

    dx = Lx/Nx
    dy = Ly/Ny
    
    print*,"Lx =",Lx
    print*,"Ly =",Ly
    print*,"dx =",dx
    print*,"dy =",dy
    print*,"Temps final =",Tmax
    print*,"Nombre de mailles",Nx*Ny
    
  end subroutine lit_entree

  subroutine initialisation
    real(PR) :: x,y
    integer  :: i,j
    
    time = 0._PR
    iter = 0
    gamma = 1.4_PR
    xorder = 1

    do i=1,Nx
       x = (i-1)*dx !x(i-1/2)
       do j=1,Ny
          y = (j-1)*dy !y(j-1/2)
          if(x>0.5_PR*Lx .and. y>0.5_PR*Ly) then
             U0(1,i,j) = rhoHD
             U0(2,i,j) = rhoHD*uxHD
             U0(3,i,j) = rhoHD*uyHD
             U0(4,i,j) = pHD/(gamma-1) + 0.5_PR*rhoHD*(uxHD**2+uyHD**2)
          else if(x>0.5_PR*Lx .and. y<=0.5_PR*Ly) then
             U0(1,i,j) = rhoBD
             U0(2,i,j) = rhoBD*uxBD
             U0(3,i,j) = rhoBD*uyBD
             U0(4,i,j) = pBD/(gamma-1) + 0.5_PR*rhoBD*(uxBD**2+uyBD**2)
          else if(x<=0.5_PR*LX .and. y>0.5_PR*Ly) then
             U0(1,i,j) = rhoHG
             U0(2,i,j) = rhoHG*uxHG
             U0(3,i,j) = rhoHG*uyHG
             U0(4,i,j) = pHG/(gamma-1) + 0.5_PR*rhoHG*(uxHG**2+uyHG**2)
          else
             U0(1,i,j) = rhoBG
             U0(2,i,j) = rhoBG*uxBG
             U0(3,i,j) = rhoBG*uyBG
             U0(4,i,j) = pBG/(gamma-1) + 0.5_PR*rhoBG*(uxBG**2+uyBG**2)
          end if
       end do
    end do
    
    nbAffichages = 100
    
  end subroutine initialisation

  subroutine ecriture
    integer :: i,j

    rho = U0(1,:,:)
    ux  = U0(2,:,:)/U0(1,:,:)
    uy  = U0(3,:,:)/U0(1,:,:)
    E   = U0(4,:,:)
    p   = (gamma-1)*(U0(4,:,:)-0.5_PR*U0(1,:,:)*(U0(2,:,:)**2+U0(3,:,:)**2))
    c   = sqrt(gamma*p/rho)

    open(1,file='density.dat')
    do i=1,Nx
       do j=1,Ny
          write(1,*) (i-0.5_PR)*dx, (j-0.5_PR)*dy, rho(i,j)
       end do
       write(1,*)
    end do
    close(1)
    
  end subroutine ecriture
  
end module es
