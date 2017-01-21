program compressibles
  use donnees
  use numerique
  use es
  
  implicit none

  integer :: i,j,ordre
  
  call lit_entree

  allocate(U0(4,Nx,Ny), U(4,Nx,Ny))
  allocate(U_g(4,2:Nx,Ny), U_d(4,2:Nx,Ny), U_b(4,Nx,2:Ny), U_h(4,Nx,2:Ny))
  allocate(ux(Nx,Ny) , uy(Nx,Ny)  , p(Nx,Ny)  , E(Nx,Ny)  , rho(Nx,Ny), c(Nx,Ny))
  allocate(bB(2:Ny)  , bH(2:Ny)   , bG(2:Nx)  , bD(2:Nx))
  
  call initialisation

  ordre = 1
  
  do while(time<Tmax)
     call calcul_dt
     if(time<dt)then
        print*, "dt =",dt
     end if
     
     if(nbAffichages*(iter/nbAffichages)==iter)then
        print*,"t =",time
     end if
     print*,iter,time

     if(ordre==1) then
        call Euler_explicite ! ordre 1
     else if(ordre==2) then
        !call RK2(0.5_PR) ! ordre 2
     else if(ordre==4) then
        !call RK4 ! ordre 4
     else
        print*,"Mauvais ordre"
        exit
     end if
     
     iter = iter + 1
     time = time + dt
     U0 = U
     
  end do

  call ecriture

  deallocate(U0, U)
  deallocate(U_g,U_d,U_b,U_h)
  deallocate(ux, uy, p, E, rho, c)
  deallocate(bB, bH, bG, bD)

end program compressibles
