program compressibles
  use donnees
  use numerique
  use es
  
  implicit none

  integer :: iter
  
  call lit_entree

  allocate(U(3,Nx), U1(3,Nx), F(3,Nx))
  allocate(ux(Nx), p(Nx), E(Nx), rho(Nx), c(Nx))
  allocate(bG(2:Nx), bD(2:Nx))
  
  call initialisation
  
  iter = 0
  do while(time<Tmax)
     call calcul_dt
     call calcul_U1
     iter = iter + 1
     time = time + dt
     U = U1
     call physique
  end do

  call ecriture

  deallocate(U, U1, F)
  deallocate(ux, p, E, rho, c)
  deallocate(bG, bD)

end program compressibles
