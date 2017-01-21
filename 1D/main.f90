program compressibles
  use donnees
  use numerique
  use es
  
  implicit none
  
  call lit_entree

  allocate(U(3,Nx), U1(3,Nx), UG(3,Nx), UD(3,Nx))
  allocate(ux(Nx), p(Nx), E(Nx), rho(Nx), c(Nx))
  allocate(bG(2:Nx), bD(2:Nx))
  allocate(theta(3))

  call config1
  
  call initialisation
  
  iter = 0
  do while(time<Tmax)
     call calcul_dt
     call calcul_U1
     iter = iter + 1
     print*,iter
     time = time + dt
     U = U1
     call reconstruction
     call vitesses(UG,UD,bg,bd)
  end do

  call ecriture

  deallocate(U, U1, UG, UD)
  deallocate(ux, p, E, rho, c)
  deallocate(bG, bD)
  deallocate(theta)

end program compressibles
