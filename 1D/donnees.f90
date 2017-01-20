module donnees
  implicit none

  Integer , Parameter :: PR=Selected_real_kind(12,70)
  
  integer                                 :: Nx
  real(PR), allocatable, dimension(:)     :: rho, ux, E, p, c
  real(PR), allocatable, dimension(:,:)   :: U, U1, F
  real(PR), allocatable, dimension(:)     :: bG, bD
  real(PR)                                :: dx, Lx
  real(PR)                                :: dt, dtsdx, Tmax, time
  real(PR)                                :: gamma, cfl
  real(PR)                                :: pD,rhoD,uxD
  real(PR)                                :: pG,rhoG,uxG
  
  integer                                 :: nbAffichages

end module donnees
