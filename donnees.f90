module donnees
  implicit none

  Integer , Parameter :: PR=Selected_real_kind(12,70)
  
  integer                                 :: Nx, Ny
  real(PR), allocatable, dimension(:,:)   :: rho, ux, uy, E, p, c
  real(PR), allocatable, dimension(:,:,:) :: U0, U
  real(PR), allocatable, dimension(:,:,:) :: U_g,U_d,U_h,U_b
  real(PR), allocatable, dimension(:)     :: bB, bH, bG, bD
  real(PR)                                :: dx, dy, Lx, Ly
  real(PR)                                :: dt, dtsdx, dtsdy, Tmax, time
  real(PR)                                :: gamma, cfl
  real(PR)                                :: pBD,rhoBD,uxBD,uyBD
  real(PR)                                :: pBG,rhoBG,uxBG,uyBG
  real(PR)                                :: pHD,rhoHD,uxHD,uyHD
  real(PR)                                :: pHG,rhoHG,uxHG,uyHG
  
  integer                                 :: nbAffichages, iter, xorder

end module donnees
