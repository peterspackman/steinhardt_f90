program steinhardt_calculator
  use iso_fortran_env, only: real64
  use steinhardt_module
  implicit none

  type(Steinhardt) :: stein

  real(real64), allocatable :: positions(:,:)
  real(real64), allocatable :: q(:), w(:)

  call stein%init(6)

  ! Cubic symmetry
  allocate(positions(3, 8))
  positions(:, 1) =  [ 1,  1,  1]
  positions(:, 2) =  [ 1,  1, -1]
  positions(:, 3) =  [ 1, -1,  1]
  positions(:, 4) =  [ 1, -1, -1]
  positions(:, 5) =  [-1,  1,  1]
  positions(:, 6) =  [-1,  1, -1]
  positions(:, 7) =  [-1, -1,  1]
  positions(:, 8) =  [-1, -1, -1]

  q = stein%compute_q(positions)
  print *, "Cubic Q"
  write(*, "(7F12.6)") q

  w = stein%compute_w(positions)
  print *, "Cubic W"
  write(*, "(7F12.6)") w

  deallocate(positions)

  ! Octahedral symmetry
  allocate(positions(3, 6))
  positions(:, 1) =  [ 1,  0,  0]
  positions(:, 2) =  [-1,  0,  0]
  positions(:, 3) =  [ 0,  1,  0]
  positions(:, 4) =  [ 0, -1,  0]
  positions(:, 5) =  [ 0,  0,  1]
  positions(:, 6) =  [ 0,  0, -1]

  q = stein%compute_q(positions)
  print *, "Octahedral Q"
  write(*, "(7F12.6)") q

  w = stein%compute_w(positions)
  print *, "Octahedral W"
  write(*, "(7F12.6)") w

  deallocate(positions)

  ! Tetrahedral symmetry
  allocate(positions(3, 4))
  positions(:,1) = [ 1.0_real64,  0.0_real64, -1.0_real64 / sqrt(2.0_real64) ]
  positions(:,2) = [-1.0_real64,  0.0_real64, -1.0_real64 / sqrt(2.0_real64) ]
  positions(:,3) = [ 0.0_real64,  1.0_real64,  1.0_real64 / sqrt(2.0_real64) ]
  positions(:,4) = [ 0.0_real64, -1.0_real64,  1.0_real64 / sqrt(2.0_real64) ]

  q = stein%compute_q(positions)
  print *, "Tetrahedral Q"
  write(*, "(7F12.6)") q

  w = stein%compute_w(positions)
  print *, "Tetrahedral W"
  write(*, "(7F12.6)") w

  deallocate(positions)

  ! Icosahedral symmetry
  allocate(positions(3, 12))
  positions(:,1)  = [ 0.0_real64,  1.0_real64,  0.5_real64 * (1.0_real64 + sqrt(5.0_real64)) ]
  positions(:,2)  = [ 0.0_real64,  1.0_real64, -0.5_real64 * (1.0_real64 + sqrt(5.0_real64)) ]
  positions(:,3)  = [ 0.0_real64, -1.0_real64,  0.5_real64 * (1.0_real64 + sqrt(5.0_real64)) ]
  positions(:,4)  = [ 0.0_real64, -1.0_real64, -0.5_real64 * (1.0_real64 + sqrt(5.0_real64)) ]
  positions(:,5)  = [ 1.0_real64,  0.5_real64 * (1.0_real64 + sqrt(5.0_real64)), 0.0_real64 ]
  positions(:,6)  = [ 1.0_real64, -0.5_real64 * (1.0_real64 + sqrt(5.0_real64)), 0.0_real64 ]
  positions(:,7)  = [-1.0_real64,  0.5_real64 * (1.0_real64 + sqrt(5.0_real64)), 0.0_real64 ]
  positions(:,8)  = [-1.0_real64, -0.5_real64 * (1.0_real64 + sqrt(5.0_real64)), 0.0_real64 ]
  positions(:,9)  = [ 0.5_real64 * (1.0_real64 + sqrt(5.0_real64)), 0.0_real64,  1.0_real64 ]
  positions(:,10) = [ 0.5_real64 * (1.0_real64 + sqrt(5.0_real64)), 0.0_real64, -1.0_real64 ]
  positions(:,11) = [-0.5_real64 * (1.0_real64 + sqrt(5.0_real64)), 0.0_real64,  1.0_real64 ]
  positions(:,12) = [-0.5_real64 * (1.0_real64 + sqrt(5.0_real64)), 0.0_real64, -1.0_real64 ]

  q = stein%compute_q(positions)
  print *, "Icosahedral Q"
  write(*, "(7F12.6)") q

  w = stein%compute_w(positions)
  print *, "Icosahedral W"
  write(*, "(7F12.6)") w

  deallocate(positions)

end program steinhardt_calculator
