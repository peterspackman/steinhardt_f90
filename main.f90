program main
  use iso_fortran_env, only: real64
  use spherical_harmonics
  use wigner3j_module
  implicit none

  integer, parameter :: lmax = 5
  real(real64), parameter :: pi = acos(-1.0_real64)
  type(SphericalHarmonics) :: sh
  complex(real64), allocatable :: coeffs(:)
  real(real64), allocatable :: w3ja(:)
  real(real64) :: theta, phi, w3j
  real(real64) :: pos(3)
  integer :: l, m, i

  ! Initialize the spherical harmonics evaluator
  sh = SphericalHarmonics(lmax)
  allocate(coeffs(sh%nlm()))

  ! Evaluate spherical harmonics at a specific theta and phi
  theta = pi / 3.0_real64
  phi = pi / 4.0_real64
  coeffs = sh%evaluate_cvec_theta_phi(theta, phi)

  write(*, '(A)') "Spherical Harmonics Coefficients (theta, phi):"
  do l = 0, lmax
    do m = -l, l
      i = l * (l + 1) + m + 1
      write(*, '(A, I2, A, I2, A, 2F12.6)') "Y(", l, ",", m, ") = ", coeffs(i)
    end do
  end do

  ! Evaluate spherical harmonics at a specific position vector
  pos = [0.5_real64, 0.5_real64, sqrt(2.0_real64) / 2.0_real64]
  coeffs = sh%evaluate_cvec_pos(pos)

  write(*, '(A)') "Spherical Harmonics Coefficients (position vector):"
  do l = 0, lmax
    do m = -l, l
      i = l * (l + 1) + m + 1
      write(*, '(A, I2, A, I2, A, 2F12.6)') "Y(", l, ",", m, ") = ", coeffs(i)
    end do
  end do

  w3j = wigner3j_single(2.0_real64, 2.0_real64, 2.0_real64, 0.0_real64, 0.0_real64,  0.0_real64)
  write(*, '(A, F12.6, A, F12.6)') "Expected w3j", -sqrt(2.0 / 35.0), " found: ", w3j
  w3ja = wigner3j(2.0_real64, 2.0_real64, 0.0_real64, 0.0_real64,  0.0_real64)
  print *, w3ja
      

end program main
