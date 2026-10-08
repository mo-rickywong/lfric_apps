!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
module calc_exner_pointwise_mod

use constants_mod,     only : r_def

implicit none

private

public :: calc_exner_pointwise
public :: linear_calc_exner_pointwise
public :: calc_pressure_pointwise

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> @brief Function to compute the exner pressure from the nonlinear equation of state
!> @details Compute the exner pressure from the equation of state:
!>          exner = ( Rd/p0 * rho * theta ) ^ (  kappa / ( 1 - kappa ) )
!! @param[in] rho   Density perturbation
!! @param[in] theta Potential temperature perturbation
!! @param[in] rd
!! @param[in] kappa
!! @param[in] p_zero
!! @return    exner Pressure perturbation
function calc_exner_pointwise(rho, theta, rd, kappa, p_zero) result(exner)

  implicit none

  real(kind=r_def)              :: exner
  real(kind=r_def), intent(in)  :: Rd
  real(kind=r_def), intent(in)  :: kappa
  real(kind=r_def), intent(in)  :: p_zero
  real(kind=r_def), intent(in)  :: rho, theta

   exner = ( Rd * rho * theta/p_zero ) ** (  kappa / ( 1.0_r_def - kappa ) )

end function calc_exner_pointwise

!> @brief Function to compute the exner pressure from the linear equation of state
!> @details Compute the exner pressure from the equation of state:
!>           exner = kappa / ( 1- kappa ) * exner_s * ( rho/rho_s + theta/theta_s )
!> @deprecated The usefulness of the linear model is to be revaluated at
!>             the end of the Gung-Ho project and removed if possible
!! @param[in] rho     Density perturbation
!! @param[in] theta   Potential temperature perturbation
!! @param[in] exner_s Reference exner pressure
!! @param[in] rho_s   Reference density
!! @param[in] theta_s Reference potential temperature
!! @return    exner   Pressure perturbation
function linear_calc_exner_pointwise(rho, theta, exner_s, rho_s, theta_s, kappa) result(exner)

  implicit none

  real(kind=r_def)              :: exner
  real(kind=r_def), intent(in)  :: rho, theta, exner_s, rho_s, theta_s, kappa

  exner = kappa / ( 1.0_r_def - kappa ) * exner_s * ( rho/rho_s + theta/theta_s )

end function linear_calc_exner_pointwise

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> @brief Function to compute the full pressure from the nonlinear equation of state
!> @details Compute the full pressure from the equation of state:
!>          pressure = ( Rd/p0 * rho * theta ) ^ (  kappa / ( 1 - kappa ) )
!! @param[in] rho   Density perturbation
!! @param[in] theta Potential temperature perturbation
!! @param[in] rd
!! @param[in] kappa
!! @param[in] p_zero
!! @return    pressure Pressure perturbation
function calc_pressure_pointwise(rho, theta, rd, kappa, p_zero) result(pressure)

  implicit none

  real(kind=r_def)              :: pressure
  real(kind=r_def), intent(in)  :: rho, theta
  real(kind=r_def), intent(in)  :: Rd
  real(kind=r_def), intent(in)  :: kappa
  real(kind=r_def), intent(in)  :: p_zero

  pressure = p_zero * ( Rd/p_zero * rho * theta ) ** (  1.0_r_def / ( 1.0_r_def - kappa ) )

end function calc_pressure_pointwise

end module calc_exner_pointwise_mod
