!------------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Routines for computing horizontal Nirvana and PPM reconstructions
!!          of fields, for use in FFSL
!------------------------------------------------------------------------------
module subgrid_horizontal_support_mod

use constants_mod,                  only: i_def, r_tran, l_def, EPS_R_TRAN
use transport_enumerated_types_mod, only: horizontal_monotone_strict,   &
                                          horizontal_monotone_relaxed,  &
                                          horizontal_monotone_positive, &
                                          horizontal_monotone_qm_pos,   &
                                          horizontal_monotone_none

implicit none

private

! Edge interpolation routines
public :: fourth_order_horizontal_edge
public :: nirvana_horizontal_edge
public :: monotonic_horizontal_edge
public :: subgrid_quadratic_recon
public :: bound_field
! Routines for special edge treatment
public :: horizontal_nirvana_case
public :: horizontal_ppm_case
public :: nirvana_special_edge
public :: linear_special_edge
public :: fourth_order_special_edge
public :: fourth_nirvana_special_edge

contains

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal PPM to estimate the quadratic subgrid representation
  !!         of a field. The function is passed five field values from consecutive cells
  !!         (which all lie in the same direction) with the dofmap
  !!         | 1 | 2 | 3 | 4 | 5 | and returns the estimated field value between
  !!         cells 2 and 3 (edge_left) and cells 3 and 4 (edge_right).
  !!         The cells are assumed to be uniform in spacing.
  !!
  !! @param[in]   field            Has dof map of the form | 1 | 2 | 3 | 4 | 5 |
  !! @param[out]  edge_left        Field value at left edge of cell
  !! @param[out]  edge_right       Field value at right edge of cell
  !----------------------------------------------------------------------------
  subroutine fourth_order_horizontal_edge(field, edge_left, edge_right)

    implicit none

    real(kind=r_tran), intent(in)  :: field(1:5)
    real(kind=r_tran), intent(out) :: edge_left, edge_right

    ! As the cell widths are assumed to be constant the edge value reduces to that given in
    ! Colella and Woodward, JCP, 54, 1984, equation (1.9)
    edge_left  = ( 7.0_r_tran/12.0_r_tran ) * ( field(2) + field(3) ) &
                -( 1.0_r_tran/12.0_r_tran ) * ( field(1) + field(4) )
    edge_right = ( 7.0_r_tran/12.0_r_tran ) * ( field(3) + field(4) ) &
                -( 1.0_r_tran/12.0_r_tran ) * ( field(2) + field(5) )

  end subroutine fourth_order_horizontal_edge

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal Nirvana to estimate the quadratic subgrid representation
  !!         of a field. The function is passed three field values from consecutive cells
  !!         (which all lie in the same direction) with the dofmap
  !!         | 1 | 2 | 3 | and returns the estimated field value between
  !!         cells 1 and 2 (edge_left) and cells 2 and 3 (edge_right).
  !!         The cells are assumed to be uniform in spacing.
  !!
  !! @param[in]   field            Has dof map of the form | 1 | 2 | 3 |
  !! @param[out]  edge_left        Field value at left edge of cell
  !! @param[out]  edge_right       Field value at right edge of cell
  !----------------------------------------------------------------------------
  subroutine nirvana_horizontal_edge(field, edge_left, edge_right)

    implicit none

    real(kind=r_tran), intent(in)  :: field(1:3)
    real(kind=r_tran), intent(out) :: edge_left, edge_right

    ! The Nirvana edges can be derived from Leonard et al. 1995.
    edge_left  = ( -field(3) + 5.0_r_tran*field(2) + 2.0_r_tran*field(1) )/6.0_r_tran
    edge_right = ( -field(1) + 5.0_r_tran*field(2) + 2.0_r_tran*field(3) )/6.0_r_tran

  end subroutine nirvana_horizontal_edge

  !----------------------------------------------------------------------------
  !> @brief  Applies monotonicty constraints to the edge values needed
  !!         for the quadratic subgrid reconstruction (i.e. PPM or Nirvana).
  !!         The three field values have dof map | 1 | 2 | 3 |
  !!         and limit the edge values on either side of cell 2.
  !!
  !! @param[in]     field          Has dof map of the form | 1 | 2 | 3 |
  !! @param[in]     monotone       Monotone option to ensures no over/undershoots
  !> @param[in]     min_val        Minimum value to enforce edge value to be for
  !!                               quasi-monotone limiter
  !! @param[inout]  edge_left      Field value at left edge of cell
  !! @param[inout]  edge_right     Field value at right edge of cell
  !----------------------------------------------------------------------------
  subroutine monotonic_horizontal_edge(field, monotone, min_val, edge_left, edge_right)

    implicit none

    real(kind=r_tran),   intent(in)    :: field(1:3)
    integer(kind=i_def), intent(in)    :: monotone
    real(kind=r_tran),   intent(in)    :: min_val
    real(kind=r_tran),   intent(inout) :: edge_left, edge_right

    real(kind=r_tran) :: t1

    if ( monotone == horizontal_monotone_strict .OR. &
         monotone == horizontal_monotone_relaxed ) then
      ! Limit left edge value to be bounded by neighbouring values
      t1 = ( edge_left - field(1) )*( field(2) - edge_left )
      if ( t1 < 0.0_r_tran ) then
         call bound_field(edge_left, field(1), field(2))
      end if
      ! Limit right edge value to be bounded by neighbouring values
      t1 = ( edge_right - field(2) )*( field(3) - edge_right )
      if ( t1 < 0.0_r_tran ) then
         call bound_field(edge_right, field(2), field(3))
      end if
    else if ( monotone == horizontal_monotone_positive ) then
      ! Positivity - ensure edge value is positive
      edge_left  = max(edge_left, 0.0_r_tran)
      edge_right = max(edge_right, 0.0_r_tran)
    else if ( monotone == horizontal_monotone_qm_pos ) then
      ! Quasi-monotone limiter based on looking at successive gradients
      t1 =( field(3)-field(2) )*(field(2)-field(1))
      if ( t1 < 0.0_r_tran) then
         call bound_field(edge_left, field(1), field(2))
         call bound_field(edge_right, field(2), field(3))
      end if
      ! Ensure edge value is greater than min_val
      edge_left  = max(edge_left, min_val)
      edge_right = max(edge_right, min_val)
    end if

  end subroutine monotonic_horizontal_edge

  !----------------------------------------------------------------------------
  !> @brief  Returns the horizontal subgrid quadratic reconstruction.
  !!         This can be used to compute the flux as:
  !!         flux = u * reconstruction
  !!         The reconstruction depends upon the cell edge values,
  !!         the field value, and the departure distance.
  !!         The reconstruction is third-order, and is based on the quadratic
  !!         subgrid reconstruction of PPM and Nirvana.
  !!
  !! @param[out]  reconstruction The PPM reconstruction
  !! @param[in]   dep            The fractional departure distance for the
  !!                             reconstruction point.
  !!                             For dep>0 the recon is on the right edge of the cell
  !!                             For dep<0 the recon is on the left edge of the cell
  !! @param[in]   field          Field value in the cell
  !! @param[in]   edge_left      Field value at left edge of cell
  !! @param[in]   edge_right     Field value at right edge of cell
  !! @param[in]   monotone       Monotone option to ensures no over/undershoots
  !----------------------------------------------------------------------------
  subroutine subgrid_quadratic_recon(reconstruction, dep, field, edge_left, &
                                     edge_right, monotone)

    implicit none

    real(kind=r_tran),    intent(out) :: reconstruction
    real(kind=r_tran),    intent(in)  :: dep
    real(kind=r_tran),    intent(in)  :: field
    real(kind=r_tran),    intent(in)  :: edge_left
    real(kind=r_tran),    intent(in)  :: edge_right
    integer(kind=i_def),  intent(in)  :: monotone

    real(kind=r_tran) :: cm, cc, cp
    real(kind=r_tran) :: t1, t2, t3, aa, bb

    ! Compute reconstruction weights
    if (dep >= 0.0_r_tran) then
      cp = 1.0_r_tran - 2.0_r_tran*dep + dep**2
      cc = 3.0_r_tran*dep - 2.0_r_tran*dep**2
      cm = - dep + dep**2
    else
      cp = dep + dep**2.0_r_tran
      cc = -3.0_r_tran*dep - 2.0_r_tran*dep**2
      cm = 1.0_r_tran + 2.0_r_tran*dep + dep**2
    end if

    ! Apply weights to field and field edge values
    reconstruction = cm*edge_left + cc*field + cp*edge_right

    ! Apply monotonicity if needed
    if ( monotone == horizontal_monotone_strict ) then
      ! Strict monotonicity
      t1 = (2.0_r_tran*edge_left + edge_right - 3.0_r_tran*field) &
           / (3.0_r_tran*edge_left + 3.0_r_tran*edge_right - 6.0_r_tran*field + EPS_R_TRAN)
      if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
        ! If subgrid reconstruction has extrema in the cell then revert to constant reconstruction
        reconstruction = field
      end if
    else if ( monotone == horizontal_monotone_relaxed .or. &
              monotone == horizontal_monotone_qm_pos ) then
      ! Relaxed monotonicity
      t1 = (2.0_r_tran*edge_left + edge_right - 3.0_r_tran*field) &
           / (3.0_r_tran*edge_left + 3.0_r_tran*edge_right - 6.0_r_tran*field + EPS_R_TRAN)
      if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
        ! If subgrid reconstruction has extrema in the cell then check smoothness of field
        t2 = (edge_right - field)*(field - edge_left)
        t3 = abs(field - edge_left) - abs(edge_right - field)
        if ( t2 < 0.0_r_tran ) then
          ! Revert to constant reconstruction
          reconstruction = field
        else
          ! Ensure subgrid reconstruction is bounded by edge values
          if ( t3 < 0.0_r_tran ) then
            if (dep >= 0.0_r_tran) then
              cp = 0.0_r_tran
              cc = 3.0_r_tran - 3.0_r_tran*dep + dep**2
              cm = -2.0_r_tran + 3.0_r_tran*dep - dep**2
            else
              cp = 0.0_r_tran
              cc = dep**2
              cm = 1.0_r_tran - dep**2
            end if
          else
            if (dep >= 0.0_r_tran) then
              cp = 1.0_r_tran - dep**2
              cc = dep**2.0_r_tran
              cm = 0.0_r_tran
            else
              cp = -2.0_r_tran - 3.0_r_tran*dep - dep**2
              cc = 3.0_r_tran + 3.0_r_tran*dep + dep**2
              cm = 0.0_r_tran
            end if
          end if
          reconstruction = cm*edge_left + cc*field + cp*edge_right
        end if
      end if
    else if ( monotone == horizontal_monotone_positive ) then
      ! Positive definite limiter
      aa = -4.0_r_tran*edge_left - 2.0_r_tran*edge_right + 6.0_r_tran*field
      bb = 3.0_r_tran*edge_left + 3.0_r_tran*edge_right - 6.0_r_tran*field
      ! Find stationary point of parabolic subgrid reconstruction
      t1 = -0.5_r_tran*aa / (bb + EPS_R_TRAN)
      if ( ( t1 + EPS_R_TRAN ) * ( 1.0_r_tran + EPS_R_TRAN - t1 ) > 0.0_r_tran ) then
        ! If stationary point lies within the grid cell and makes the subgrid reconstruction
        ! negative, we revert to constant reconstruction
        t2 = edge_left + aa * t1 + bb * t1 * t1
        if ( t2 .le. 0.0_r_tran ) then
          reconstruction = field
        end if
      end if
    end if

  end subroutine subgrid_quadratic_recon

  !----------------------------------------------------------------------------
  !> @brief  Returns the field value bounded by two other field values such that
  !!         min(field_one,field_two) <= field <= max(field_one,field_two)
  !!
  !! @param[inout] field       The field value to be bounded
  !! @param[in]    field_one   The first field value
  !! @param[in]    field_two   The second field value
  !----------------------------------------------------------------------------
  subroutine bound_field(field,field_one,field_two)

    implicit none

    ! Arguments
    real(kind=r_tran), intent(inout) :: field
    real(kind=r_tran), intent(in)    :: field_one
    real(kind=r_tran), intent(in)    :: field_two

    ! Min/max values
    real(kind=r_tran) :: fmin, fmax

    fmin = min(field_one,field_two)
    fmax = max(field_one,field_two)
    field = min( fmax, max(field, fmin) )

  end subroutine bound_field

  !----------------------------------------------------------------------------
  !> @brief Identify which case for Nirvana reconstruction
  !! This subroutine identifies if it leaves the reconstruction centred around
  !! cell 3 (default=case 1) away from panel edges, or shift the stencil left
  !! (case 2) or shift right (case 3)
  !! case 1 = all points on same panel
  !! case 2 = points 1 and 2 on same panel and 3 on different one (shift left)
  !! case 3 = points 2 and 3 on same panel and 1 on different one (shift right)
  !! @param[in]   ipanel  Panel IDs
  !! @param[out]  case    Identifies which cases 1, 2, or 3
  !----------------------------------------------------------------------------
  subroutine horizontal_nirvana_case(ipanel,case)

    implicit none

    integer(kind=i_def),  intent(in)   :: ipanel(1:3)
    integer(kind=i_def),  intent(out)  :: case
    integer(kind=i_def) :: test1, test2, test3

    test1 = ipanel(2) - ipanel(1)
    test2 = ipanel(2) - ipanel(3)
    test3 = test1+test2

    if ( test3 == 0_i_def ) then
      ! 3 points on same panel
      case = 1_i_def
    else
      ! Not all 3 points on same panel
      if ( test1 == 0_i_def ) then
        ! point 3 is on different panel
        case = 2_i_def
      else
        ! point 1 is on different panel
        case = 3_i_def
      end if
    end if
  end subroutine horizontal_nirvana_case

  !----------------------------------------------------------------------------------------
  !> @brief Identify which case for PPM reconstruction
  !! This subroutine identifies if it leaves the reconstruction centred around
  !! cell 3 (default=case 1) away from panel edges, or shifts the stencil left
  !! (by either 1 or 2 cells) or shifts right (by either 1 or 2 cells)
  !! case 1 = all points on same panel
  !! case 2 = points 1,2,3,4 on same panel and 5 on different one (shift left by 1 cell)
  !! case 3 = points 1,2,3 on same panel and 4, 5 on different one (shift left by 2 cells)
  !! case 4 = points 2,3,4,5 on same panel and 1 on different one (shift right by 1 cell)
  !! case 5 = points 3,4,5 on same panel and 1, 2 on different one (shift right by 2 cells)
  !! @param[in]   ipanel  Panel IDs
  !! @param[out]  case    Identifies which cases 1:5
  !----------------------------------------------------------------------------------------
  subroutine horizontal_ppm_case(ipanel,case)

    implicit none

    integer(kind=i_def),  intent(in)   :: ipanel(1:5)
    integer(kind=i_def),  intent(out)  :: case

    integer(kind=i_def) :: test(1:4)
    integer(kind=i_def) :: test1, test2

    test(1) = ipanel(3) - ipanel(2)
    test(2) = ipanel(3) - ipanel(1)
    test(3) = ipanel(3) - ipanel(4)
    test(4) = ipanel(3) - ipanel(5)
    test1 = test(1)+test(2)+test(3)+test(4)

    if ( test1 == 0_i_def ) then
      ! All points on same panel
      case = 1_i_def
    else
      ! Not all points on same panel
      test2 = test(1)+test(2)
      if (test2 == 0_i_def ) then
        ! Points 4 or 5 on different panel
          if ( test(3) == 0_i_def ) then
            ! Point 5 on different panel
             case = 2_i_def
          else
            ! Points 4 and 5 on different panel
             case = 3_i_def
          end if
      else
        ! Points 1 or 2 on different panel
         if ( test(1) == 0_i_def ) then
           ! Point 1 on different panel
            case = 4_i_def
         else
           ! Points 1 and 2 on different panel
            case = 5_i_def
         end if
      end if
    end if
  end subroutine horizontal_ppm_case

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal Nirvana to estimate the quadratic subgrid representation
  !!         of a field. The special edge treatment is used to shift the stencil left
  !!         or right depending on where the panel edge lies. As the interpolation is
  !!         third order, even with special edges, this requires a stencil of
  !!         size 5 instead of the usual 3 for Nirvana. The field dofmap is
  !!         | 1 | 2 | 3 | 4 | 5 | and this returns the estimated field value at
  !!         the edges of cell 3. The cells are assumed to be uniform in spacing.
  !!
  !! @param[in]   field            Has dof map of the form | 1 | 2 | 3 | 4 | 5 |
  !! @param[in]   spt_case         Special edges case to shift stencil
  !! @param[out]  edge_left        Field value at left edge of cell
  !! @param[out]  edge_right       Field value at right edge of cell
  !----------------------------------------------------------------------------
  subroutine nirvana_special_edge(field, spt_case, edge_left, edge_right)

    implicit none

    real(kind=r_tran),   intent(in)  :: field(1:5)
    integer(kind=i_def), intent(in)  :: spt_case
    real(kind=r_tran),   intent(out) :: edge_left
    real(kind=r_tran),   intent(out) :: edge_right

    select case ( spt_case )
    case ( 1 )
      ! All points on same panel
      edge_left  = ( -field(4) + 5.0_r_tran*field(3) + 2.0_r_tran*field(2) )/6.0_r_tran
      edge_right = ( -field(2) + 5.0_r_tran*field(3) + 2.0_r_tran*field(4) )/6.0_r_tran
    case ( 2 )
      ! Shift left
      edge_left  = ( -field(1) + 5.0_r_tran*field(2)            &
                   + 2.0_r_tran*field(3) )/6.0_r_tran
      edge_right = ( 2.0_r_tran*field(1) - 7.0_r_tran*field(2)  &
                   + 11.0_r_tran*field(3) )/6.0_r_tran
    case ( 3 )
      ! Shift right
      edge_left  = ( 11.0_r_tran*field(3) - 7.0_r_tran*field(4) &
                   + 2.0_r_tran*field(5) )/6.0_r_tran
      edge_right = ( 2.0_r_tran*field(3) + 5.0_r_tran*field(4)  &
                   - field(5) )/6.0_r_tran
    end select

  end subroutine nirvana_special_edge

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal Nirvana to estimate the quadratic subgrid representation
  !!         of a field. The special edge treatment is used to shift the stencil left
  !!         or right depending on where the panel edge lies. At panel edges the
  !!         interpolation reduces to linear (it is quadratic elsewhere) and so
  !!         this requires a stencil of size 3. The field dofmap is
  !!         | 1 | 2 | 3 | and this returns the estimated field value at
  !!         the edges of cell 2. The cells are assumed to be uniform in spacing.
  !!
  !! @param[in]   field            Has dof map of the form | 1 | 2 | 3 |
  !! @param[in]   spt_case         Special edges case to shift stencil
  !! @param[out]  edge_left        Field value at left edge of cell
  !! @param[out]  edge_right       Field value at right edge of cell
  !----------------------------------------------------------------------------
  subroutine linear_special_edge(field, spt_case, edge_left, edge_right)

    implicit none

    real(kind=r_tran),   intent(in)  :: field(1:3)
    integer(kind=i_def), intent(in)  :: spt_case
    real(kind=r_tran),   intent(out) :: edge_left
    real(kind=r_tran),   intent(out) :: edge_right

    select case ( spt_case )
    case ( 1 )
      ! All points on same panel
      edge_left  = (-field(3) + 5.0_r_tran*field(2) + 2.0_r_tran*field(1) )/6.0_r_tran
      edge_right = (-field(1) + 5.0_r_tran*field(2) + 2.0_r_tran*field(3) )/6.0_r_tran
    case ( 2 )
      ! Shift left reverting to linear
      edge_left  = ( field(2) + field(1) )/2.0_r_tran
      edge_right = ( 3.0_r_tran*field(2) - field(1) )/2.0_r_tran
    case ( 3 )
      ! Shift right reverting to linear
      edge_left  = ( 3.0_r_tran*field(2) - field(3) )/2.0_r_tran
      edge_right = ( field(2) + field(3) )/2.0_r_tran
    end select

  end subroutine linear_special_edge

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal PPM to estimate the quadratic subgrid representation
  !!         of a field. The special edge treatment is used to shift the stencil left
  !!         or right depending on where the panel edge lies. Fourth order interpolation
  !!         is used and so the stencil size increases to 7 instead of 5 for usual fourth
  !!         order interpolation. The field dofmap is | 1 | 2 | 3 | 4 | 5 | 6 | 7 | and
  !!         this returns the estimated field value at the edges of cell 4.
  !!         The cells are assumed to be uniform in spacing.
  !!
  !! @param[in]   field            Has dof map of the form | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
  !! @param[in]   spt_case         Special edges case to shift stencil
  !! @param[out]  edge_left        Field value at left edge of cell
  !! @param[out]  edge_right       Field value at right edge of cell
  !----------------------------------------------------------------------------
  subroutine fourth_order_special_edge(field, spt_case, edge_left, edge_right)

    implicit none

    real(kind=r_tran),   intent(in)  :: field(1:7)
    integer(kind=i_def), intent(in)  :: spt_case
    real(kind=r_tran),   intent(out) :: edge_left
    real(kind=r_tran),   intent(out) :: edge_right

    select case ( spt_case )
    case ( 1 )
      ! All points on same panel
      edge_left  = ( 7.0_r_tran/12.0_r_tran ) * ( field(3) + field(4) ) &
                  -( 1.0_r_tran/12.0_r_tran ) * ( field(2) + field(5) )
      edge_right = ( 7.0_r_tran/12.0_r_tran ) * ( field(4) + field(5) ) &
                  -( 1.0_r_tran/12.0_r_tran ) * ( field(3) + field(6) )
    case ( 2 )
      ! Edge between cells 5 and 6 so shift left by 1
      edge_left  = ( 7.0_r_tran/12.0_r_tran ) * ( field(3) + field(4) ) &
                  -( 1.0_r_tran/12.0_r_tran ) * ( field(2) + field(5) )
      edge_right = ( field(2) - 5.0_r_tran * field(3)                   &
                   + 13.0_r_tran * field(4) + 3.0_r_tran * field(5) )/12.0_r_tran
    case ( 3 )
      ! Shift left by 2
      edge_left  = ( field(1) - 5.0_r_tran * field(2)                   &
                   + 13.0_r_tran * field(3) + 3.0_r_tran * field(4) )/12.0_r_tran
      edge_right = ( -3.0_r_tran * field(1)  + 13.0_r_tran * field(2)   &
                   - 23.0_r_tran * field(3) + 25.0_r_tran * field(4) )/12.0_r_tran
    case ( 4 )
      ! Edge between cells 2 and 3 so shift right by 1
      edge_left  = ( 3.0_r_tran * field(3) + 13.0_r_tran * field(4)     &
                   - 5.0_r_tran * field(5) + field(6) )/12.0_r_tran
      edge_right = ( 7.0_r_tran/12.0_r_tran ) * ( field(4) + field(5) ) &
                  -( 1.0_r_tran/12.0_r_tran ) * ( field(3) + field(6) )
    case ( 5 )
      ! Shift right by 2
      edge_left  = ( -3.0_r_tran * field(7) + 13.0_r_tran * field(6)    &
                   - 23.0_r_tran * field(5) + 25.0_r_tran * field(4) )/12.0_r_tran
      edge_right = ( field(7) - 5.0_r_tran * field(6)                   &
                   + 13.0_r_tran * field(5) + 3.0_r_tran * field(4) )/12.0_r_tran
    end select

  end subroutine fourth_order_special_edge

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal PPM to estimate the quadratic subgrid representation
  !!         of a field. The special edge treatment is used to shift the stencil left
  !!         or right depending on where the panel edge lies, and reverts from fourth
  !!         order interpolation to Nirvana (third-order) interpolation at panel edges.
  !!         This does not require a larger stencil than usual fourth-order interpolation.
  !!         The field dofmap is | 1 | 2 | 3 | 4 | 5 | and this returns the estimated
  !!         field value at the edges of cell 3. The cells are assumed to be uniform in spacing.
  !!
  !! @param[in]   field            Has dof map of the form | 1 | 2 | 3 | 4 | 5 |
  !! @param[in]   spt_case         Special edges case to shift stencil
  !! @param[out]  edge_left        Field value at left edge of cell
  !! @param[out]  edge_right       Field value at right edge of cell
  !----------------------------------------------------------------------------
  subroutine fourth_nirvana_special_edge(field, spt_case, edge_left, edge_right)

    implicit none

    real(kind=r_tran),   intent(in)  :: field(1:5)
    integer(kind=i_def), intent(in)  :: spt_case
    real(kind=r_tran),   intent(out) :: edge_left
    real(kind=r_tran),   intent(out) :: edge_right

    select case ( spt_case )
    case ( 1 )
      ! All points on same panel
      edge_left  = ( 7.0_r_tran/12.0_r_tran ) * ( field(2) + field(3) ) &
                  -( 1.0_r_tran/12.0_r_tran ) * ( field(1) + field(4) )
      edge_right = ( 7.0_r_tran/12.0_r_tran ) * ( field(3) + field(4) ) &
                  -( 1.0_r_tran/12.0_r_tran ) * ( field(2) + field(5) )
    case ( 2, 4 )
      ! Edge between cells 1 and 2 or 4 and 5
      edge_left  = ( -field(4) + 5.0_r_tran*field(3) + 2.0_r_tran*field(2) )/6.0_r_tran
      edge_right = ( -field(2) + 5.0_r_tran*field(3) + 2.0_r_tran*field(4) )/6.0_r_tran
    case ( 3 )
      ! Shift left by 2
      edge_left  = ( -field(1) + 5.0_r_tran*field(2) + 2.0_r_tran*field(3) )/6.0_r_tran
      edge_right = ( 2.0_r_tran*field(1) - 7.0_r_tran*field(2)          &
                   + 11.0_r_tran*field(3) )/6.0_r_tran
    case ( 5 )
      ! Shift right by 2
      edge_left  = ( 11.0_r_tran*field(3) - 7.0_r_tran*field(4)         &
                   + 2.0_r_tran*field(5) )/6.0_r_tran
      edge_right = ( 2.0_r_tran*field(3) + 5.0_r_tran*field(4) - field(5) )/6.0_r_tran
    end select

  end subroutine fourth_nirvana_special_edge

end module subgrid_horizontal_support_mod
