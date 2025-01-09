!>
!! @par (c) Copyright
!! This software is provided under:
!!
!! The 3-Clause BSD License
!! SPDX short identifier: BSD-3-Clause
!! See https://opensource.org/licenses/BSD-3-Clause
!!
!! (c) Copyright 2016-2021 MPI-M, Joeran Maerz, Irene Stemmler;
!!     first published 2020
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright notice,
!!    this list of conditions and the following disclaimer.
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!!    this list of conditions and the following disclaimer in the documentation
!!    and/or other materials provided with the distribution.
!! 3. Neither the name of the copyright holder nor the names of its contributors
!!    may be used to endorse or promote products derived from this software
!!    without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.[7]
!!
!!
!! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
!! @file mo_m4ago_core.F90
!! @brief Module for Marine Aggregates:
!!        The Microstructure, Multiscale, Mechanistic, Marine Aggregates
!!        in the Global Ocean (M4AGO) sinking scheme
!!
!! The mo_m4ago_core module contains routines to calculate:
!!      - aggregate properties
!!      - mean sinking velocity of aggregates
!!
!! See:
!! Maerz et al. 2020: Microstructure and composition of marine aggregates
!!                    as co-determinants for vertical particulate organic
!!                    carbon transfer in the global ocean.
!!                    Biogeosciences, 17, 1765-1803,
!!                    https://doi.org/10.5194/bg-17-1765-2020
!!
!! This module is written within the project:
!! Multiscale Approach on the Role of Marine Aggregates (MARMA)
!! funded by the Max Planck Society (MPG)
!!
!! @author: joeran maerz (joeran.maerz@mpimet.mpg.de), MPI-M, HH
!! 2019, June, revised by Irene Stemmler (refactoring, cleaning), MPI-M, HH
!!
!! 2023 adopted to iHAMOCC by joeran maerz, UiB, Bergen
!! 2024 packaged as individual module (initially for iHAMOCC) by joeran maerz, UiB, Bergen
!! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
!!
!!


module mo_m4ago_core

  use mo_m4ago_kind,   only: wp
  use mo_m4ago_types,  only: aggregates, agg_environment
  use mo_m4ago_params, only: grav_acc_const,AJ1,AJ2,AJ3,BJ1,BJ2,BJ3,ONE_SIXTH,PI

  implicit none

  private

  ! Public subroutines & functions
  public :: init_m4ago_core_parameters        ! Initialization of module parameters
  public :: aggregate_properties              ! calculation of aggregate properties from primary particles information
  public :: mean_aggregate_sinking_speed
  public :: ws_Re_approx                      ! mass concentration-weighted mean sinking velocity
  public :: volweighted_agg_density           ! Aggregate volume-weighted mean aggregate density (diagnostic)
  public :: volweighted_agg_porosity          ! Aggregate Volume-weighted mean aggregate porosity (diagnostic)
  public :: conc_weighted_mean_agg_diameter   ! mass concentration-weighted mean aggregate diameter (diagnostic)

  ! Core parameters for M4AGO
  real(wp), protected :: agg_Re_crit                       ! critical diameter-based particle Reynolds number for fragmentation
  real(wp), protected :: agg_df_min,agg_df_max             ! minimum and maximum fractal dim of aggregates
  real(wp), protected :: df_slope                          ! slope of df versus stickiness mapping
  real(wp), protected :: stickiness_min,stickiness_max     ! minimum and maximum stickiness of marine aggregates

contains

  !=================================================================================================
  subroutine init_m4ago_core_parameters(Re_crit,df_min,df_max, stick_min, stick_max)
    !>
    !! Initializing fixed core parameters that don't change over run time of the model
    !!

    implicit none

    real(wp), intent(in) :: Re_crit
    real(wp), intent(in) :: df_min
    real(wp), intent(in) :: df_max
    real(wp), intent(in) :: stick_min
    real(wp), intent(in) :: stick_max

    agg_Re_crit    = Re_crit
    agg_df_min     = df_min
    agg_df_max     = df_max
    stickiness_min = stick_min
    stickiness_max = stick_max
    df_slope       = log(agg_df_min / agg_df_max)

  end subroutine init_m4ago_core_parameters

  !=================================================================================================
  subroutine aggregate_properties(aggs, agg_env)
    !>
    !! Calculate aggregate properties from individual primary particle information
    !!
    !!  - stickiness_agg: aggregate stickiness - range : 0-1 (-)
    !!  - df_agg:         aggregate fractal dimension - range: agg_df_min-agg_df_max (-)
    !!  - b_agg:          aggregate number distribution slope (-)
    !!  - av_dp:          mean primary particle diameter (m)
    !!  - av_rho_p:       mean primary particle density (kg/m3)
    !!  - dmax_agg:       maximum aggregate diameter (m)

    implicit none
    type(aggregates),intent(inout)   :: aggs
    type(agg_environment),intent(in) :: agg_env

    integer :: ipp
    real(wp)    :: stickiness_mapped
    real(wp)    :: A_total
    real(wp)    :: V_solid
    real(wp)    :: Vdpfrac

    aggs%av_dp          = 0._wp
    aggs%av_rho_p       = 0._wp
    aggs%stickiness_agg = 0._wp
    aggs%dmax_agg       = 0._wp
    aggs%b_agg          = 0._wp
    aggs%df_agg         = 0._wp
    aggs%Re_crit_agg    = 0._wp
    aggs%n_pptotal      = 0._wp

    A_total        = 0._wp
    V_solid        = 0._wp
    Vdpfrac        = 0._wp

    ! ------ calc mean aggregate stickiness
    do ipp = 1,aggs%NPrimPartTypes
       aggs%n_pptotal    = aggs%n_pptotal + aggs%n_pp(ipp)
       A_total = A_total + aggs%A_pp(ipp)
       aggs%stickiness_agg = aggs%stickiness_agg + aggs%A_pp(ipp)*aggs%stickiness_pp(ipp)
    enddo

    ! Only if particles exist, perform the calculation
    if (aggs%n_pptotal > 0._wp) then
      aggs%stickiness_agg = aggs%stickiness_agg/A_total

      ! primary particles surface weighted stickiness is mapped
      ! on range between 0 and 1
      stickiness_mapped = (aggs%stickiness_agg - stickiness_min) / (stickiness_max - stickiness_min)

      ! ------ calc fractal dimension
      ! fractal dimension of aggregates is based on that mapped stickiness
      aggs%df_agg = agg_df_max*exp(df_slope*stickiness_mapped)

      ! ------ calc the aggregate number distribution slope
      ! number distribution slope b is based on df
      ! Slope is here positive defined (as n(d)~d^-b), so *-1 of
      ! Jiang & Logan 1991: Fractal dimensions of aggregates
      ! determined from steady-state size distributions.
      ! Environ. Sci. Technol. 25, 2031-2038.
      !
      ! See also:
      ! Hunt 1980: Prediction of oceanic particle size distributions
      !            from coagulation and sedimentation mechanisms.
      !
      ! Additional assumptions made here:
      ! b in Jiang & Logan     (used for       Re <   0.1: b=1
      !                              for 0.1 < Re <  10  : b=0.871
      !                              for 10  < Re < 100  : b=0.547)
      ! is set to 0.871 as an 'average for our range of 0<Re<Re_crit'
      ! D2=min(2,df(3d)) (Meakin 1988)
      !
      ! => Formulation in Jiang & Logan 1991:
      ! slope = -0.5*(3+df+(2+df-D2)/(2-b)) reduces to:
      !
      ! careful: for df=1.5904: b_agg=2*df where w_s is undefined.
      aggs%b_agg = 0.5_wp*(3._wp + aggs%df_agg                                                     &
                        & + (2._wp + aggs%df_agg - min(2._wp, aggs%df_agg))/(2._wp - BJ2))


      ! ----- calc primary particle mean diameter and mean density
      ! primary particle mean diameter according to Bushell & Amal 1998, 2000
      ! sum(n_i) not changing
      ! previously introduced numerical factor (NUMFAC) can be pulled out and thus cancels out
      do ipp = 1,aggs%NPrimPartTypes
        aggs%av_dp   = aggs%av_dp   + aggs%n_pp(ipp)*aggs%dp_pp(ipp)**3
        Vdpfrac = Vdpfrac + aggs%n_pp(ipp)*aggs%dp_pp(ipp)**aggs%df_agg

        aggs%av_rho_p = aggs%av_rho_p + aggs%V_pp(ipp)*aggs%rho_pp(ipp)
        V_solid  = V_solid  + aggs%V_pp(ipp)
      enddo
      aggs%av_dp    = (aggs%av_dp/Vdpfrac)**(1._wp/(3._wp - aggs%df_agg))
      aggs%av_rho_p = aggs%av_rho_p/V_solid

      ! init Re_crit_agg - with a global value
      aggs%Re_crit_agg = agg_Re_crit

      ! Calculate the maximum aggregate diameter (based on critical particle Reynolds number)
      call max_agg_diam(aggs,agg_env)
    endif
  end subroutine aggregate_properties

  !=================================================================================================
  subroutine mean_aggregate_sinking_speed(aggs,agg_env)
    !
    !>
    !! mean_aggregate_sinking_speed: routine that provides the sinking velocity
    !! of aggregates to the calling BGC model
    !!

    implicit none

    type(aggregates),intent(inout)   :: aggs
    type(agg_environment),intent(in) :: agg_env

    call ws_Re_approx(aggs,agg_env)
  end subroutine

  !=================================================================================================
  subroutine ws_Re_approx(aggs,agg_env)
    !-----------------------------------------------------------------------
    !>
    !! ws_Re_approx:  distribution integrated to Lmax (Re crit dependent maximum agg size)
    !! Renolds number-dependent sinking velocity.
    !! Approximation for c_D-value taken from Jiang & Logan 1991:
    !! c_D=a*Re^-b
    !!

    implicit none

    type(aggregates),intent(inout)   :: aggs
    type(agg_environment),intent(in) :: agg_env

    if (aggs%n_pptotal > 0._wp) then
      aggs%ws_aggregates = ws_Re(aggs,agg_env)
    else
      aggs%ws_aggregates = 0._wp
    endif
  end subroutine ws_Re_approx

  !=================================================================================================
  real(wp) function get_dRe(aggs,AJ, BJ, Re,agg_env)
    !------------------------------------------------------------------------
    !>
    !! get the diameter of particles that feature a certain particle Reynolds number
    !!

    implicit none

    ! Arguments
    type(aggregates),intent(in)      :: aggs
    type(agg_environment),intent(in) :: agg_env
    real(wp), intent(in) :: AJ
    real(wp), intent(in) :: BJ
    real(wp), intent(in) :: Re

    ! Local variables

    real(wp) :: nu_vis

    nu_vis =  agg_env%mu/agg_env%rho_aq

    get_dRe = (Re*nu_vis)**((2._wp - BJ)/aggs%df_agg)/(4._wp/3._wp*(aggs%av_rho_p - agg_env%rho_aq)&
        &                                                            /agg_env%rho_aq &
        & *aggs%av_dp**(3._wp - aggs%df_agg)*grav_acc_const/(AJ*nu_vis**(BJ)))**(1._wp/aggs%df_agg)

  end function get_dRe

  !=================================================================================================
  real(wp) function get_ws_agg_integral(aggs,AJ, BJ, lower_bound, upper_bound,agg_env)
    !------------------------------------------------------------------------
    !>
    !! Calculate piecewise defined integral
    !!
    implicit none

    type(aggregates),intent(in) :: aggs
    type(agg_environment),intent(in) :: agg_env
    real(wp), intent(in) :: AJ
    real(wp), intent(in) :: BJ
    real(wp), intent(in) :: upper_bound
    real(wp), intent(in) :: lower_bound

    ! Local variables
    real(wp) :: nu_vis

    nu_vis =  agg_env%mu / agg_env%rho_aq
    get_ws_agg_integral = (4._wp/3._wp*(aggs%av_rho_p - agg_env%rho_aq)/agg_env%rho_aq             &
                   & *aggs%av_dp**(3._wp - aggs%df_agg)*grav_acc_const                             &
                   & /(AJ*nu_vis**BJ))**(1._wp/(2._wp - BJ))                                       &
                   & *(upper_bound**(1._wp - aggs%b_agg + aggs%df_agg                              &
                   & + (BJ + aggs%df_agg - 2._wp)/(2._wp - BJ))                                    &
                   & /(1._wp - aggs%b_agg + aggs%df_agg + (BJ + aggs%df_agg - 2._wp)/(2._wp - BJ)) &
                   & - lower_bound**(1._wp - aggs%b_agg + aggs%df_agg + (BJ + aggs%df_agg -2._wp)  &
                   & /(2._wp - BJ))                                                                &
                   & /(1._wp - aggs%b_agg + aggs%df_agg + (BJ + aggs%df_agg - 2._wp)/(2._wp - BJ)))

  end function get_ws_agg_integral

  !===================================================================================== ws_Re
  real(wp) function ws_Re(aggs,agg_env)
    !-----------------------------------------------------------------------
    !>
    !! ws_Re:  distribution integrated to Lmax (Re crit dependent maximum agg size)
    !! Reynolds number-dependent sinking velocity.
    !! Approximation for c_D-value taken from Jiang & Logan 1991:
    !! c_D=a*Re^-b
    !! written in such a way that we check the critical Reynolds
    !! number (in case that we extend the maximum size by shear-
    !! driven break-up).
    !!

    implicit none

    type(aggregates),intent(inout)   :: aggs
    type(agg_environment),intent(in) :: agg_env

    ! Local
    real(wp) :: d_Re01, d_Re10,wsints

    ! for Re-dependent, it should always be Re_crit_agg>10 (assumption in dmax calculation)
    ! for shear-driven break-up, check against integration bounds
    ! calc integration limits (diameters) for Re-dependent sinking:
    ! Re=0.1
    d_Re01 = get_dRe(aggs,AJ1, BJ1, 0.1_wp,agg_env)
    ! Re=10
    d_Re10 = get_dRe(aggs,AJ2, BJ2, 10._wp,agg_env)

    wsints = 0._wp

    ! assuming dmax_agg > av_dp and sum up the numerator integrals in Eq. 20
    wsints = get_ws_agg_integral(aggs,AJ1, BJ1, aggs%av_dp, min(aggs%dmax_agg,d_Re01),agg_env)
    if (aggs%dmax_agg > d_Re01) then
      wsints = wsints +get_ws_agg_integral(aggs,AJ2, BJ2, d_Re01, min(aggs%dmax_agg,d_Re10),agg_env)
    endif
    if (aggs%dmax_agg > d_Re10) then
      wsints = wsints +get_ws_agg_integral(aggs,AJ3, BJ3, d_Re10,aggs%dmax_agg,agg_env)
    endif

    ! concentration-weighted mean sinking velocity
    ws_Re = wsints                                                                                 &
            & /((aggs%dmax_agg**(1._wp + aggs%df_agg - aggs%b_agg)                                 &
            & - aggs%av_dp**(1._wp + aggs%df_agg - aggs%b_agg))                                    &
            & / (1._wp + aggs%df_agg - aggs%b_agg))  ! (m/s)

  end function ws_Re


  !=================================================================================================
  subroutine max_agg_diam(aggs,agg_env)
    !-----------------------------------------------------------------------
    !>
    !! max_agg_diam calculates the maximum aggregate diameter of the aggregate
    !! number distribution, assumes Re_crit_agg > 10
    !!
    type(aggregates),intent(inout)   :: aggs
    type(agg_environment),intent(in) :: agg_env

    ! based on analytical Jiang approximation
    aggs%dmax_agg   = max_agg_diam_white(aggs,agg_env)

  end subroutine max_agg_diam

  !================================================ maximum diameter of agg in non-stratified fluid
  real(wp)  function max_agg_diam_white(aggs,agg_env)
    !-------------------------------------------------------------------------
    !>
    !! maximum aggregate diameter in a non-stratified fluid - following the
    !! White drag approaximation by Jiang & Logan 1991, assuming Re_crit_agg > 10
    !! (otherwise AJX,BJX needs to be adjusted)
    !!

    implicit none

    type(aggregates),intent(inout)   :: aggs
    type(agg_environment),intent(in) :: agg_env
    real(wp)        :: nu_vis

    nu_vis  =  agg_env%mu / agg_env%rho_aq
    max_agg_diam_white = (aggs%Re_crit_agg*nu_vis)**((2._wp - BJ3)/aggs%df_agg)                    &
                        & /((4._wp/3._wp)*(aggs%av_rho_p - agg_env%rho_aq)/agg_env%rho_aq          &
                        & *aggs%av_dp**(3._wp - aggs%df_agg)*grav_acc_const                        &
                        & /(AJ3*nu_vis**BJ3))**(1._wp/aggs%df_agg)

  end function max_agg_diam_white

  !=================================================================================================
  !=================================================================================================
  !=================================================================================================
  ! DIAGNOSTICS

  real(wp) function volweighted_agg_density(aggs,agg_env)

    type(aggregates),intent(in)      :: aggs
    type(agg_environment),intent(in) :: agg_env

    if (aggs%n_pptotal > 0._wp) then
      ! Volume-weighted mean aggregate density
      volweighted_agg_density = (aggs%av_rho_p-agg_env%rho_aq)*aggs%av_dp**(3._wp-aggs%df_agg)     &
                            & *(4._wp-aggs%b_agg)*(aggs%dmax_agg**(1._wp+aggs%df_agg-aggs%b_agg)   &
                            &                  - aggs%av_dp**(1._wp+aggs%df_agg-aggs%b_agg))       &
                            &  / ((1._wp+aggs%df_agg-aggs%b_agg)                                   &
                            & *(aggs%dmax_agg**(4._wp-aggs%b_agg) -aggs%av_dp**(4._wp-aggs%b_agg)))&
                            & + agg_env%rho_aq
    else
      volweighted_agg_density = 0._wp
    endif

  end function volweighted_agg_density

  !=================================================================================================
  real(wp) function volweighted_agg_porosity(aggs)

    type(aggregates),intent(in) :: aggs

    if (aggs%n_pptotal > 0._wp) then
      ! Volume-weighted mean aggregate porosity
      volweighted_agg_porosity =  1._wp - ((4._wp-aggs%b_agg)*aggs%av_dp**(3._wp-aggs%df_agg)      &
                             &         *(aggs%dmax_agg**(1._wp+aggs%df_agg-aggs%b_agg)             &
                             &                     - aggs%av_dp**(1._wp+aggs%df_agg-aggs%b_agg)))  &
                             &       /((1._wp+aggs%df_agg-aggs%b_agg)                              &
                             &                                 *(aggs%dmax_agg**(4._wp-aggs%b_agg) &
                             &                                   - aggs%av_dp**(4._wp-aggs%b_agg)))
    else
      volweighted_agg_porosity = 1._wp
    endif

  end function volweighted_agg_porosity

  !=================================================================================================
  real(wp) function conc_weighted_mean_agg_diameter(aggs)

    type(aggregates),intent(in) :: aggs

    if (aggs%n_pptotal > 0._wp) then
      conc_weighted_mean_agg_diameter =  (1._wp + aggs%df_agg - aggs%b_agg)                        &
                          &             / (2._wp + aggs%df_agg - aggs%b_agg)                       &
                          & *(aggs%dmax_agg**(2._wp + aggs%df_agg - aggs%b_agg)                    &
                          &                    - aggs%av_dp**(2._wp + aggs%df_agg - aggs%b_agg))   &
                          & / (aggs%dmax_agg**(1._wp+aggs%df_agg-aggs%b_agg)                       &
                          &                    - aggs%av_dp**(1._wp + aggs%df_agg-aggs%b_agg))
    else
      conc_weighted_mean_agg_diameter = 0._wp
    endif
  end function conc_weighted_mean_agg_diameter





  !=================================================================================================
  !=================================================================================================
  !=================================================================================================
  !=================================================================================================
  ! CURRENTLY UN-USED FUNCTIONS

  real(wp)  function mass_factor(dp,df,rhop)
    !-----------------------------------------------------------------------
    !>
    !! mass_factor calculates the mass factor for the mass of a single
    !! aggregate
    !!
    implicit none

    real(wp), intent(in) :: dp
    real(wp), intent(in) :: df
    real(wp), intent(in) :: rhop

    ! mass factor
    mass_factor = ONE_SIXTH * PI * dp**(3._wp - df) * rhop

  end function mass_factor

end module mo_m4ago_core

