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

  implicit none

  private

  ! Public subroutines & functions
  public :: init_m4ago_core_parameters        ! Initialization of module parameters
  public :: aggregate_properties              ! calculation of aggregate properties from primary particles information
  public :: ws_Re_approx                      ! mass concentration-weighted mean sinking velocity
  public :: volweighted_agg_density           ! Aggregate volume-weighted mean aggregate density (diagnostic)
  public :: volweighted_agg_porosity          ! Aggregate Volume-weighted mean aggregate porosity (diagnostic)
  public :: conc_weighted_mean_agg_diameter   ! mass concentration-weighted mean aggregate diameter (diagnostic)

  ! Public values
  public :: av_dp,av_rho_p,df_agg,b_agg,Lmax_agg,rho_aq,stickiness_agg,ONE_SIXTH,PI,ws_aggregates

  real, protected :: av_dp                             ! mean primary particle diameter
  real, protected :: av_rho_p                          ! mean primary particle density
  real, protected :: df_agg                            ! fractal dimesnion of aggregates
  real, protected :: b_agg                             ! aggregate number distribution slope
  real, protected :: Lmax_agg                          ! maximum diamater of aggregates
  real, protected :: agg_Re_crit                       ! critical diameter-based particle Reynolds number for fragmentation
  real, protected :: agg_df_min,agg_df_max             ! minimum and maximum fractal dim of aggregates
  real, protected :: df_slope                          ! slope of df versus stickiness mapping
  real, protected :: stickiness_min,stickiness_max     ! minimum and maximum stickiness of marine aggregates
  real, protected :: stickiness_agg                    ! mean aggregate stickiness
  real, protected :: ws_aggregates                     ! mean mass concentration-weighted sinking velocity of aggregates
  !$OMP THREADPRIVATE(av_dp,av_rho_p,df_agg,b_agg,Lmax_agg,rho_aq,stickiness_agg,ONE_SIXTH,PI,ws_aggregates)


  real, parameter :: rho_aq         = 1025.            ! water reference density  (1025 kg/m^3)
  real, parameter :: grav_acc_const = 9.81             ! gravitational acceleration constant

  ! constants for the drag coefficient CD according to Ji & Logan 1991
  real, parameter :: AJ1 = 24.00
  real, parameter :: AJ2 = 29.03
  real, parameter :: AJ3 = 14.15
  real, parameter :: BJ1 = 1.0
  real, parameter :: BJ2 = 0.871
  real, parameter :: BJ3 = 0.547


  ! Helping parameters
  real, parameter :: EPS_ONE   = EPSILON(1.0D0)
  real, parameter :: ONE_SIXTH = 1./6.
  real, parameter :: PI        = 3.141592654

contains

  !=================================================================================================
  subroutine init_m4ago_core_parameters(Re_crit,df_min,df_max, stick_min, stick_max)
    !>
    !! Initializing fixed core parameters that don't change over run time of the model
    !!

    implicit none

    real, intent(in) :: Re_crit
    real, intent(in) :: df_min
    real, intent(in) :: df_max
    real, intent(in) :: stick_min
    real, intent(in) :: stick_max

    agg_Re_crit    = Re_crit
    agg_df_min     = df_min
    agg_df_max     = df_max
    stickiness_min = stick_min
    stickiness_max = stick_max
    df_slope       = log(agg_df_min / agg_df_max)

  end subroutine init_m4ago_core_parameters

  !=================================================================================================
  subroutine aggregate_properties(Ntypes,dp_pp,rho_pp,n_pp,A_pp,V_pp,stickiness_pp, mu)
    !>
    !! Calculate aggregate properties from individual primary particle information
    !!
    !!  - stickiness_agg: aggregate stickiness - range : 0-1 (-)
    !!  - df_agg:         aggregate fractal dimension - range: agg_df_min-agg_df_max (-)
    !!  - b_agg:          aggregate number distribution slope (-)
    !!  - av_dp:          mean primary particle diameter (m)
    !!  - av_rho_p:       mean primary particle density (kg/m3)
    !!  - Lmax_agg:       maximum aggregate diameter (m)

    implicit none

    integer,                 intent(in) :: Ntypes        ! Number of individual primary particle types
    real, dimension(Ntypes), intent(in) :: dp_pp         ! Primary particle diameters (m)
    real, dimension(Ntypes), intent(in) :: rho_pp        ! Primary particle densities (kg/m3)
    real, dimension(Ntypes), intent(in) :: n_pp          ! Number of individual primary particles per unit volume (#/m3)
    real, dimension(Ntypes), intent(in) :: A_pp          ! Total surface area of individual primary particles per unit volume (m2/m3)
    real, dimension(Ntypes), intent(in) :: V_pp          ! Total volume of individual primary particles per unit volume (m3/m3)
    real, dimension(Ntypes), intent(in) :: stickiness_pp ! Stickiness of individual primary particles - value range: 0-1 (-)
    real,                    intent(in) :: mu            ! molecular dynamic viscosity of sea water (kg/(m*s))

    integer :: ipp
    real    :: stickiness_mapped
    real    :: A_total
    real    :: V_solid
    real    :: Vdpfrac

    av_dp          = 0.
    av_rho_p       = 0.
    stickiness_agg = 0.
    A_total        = 0.
    V_solid        = 0.
    Vdpfrac        = 0.

    ! ------ calc mean aggregate stickiness
    do ipp = 1,Ntypes
       A_total = A_total + A_pp(ipp)
       stickiness_agg = stickiness_agg + A_pp(ipp)*stickiness_pp(ipp)
    enddo
    stickiness_agg = stickiness_agg/(A_total+EPS_ONE)

    ! primary particles surface weighted stickiness is mapped
    ! on range between 0 and 1
    stickiness_mapped = (stickiness_agg - stickiness_min) / (stickiness_max - stickiness_min)

    ! ------ calc fractal dimension
    ! fractal dimension of aggregates is based on that mapped stickiness
    df_agg = agg_df_max*exp(df_slope*stickiness_mapped)

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
    b_agg = 0.5*(3. + df_agg + (2. + df_agg - min(2., df_agg))/(2. - BJ2))


    ! ----- calc primary particle mean diameter and mean density
    ! primary particle mean diameter according to Bushell & Amal 1998, 2000
    ! sum(n_i) not changing - can be pulled out and thus cancels out
    do ipp = 1,Ntypes
      av_dp   = av_dp   + n_pp(ipp)*dp_pp(ipp)**3.
      Vdpfrac = Vdpfrac + n_pp(ipp)*dp_pp(ipp)**df_agg

      av_rho_p = av_rho_p + V_pp(ipp)*rho_pp(ipp)
      V_solid  = V_solid  + V_pp(ipp)
    enddo
    av_dp    = (av_dp/Vdpfrac)**(1./(3. - df_agg))
    av_rho_p = av_rho_p/V_solid
    !    av_dp    = (av_dp/(Vdpfrac+EPS_ONE))**(1./(3. - df_agg))
    !    av_rho_p = av_rho_p/(V_solid+EPS_ONE)

    ! Calculate the maximum aggregate diameter (based on critical particle Reynolds number)
    call max_agg_diam(mu)

  end subroutine aggregate_properties

  !=================================================================================================
  subroutine ws_Re_approx(mu)
    !-----------------------------------------------------------------------
    !>
    !! ws_Re_approx:  distribution integrated to Lmax (Re crit dependent maximum agg size)
    !! Renolds number-dependent sinking velocity.
    !! Approximation for c_D-value taken from Jiang & Logan 1991:
    !! c_D=a*Re^-b
    !!

    implicit none

    real, intent(in) :: mu

    ws_aggregates = ws_Re(Lmax_agg,mu)

  end subroutine ws_Re_approx

  !=================================================================================================
  real function get_dRe(AJ, BJ, Re,mu)
    !------------------------------------------------------------------------
    !>
    !! get the diameter of particles that feature a certain particle Reynolds number
    !!

    implicit none
    ! Arguments
    real, intent(in) :: AJ
    real, intent(in) :: BJ
    real, intent(in) :: Re
    real, intent(in) :: mu

    ! Local variables

    real :: nu_vis

    nu_vis =  mu/rho_aq

    get_dRe = (Re*nu_vis)**((2. - BJ)/df_agg)/(4./3.*(av_rho_p - rho_aq)/rho_aq                    &
           *av_dp**(3. - df_agg)*grav_acc_const/(AJ*nu_vis**(BJ)))**(1./df_agg)

  end function get_dRe

  !=================================================================================================
  real function get_ws_agg_integral(AJ, BJ, lower_bound, upper_bound,mu)
    !------------------------------------------------------------------------
    !>
    !! Calculate piecewise defined integral
    !!
    implicit none

    real, intent(in) :: AJ
    real, intent(in) :: BJ
    real, intent(in) :: upper_bound
    real, intent(in) :: lower_bound
    real, intent(in) :: mu

    ! Local variables
    real :: nu_vis

    nu_vis =  mu/rho_aq
    get_ws_agg_integral = (4./3.*(av_rho_p - rho_aq)/rho_aq                                        &
                        & *av_dp**(3. - df_agg)*grav_acc_const                                     &
                        & /(AJ*nu_vis**BJ))**(1./(2. - BJ))                                        &
                        & *(upper_bound**(1. - b_agg + df_agg                                      &
                        & + (BJ + df_agg - 2.)/(2. - BJ))                                          &
                        & /(1. - b_agg + df_agg + (BJ + df_agg - 2.)/(2. - BJ))                    &
                        & - lower_bound**(1. - b_agg + df_agg + (BJ + df_agg -2.)                  &
                        & /(2. - BJ))                                                              &
                        & /(1. - b_agg + df_agg + (BJ + df_agg - 2.)/(2. - BJ)))

  end function get_ws_agg_integral

  !===================================================================================== ws_Re
  real function ws_Re(dmax_agg,mu)
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

    real, intent(in) :: dmax_agg
    real, intent(in) :: mu

    ! Local
    real :: d_Re01, d_Re10, d_low, ws_agg_ints

    ! for Re-dependent, it should always be agg_Re_crit>10
    ! for shear-driven break-up, check against integration bounds
    ! calc integration limits for Re-dependent sinking:
    ! Re=0.1
    d_Re01 = get_dRe(AJ1, BJ1, 0.1,mu)
    ! Re=10
    d_Re10 = get_dRe(AJ2, BJ2, 10.,mu)
    d_low = av_dp

    ws_agg_ints = 0.
    if(dmax_agg >= d_Re01)then ! Re > 0.1
                                       ! - collect full range up to
                                       ! 0.1, (dp->d_Re1) and set lower bound to
                                       ! Re=0.1 val
                                       ! aj=AJ1, bj=1
        ws_agg_ints = get_ws_agg_integral(AJ1, BJ1, av_dp, d_Re01,mu)
        d_low = d_Re01
    endif

    if(dmax_agg >= d_Re10)then ! Re > 10
                                         ! - collect full range Re=0.1-10 (d_Re1-> d_Re2)
                                         ! and set lower bound to
                                         ! Re=10 val
                                         ! aj=AJ2, bj=0.871
        ws_agg_ints = ws_agg_ints  + get_ws_agg_integral(AJ2, BJ2, d_Re01, d_Re10, mu)
        d_low = d_Re10
    endif

    if(d_low < d_Re01)then ! Re<0.1 and Lmax < d_Re1
        ws_agg_ints = get_ws_agg_integral(AJ1, BJ1, av_dp, dmax_agg, mu)
    else ! Re > 10, aj=AJ3, bj=BJ3
        ws_agg_ints = ws_agg_ints + get_ws_agg_integral(AJ3, BJ3, d_low, dmax_agg,mu)
    endif

    ! concentration-weighted mean sinking velocity
    ws_Re = (ws_agg_ints                                                                           &
            & /((dmax_agg**(1. + df_agg - b_agg)                                                   &
            & - av_dp**(1. + df_agg - b_agg))                                                      &
            & / (1. + df_agg - b_agg)))  ! (m/s)

  end function ws_Re


  !=================================================================================================
  subroutine max_agg_diam(mu)
    !-----------------------------------------------------------------------
    !>
    !! max_agg_diam calculates the maximum aggregate diameter of the aggregate
    !! number distribution, assumes agg_Re_crit > 10
    !!
    real, intent(in)  :: mu                  !< molecular dynamic viscosity

    ! base on analytical Jiang approximation
    Lmax_agg   = max_agg_diam_white(mu)
  end subroutine max_agg_diam

  !================================================ maximum diameter of agg in non-stratified fluid
  real  function max_agg_diam_white(mu)
    !-------------------------------------------------------------------------
    !>
    !! maximum aggregate diameter in a non-stratified fluid - following the
    !! White drag approaximation by Jiang & Logan 1991, assuming agg_Re_crit > 10
    !! (otherwise AJX,BJX needs to be adjusted)
    !!

    implicit none

    real,intent(in) :: mu
    real        :: nu_vis

    nu_vis  =  mu/rho_aq
    max_agg_diam_white = (agg_Re_crit*nu_vis)**((2. - BJ3)/df_agg)                                 &
                        & /((4./3.)*(av_rho_p - rho_aq)/rho_aq                                     &
                        & *av_dp**(3. - df_agg)*grav_acc_const                                     &
                        & /(AJ3*nu_vis**BJ3))**(1./df_agg)

  end function max_agg_diam_white

  !=================================================================================================
  !=================================================================================================
  !=================================================================================================
  ! DIAGNOSTICS

  real function volweighted_agg_density()

    ! Volume-weighted mean aggregate density
    volweighted_agg_density = (av_rho_p-rho_aq)*av_dp**(3.-df_agg)                                 &
                            & *(4.-b_agg)*(Lmax_agg**(1.+df_agg-b_agg) - av_dp**(1.+df_agg-b_agg)) &
                            &  / ((1.+df_agg-b_agg)                                                &
                            & *(Lmax_agg**(4.-b_agg) - av_dp**(4.-b_agg))) + rho_aq

  end function volweighted_agg_density

  !=================================================================================================
  real function volweighted_agg_porosity()

    ! Volume-weighted mean aggregate porosity
    volweighted_agg_porosity =  1. - ((4.-b_agg)*av_dp**(3.-df_agg)                                &
                             &         *(Lmax_agg**(1.+df_agg-b_agg) - av_dp**(1.+df_agg-b_agg)))  &
                             &       /((1.+df_agg-b_agg)*(Lmax_agg**(4.-b_agg) - av_dp**(4.-b_agg)))

  end function volweighted_agg_porosity

  !=================================================================================================
  real function conc_weighted_mean_agg_diameter()
    conc_weighted_mean_agg_diameter =  (1. + df_agg - b_agg) / (2. + df_agg - b_agg)               &
                          & *(Lmax_agg**(2. + df_agg - b_agg) - av_dp**(2. + df_agg - b_agg))      &
                          & / (Lmax_agg**(1.+df_agg-b_agg) - av_dp**(1. + df_agg-b_agg))

  end function conc_weighted_mean_agg_diameter





  !=================================================================================================
  !=================================================================================================
  !=================================================================================================
  !=================================================================================================
  ! CURRENTLY UN-USED FUNCTIONS

  real  function mass_factor(dp,df,rhop)
    !-----------------------------------------------------------------------
    !>
    !! mass_factor calculates the mass factor for the mass of a single
    !! aggregate
    !!
    implicit none

    real, intent(in) :: dp
    real, intent(in) :: df
    real, intent(in) :: rhop

    ! mass factor
    mass_factor = ONE_SIXTH * PI * dp**(3. - df) * rhop

  end function mass_factor


  !=================================================================================================
  real function rho_agg(d,rhop,dp,df,rho)
    !-----------------------------------------------------------------------
    !>
    !! rho_agg provides the aggregate density
    !!

    implicit none

    real, intent(in) :: d
    real, intent(in) :: rhop
    real, intent(in) :: dp
    real, intent(in) :: df
    real, intent(in) :: rho

    rho_agg =  (rhop - rho)*(dp/d)**(3. - df) + rho

  end function rho_agg

  !=================================================================================================
  real function Re_fun(ws,d,mu,rho)
    !-----------------------------------------------------------------------
    !>
    !! Particle Reynolds number for settling particles (based on diameter)
    !!

    implicit none

    real,intent(in) :: ws,d,mu,rho

    Re_fun = abs(ws*d*rho/mu)

  end function Re_fun

end module mo_m4ago_core

