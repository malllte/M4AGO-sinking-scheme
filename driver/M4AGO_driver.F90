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
!! @file M4AGO_driver.F90
!! @brief Module for Marine Aggregates:
!!        The Microstructure, Multiscale, Mechanistic, Marine Aggregates
!!        in the Global Ocean (M4AGO) sinking scheme
!!
!! The M4AGO_driver module provides some test routines to improve and check
!! the M4AGO code for robustness
!!
!! 2024 packaged as individual module (initially for iHAMOCC) by joeran maerz, UiB, Bergen
!!
!! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
!!
!!

program M4AGO_driver

use mo_m4ago_kind,    only: wp
use mo_m4ago_types,   only: aggregates,agg_environment
use mo_m4ago_params,  only: rho_aq
use mo_m4ago_core,    only: ws_Re_approx,volweighted_agg_density,                                  &
                          & volweighted_agg_porosity,conc_weighted_mean_agg_diameter,              &
                          & aggregate_properties, init_m4ago_core_parameters
use driver,           only: agg_df_max,agg_df_min,agg_Re_crit,dynvis,NPrimPartTypes,stickiness_max,&
                          & stickiness_min, init_m4ago_params,prepare_primary_particles,           &
                          & print_information,init_m4ago_nml_params
!use mo_m4ago_physics, only: mol_dyn_vis

  implicit none

  ! Parameter for M4AGO core
  type(agg_environment) :: agg_env
  type(aggregates)      :: aggs

  ! For concentrations
  real(wp), allocatable :: C_det(:)
  real(wp), allocatable :: C_opal(:)
  real(wp), allocatable :: C_calc(:)
  real(wp), allocatable :: C_dust(:)

  ! Choose the test case to run
  ! Available cases:
  !    - single   : just test single concentration values
  !    - Recrit   : test single concentration values with varying Recrit values
  character(100) :: testcase = 'single'

  integer i

  call init_m4ago_nml_params
  call init_m4ago_params
  allocate(aggs%dp_pp(NPrimPartTypes))
  allocate(aggs%rho_pp(NPrimPartTypes))
  allocate(aggs%stickiness_pp(NPrimPartTypes))
  allocate(aggs%n_pp(NPrimPartTypes))
  allocate(aggs%A_pp(NPrimPartTypes))
  allocate(aggs%V_pp(NPrimPartTypes))

  select case (testcase)

    case ('single')
      print*, '========================  Running test case: single ============='
      allocate(C_det(1))
      allocate(C_opal(1))
      allocate(C_calc(1))
      allocate(C_dust(1))

      C_det   = 1e-7_wp
      C_opal  = 1e-8_wp
      C_calc  = 0._wp !1e-12_wp
      C_dust  = 0._wp !1e-11_wp
      ! Provide aggregates environment
      agg_env%rho_aq = rho_aq
      agg_env%mu     = dynvis
      ! ------ prepare primary particle information
      call prepare_primary_particles(C_det(1),C_opal(1),C_calc(1),C_dust(1),aggs,agg_env)

      ! ------ calculate aggregate properties from individual primary particle information
      call aggregate_properties(aggs, agg_env)
      ! ======== calculate the mean sinking velocity of aggregates =======
      call ws_Re_approx(aggs, agg_env)

      call print_information(aggs,agg_env)

    case ('Recrit')
      print*, '========================  Running test case: Recrit ============='
      allocate(C_det(1))
      allocate(C_opal(1))
      allocate(C_calc(1))
      allocate(C_dust(1))

      C_det   = 1e-7
      C_opal  = 1e-8
      C_calc  = 0. !1e-12
      C_dust  = 0. !1e-11

      ! Provide aggregates environment
      agg_env%rho_aq = rho_aq
      agg_env%mu     = dynvis

      ! ------ prepare primary particle information
      call prepare_primary_particles(C_det(1),C_opal(1),C_calc(1),C_dust(1),aggs,agg_env)

      do i=0,20
        ! Change the critical aggregate Reynolds number for break up
        ! and calculate the sinking velocity, etc.
        agg_Re_crit = 0.099_wp + float(i)/(20._wp/(20._wp-0.099_wp))
        call init_m4ago_core_parameters(agg_Re_crit,agg_df_min,agg_df_max,stickiness_min,stickiness_max)

        ! ------ calculate aggregate properties from individual primary particle information
        call aggregate_properties(aggs, agg_env)
        ! ======== calculate the mean sinking velocity of aggregates =======
        call ws_Re_approx(aggs, agg_env)

        print*,aggs%Re_crit_agg,aggs%dmax_agg*100._wp,aggs%ws_aggregates*86400._wp
      end do

      call print_information(aggs,agg_env)

    case default
      print*, '======================== Invalid test case ======================'

   end select

end program M4AGO_driver



! ==================================================================================================
! ==================================================================================================
! ==================================================================================================
! ==================================================================================================

module driver
  use mo_m4ago_kind,    only: wp
  use mo_m4ago_params,  only: rho_aq,ONE_SIXTH,PI
  use mo_m4ago_kind,    only: wp
  use mo_m4ago_types,   only: aggregates,agg_environment
  use mo_m4ago_params,  only: rho_aq,ONE_SIXTH,PI
  use mo_m4ago_core,    only: ws_Re_approx,volweighted_agg_density,                                &
                            & volweighted_agg_porosity,conc_weighted_mean_agg_diameter,            &
                            & aggregate_properties, init_m4ago_core_parameters

  implicit none

  private

  public :: agg_df_max,agg_df_min,agg_Re_crit,dynvis,NPrimPartTypes,stickiness_max,                &
          & stickiness_min, init_m4ago_params,prepare_primary_particles,                           &
          & print_information,init_m4ago_nml_params

  ! ------------------------------------------------------------------------------------------------
  ! biogeochemistry model-specific parameters
  ! ------------------------------------------------------------------------------------------------
  ! primary particle diameter for POM & PIM species involved in parametrized aggregation (m)
  real(wp), protected :: dp_dust ! primary particle diameter dust
  real(wp), protected :: dp_det  ! primary particle diameter detritus
  real(wp), protected :: dp_calc ! primary particle diameter calc
  real(wp), protected :: dp_opal ! primary particle diameter opal

  ! Stickiness of primary particles
  real(wp), protected :: stickiness_TEP  ! stickiness of TEP (related to opal frustules)
  real(wp), protected :: stickiness_det  ! normal detritus stickiness
  real(wp), protected :: stickiness_opal ! stickiness of opal (without TEP - just normal coating)
  real(wp), protected :: stickiness_calc ! stickiness of calc particles (coated with organics)
  real(wp), protected :: stickiness_dust ! stickiness of dust particles (coated with organics)

  real(wp), protected :: agg_df_max      ! maximum fractal dimension of aggregates (~2.5)
  real(wp), protected :: agg_df_min      ! minimum fractal dimension of aggregates (~1.2 - 1.6)
  real(wp), protected :: rho_TEP         ! density of TEP particles
  real(wp), protected :: agg_org_dens    ! organic detritus density (alternative to orgdens to avoid negative ws)
  real(wp)            :: agg_Re_crit     ! critical particle Reynolds number for fragmentation

  ! calculated model-specific parameters
  real(wp), protected :: det_mol2mass ! mol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)
  real(wp), protected :: V_dp_dust,V_dp_det,V_dp_calc,V_dp_opal   ! volumes of primary particles (L^3)
  real(wp), protected :: A_dp_dust,A_dp_det,A_dp_calc,A_dp_opal   ! surface areas of primary particles (L^2)
  real(wp), protected :: stickiness_min, stickiness_max           ! minimum and maximum stickiness of primary particles
  real(wp), protected :: rho_V_dp_dust,rho_V_dp_det,rho_V_dp_calc ! rho_V_dp_opal ! mass of primary particles (M)
  real(wp), protected :: Rm_SiP                                   ! molar mass ratio opal (SiO_2) to POM
  real(wp), protected :: thick_shell                              ! diatom frustule shell thickness (L)
  real(wp), protected :: d_frustule_inner                         ! diameter of hollow part in diatom frustule (L)
  real(wp), protected :: V_frustule_inner                         ! volume of hollow part in diatom frustule (L^3)
  real(wp), protected :: V_frustule_opal                          ! volume of opal shell material (L^3)
  real(wp), protected :: rho_V_frustule_opal                      ! mass of frustule material (M)

  real(wp), protected :: dynvis   = 0.001567_wp ! [kg/(m s)] dynamic molecular viscosity
  real(wp), protected :: calcdens = 2600._wp
  real(wp), protected :: claydens = 2600._wp
  real(wp), protected :: NUM_FAC  = 1e9_wp
  real(wp), protected :: opaldens = 2200._wp
  real(wp), protected :: opalwei  = 60._wp
  real(wp), protected :: calcwei  = 100._wp
  real(wp), protected :: ropal    = 20._wp

  integer, parameter    :: NPrimPartTypes = 4 ! Number of primary particle types generated from the biogeochemistry model

contains
  subroutine print_information(aggs,agg_env)

    implicit none

    type(aggregates),intent(in)   :: aggs
    type(agg_environment),intent(in) :: agg_env

    character(22)  :: sp = '                   ' ! just some space for pretty printing

    print*,'Primary particle types                 (#)', aggs%NPrimPartTypes
    print*,'                                             ','dust',sp,'calc',sp,'det ',sp,'opal'
    print*,'Primary particles diameter           (mum)', aggs%dp_pp*1e6_wp
    print*,'Primary particles density          (kg/m3)', aggs%rho_pp
    print*,'Number of primary particles            (#)', aggs%n_pp*NUM_FAC
    print*,'Surface area of primary particles     (m2)', aggs%A_pp
    print*,'Volume of primary particles           (m3)', aggs%V_pp/NUM_FAC
    print*,'Stickiness of primary particles        (-)', aggs%stickiness_pp
    print*,'-----------------------------------------------------------------'
    print*,'Maximum diameter                      (cm)', aggs%dmax_agg*100._wp
    print*,'Frustule stickiness                    (-)', aggs%stickiness_frustule
    print*,'Aggregate stickiness                   (-)', aggs%stickiness_agg
    print*,'Average primary particle diameter    (mum)', aggs%av_dp*1e6
    print*,'Average primary particle density   (kg/m3)', aggs%av_rho_p
    print*,'Fractal dimension                      (-)', aggs%df_agg
    print*,'Aggregate number distribution slope    (-)', aggs%b_agg
    print*,'Sinking velocity                     (m/d)', aggs%ws_aggregates*86400._wp
    print*,'-----------------------------------------------------------------'
    print*,'Conc.-weighted mean agg. diam.       (mum)', conc_weighted_mean_agg_diameter(aggs)*1e6_wp
    print*,'Volume-weighted aggregate density  (kg/m3)', volweighted_agg_density(aggs,agg_env)
    print*,'Volume-weighted aggregate porosity     (-)', volweighted_agg_porosity(aggs)
    print*,'-----------------------------------------------------------------'
    print*,'-----------------------------------------------------------------'
  end subroutine print_information

  subroutine init_m4ago_nml_params
    !>
    !! Initialization of namelist parameters
    !!

    implicit none

    ! Primary particle diameters
    dp_dust = 2.e-6_wp      ! [m] following the classical HAMOCC parametrization
    dp_det  = 4.e-6_wp      ! [m] not well defined
    dp_calc = 3.e-6_wp      ! [m] following Henderiks 2008, Henderiks & Pagani 2008
    dp_opal = 20.e-6_wp     ! [m] mean frustule diameter of diatoms

    ! Stickiness values - note that their relative values to each other matter!
    stickiness_TEP    = 0.19_wp ! [-] range 0-1
    stickiness_det    = 0.1_wp  ! [-] range 0-1
    stickiness_opal   = 0.08_wp ! [-] range 0-1
    stickiness_calc   = 0.09_wp ! [-] range 0-1
    stickiness_dust   = 0.07_wp ! [-] range 0-1

    ! Minimum and maximum aggregate fractal dimension
    agg_df_min        = 1.6_wp ! [-]
    agg_df_max        = 2.4_wp ! [-]

    ! Density of primary particles
    rho_TEP           = 800._wp  ! [kg/m^3] 700.-840. kg/m^3 Azetsu-Scott & Passow 2004
    agg_org_dens      = 1100._wp ! [kg/m^3] detritus density - don't use orgdens to avoid negative ws

    ! Critical particle Reynolds number (based on diameter) for limiting nr-distribution
    agg_Re_crit       = 20._wp ! [-]

  end subroutine init_m4ago_nml_params

  subroutine init_m4ago_params
    !>
    !! Initilization of parameters
    !!

    implicit none
    det_mol2mass   = 3166._wp  ! mol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)

    ! Volume of an individual primary particle*NUMFAC
    V_dp_dust = ONE_SIXTH*PI*dp_dust**3*NUM_FAC
    V_dp_det  = ONE_SIXTH*PI*dp_det**3 *NUM_FAC
    V_dp_calc = ONE_SIXTH*PI*dp_calc**3*NUM_FAC
    V_dp_opal = ONE_SIXTH*PI*dp_opal**3*NUM_FAC

    ! Surface area of an individual primary particle*NUMFAC
    A_dp_dust = PI*dp_dust**2*NUM_FAC
    A_dp_det  = PI*dp_det**2 *NUM_FAC
    A_dp_calc = PI*dp_calc**2*NUM_FAC
    A_dp_opal = PI*dp_opal**2*NUM_FAC

    ! Mass of an individual primary particle*NUMFAC
    rho_V_dp_dust = V_dp_dust*claydens
    rho_V_dp_det  = V_dp_det*agg_org_dens
    rho_V_dp_calc = V_dp_calc*calcdens

    Rm_SiP              = ropal*opalwei/det_mol2mass
    ! shell thickness
    thick_shell         = 0.5_wp*dp_opal*(1._wp - (opaldens/(Rm_SiP*agg_org_dens+opaldens))**(1._wp/3._wp))
    d_frustule_inner    = dp_opal - 2._wp*thick_shell
    ! volume of hollow part of frustule
    V_frustule_inner    = ONE_SIXTH* PI*d_frustule_inner**3*NUM_FAC
    ! volume of opal part of frustule
    V_frustule_opal     = ONE_SIXTH*PI*(dp_opal**3 - d_frustule_inner**3)*NUM_FAC
    rho_V_frustule_opal = V_frustule_opal*opaldens

    ! Minimum and maximum reachable stickiness
    stickiness_min      = min(stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust)
    stickiness_max      = max(stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust)

    ! Init core M4AGO parameters
    call init_m4ago_core_parameters(agg_Re_crit,agg_df_min,agg_df_max,stickiness_min,stickiness_max)

  end subroutine init_m4ago_params

  subroutine prepare_primary_particles(C_det,C_opal,C_calc,C_dust,aggs,agg_env)
    !-----------------------------------------------------------------------
    !>
    !! prepare_primary_particles
    !! calculates/provides fields with:
    !!  - primary particle diameter
    !!  - primary particle density
    !!  - number of primary particles
    !!  - surface area of primary particles
    !!  - volume of primary particles
    !!  - stickiness of the primary particles
    !! based on the driving ocean biogeochmistry model tracer field

    implicit none

    real(wp), intent(in)  :: C_det                  !< detritus concentration kmol P/m3
    real(wp), intent(in)  :: C_opal                 !< opal concentration kmol Si/m3
    real(wp), intent(in)  :: C_calc                 !< CaCO3 concentration kmol Ca/m3
    real(wp), intent(in)  :: C_dust                 !< dust kg/m3
    type(aggregates),intent(inout)   :: aggs
    type(agg_environment),intent(in) :: agg_env

    real(wp) :: n_det,n_opal,n_calc,n_dust         ! total primary particle number (#)
    real(wp) :: A_dust,A_det,A_calc,A_opal         ! total surface area of primary particles per unit volume (L^2/L^3)
    real(wp) :: V_det,V_opal,V_calc,V_dust,V_solid ! total volume of primary particles in a unit volume (L^3/L^3)

    real(wp) :: cell_det_mass                      ! mass of detritus material in diatoms
    real(wp) :: cell_pot_det_mass                  ! potential (max) mass detritus material in diatoms
    real(wp) :: free_detritus                      ! freely available detritus mass outside the frustule
    real(wp) :: V_POM_cell                         ! volume of POM in frustule
    real(wp) :: V_aq                               ! volume of water space in frustule
    real(wp) :: rho_frustule                       ! density of diatom frustule incl. opal, detritus and water
    real(wp) :: rho_diatom                         ! density of either hollow frustule or with additions of detritus and water
    real(wp) :: stickiness_frustule                ! stickiness of the diatom frustile as primary particle

    n_det   = 0._wp ! number of primary particles in a unit volume
    n_opal  = 0._wp
    n_dust  = 0._wp
    n_calc  = 0._wp

    V_det   = 0._wp ! total volume of primary particles in a unit volume
    V_opal  = 0._wp
    V_calc  = 0._wp
    V_dust  = 0._wp
    V_solid = 0._wp

    ! n_det are detritus primary particle that are
    ! NOT linked to any diatom frustule
    ! n_opal are number of frustule-like primary particles possessing
    ! a density i) different from pure opal ii) due to a mixture of
    ! opal frustule, detritus inside the frustule and potentially water
    ! inside the frustule

    ! describing diatom frustule as hollow sphere
    ! that is completely or partially filled with detritus
    ! and water
    free_detritus     = 0._wp
    rho_diatom        = 0._wp
    cell_det_mass     = 0._wp
    cell_pot_det_mass = 0._wp
    V_POM_cell        = 0._wp
    V_aq              = 0._wp
    rho_frustule      = 0._wp
    stickiness_frustule = 0._wp

    ! number of opal frustules (/NUM_FAC)
    n_opal = C_opal*opalwei/rho_V_frustule_opal

    ! maximum mass of detritus inside a frustule
    cell_pot_det_mass = n_opal*V_frustule_inner*agg_org_dens

    ! detritus mass inside frustules
    cell_det_mass = max(0._wp,min(cell_pot_det_mass,C_det*det_mol2mass))

    if (n_opal > 0._wp) then
      ! volume of detritus component in diatom cell
      V_POM_cell = (cell_det_mass/n_opal)/agg_org_dens

      ! if not enough detritus is available, water is added
      V_aq = V_frustule_inner -  V_POM_cell

      ! density of the diatom frustules incl. opal, detritus and water
      rho_frustule = (rho_V_frustule_opal + cell_det_mass/n_opal + V_aq*agg_env%rho_aq)/V_dp_opal
      rho_diatom = (rho_frustule + cell_det_mass/cell_pot_det_mass*rho_TEP)                        &
                   /(1._wp + cell_det_mass/cell_pot_det_mass)

      ! calc frustule stickiness
      stickiness_frustule = cell_det_mass/(cell_pot_det_mass)*stickiness_TEP                       &
                               & + (1._wp - cell_det_mass/(cell_pot_det_mass))                     &
                               &   *stickiness_opal
    endif

    ! mass of extra cellular detritus particles
    free_detritus = C_det*det_mol2mass  - cell_det_mass

    ! number of primary particles
    n_det  = free_detritus/rho_V_dp_det  ! includes NUM_FAC
    n_calc = C_calc*calcwei/rho_V_dp_calc
    n_dust = C_dust/rho_V_dp_dust        ! dust is in kg/m3

    ! calc total areas
    A_det   = n_det*A_dp_det
    A_opal  = n_opal*A_dp_opal
    A_calc  = n_calc*A_dp_calc
    A_dust  = n_dust*A_dp_dust

    ! total volume of primary particles
    V_det   = n_det*V_dp_det*NUM_FAC
    V_opal  = n_opal*V_dp_opal*NUM_FAC
    V_calc  = n_calc*V_dp_calc*NUM_FAC
    V_dust  = n_dust*V_dp_dust*NUM_FAC


    ! IMPORTANT: the order requires to be the same for all information
    ! NUMFAC cancels out in the subsequent calculations
    aggs%NPrimPartTypes = NPrimPartTypes
    aggs%dp_pp          = (/ dp_dust,  dp_calc,  dp_det,       dp_opal /)
    aggs%rho_pp         = (/ claydens, calcdens, agg_org_dens, rho_diatom /)
    aggs%n_pp           = (/ n_dust,   n_calc,   n_det,        n_opal  /)
    aggs%A_pp           = (/ A_dust,   A_calc,   A_det,        A_opal /)
    aggs%V_pp           = (/ V_dust,   V_calc,   V_det,        V_opal /)
    aggs%stickiness_pp  = (/ stickiness_dust, stickiness_calc, stickiness_det, stickiness_frustule /)
    aggs%stickiness_frustule = stickiness_frustule


  end subroutine prepare_primary_particles

end module driver
