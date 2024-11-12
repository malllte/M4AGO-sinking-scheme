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
!! @file mo_m4ago_climberx.F90
!! @brief Module for Marine Aggregates:
!!        The Microstructure, Multiscale, Mechanistic, Marine Aggregates
!!        in the Global Ocean (M4AGO) sinking scheme
!!
!! The mo_m4ago_climberx module contains routines to calculate:
!!      - primary particles from climberx tracers
!!      - diagnostics and enables to return a 3D sinking velocity field
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
!! 2024 adopted to CLIMBER-X (based on mo_ihammoc4m4ago) by Malte Heinemann, Kiel University
!!
!! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
!!
!!

module mo_m4ago_climberx

  ! CLIMBER-X-bgc-model-specific parameters:
  use control,          only: wp
  use bgc_grid,         only: kpke => kt
  use bgc_params,       only: dtb, dtbgc
  use bgc_params,       only: calcdens, claydens, opaldens, calcwei, opalwei, ropal
  !use bgc_def,          only: carb_t
  !use bgc_params,       only: iopal, ifdust, icalc, idet

  ! M4AGO subroutines, functions, variables and parameters:
  use mo_m4ago_core,    only: rho_aq,ONE_SIXTH,PI,aggregates,agg_environment,                      &
                            & ws_Re_approx,volweighted_agg_density,                                &
                            & volweighted_agg_porosity,conc_weighted_mean_agg_diameter,            &
                            & aggregate_properties, init_m4ago_core_parameters

  use mo_m4ago_physics, only: mol_dyn_vis

  implicit none

  private

  ! Public subroutines / called from subroutine ocprod in CLIMBER-X
  public :: climberx_mean_aggregate_sinking_speed, init_m4ago_nml_params, init_m4ago_params


  ! Public fields and parameters
  ! public :: ws_agg

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
  real(wp), protected :: agg_Re_crit     ! critical particle Reynolds number for fragmentation

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

  ! Parameter for M4AGO core
  integer, parameter :: NPrimPartTypes = 4 ! Number of primary particle types generated from the biogeochemistry model

  ! Fields
  !real(wp),allocatable :: ws_agg(:,:,:)       ! mass concentration-weighted aggregate mean sinking velocity
  !real(wp),allocatable :: dyn_vis(:,:,:)      ! molecular dynamic viscosity
  !real(wp),allocatable :: m4ago_ppo(:,:,:)    ! pressure


  ! Internally used parameters and values
  real(wp), parameter :: NUM_FAC = 1.e9_wp             ! factor to avoid numerical precision problems
  real(wp), parameter :: EPS_ONE = EPSILON(1._wp)

contains

  !===================================================================================== m4ago_init_params
  subroutine init_m4ago_nml_params
    !>
    !! Initialization of namelist parameters
    !!

    implicit none

    ! Primary particle sizes
    dp_dust = 2.e-6_wp      ! following the classical HAMOCC parametrization
    dp_det  = 4.e-6_wp      ! not well defined
    dp_calc = 3.e-6_wp      ! following Henderiks 2008, Henderiks & Pagani 2008
    dp_opal = 20.e-6_wp     ! mean frustule diameter of diatoms

    ! Stickiness values - note that their relative values to each other matter!
    stickiness_TEP    = 0.19_wp
    stickiness_det    = 0.1_wp
    stickiness_opal   = 0.08_wp
    stickiness_calc   = 0.09_wp
    stickiness_dust   = 0.07_wp

    ! Minimum and maximum aggregate fractal dimension
    agg_df_min        = 1.6_wp
    agg_df_max        = 2.4_wp

    ! Density of primary particles
    rho_TEP           = 800._wp ! 700.-840. kg/m^3 Azetsu-Scott & Passow 2004
    agg_org_dens      = 1100._wp ! detritus density - don't use orgdens to avoid negative ws

    ! Critical particle Reynolds number for limiting nr-distribution
    agg_Re_crit       = 20._wp

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

  !===================================================================================== mean_agg_ws
  !FIXME: MOVE BACK TO TRA?
  !subroutine climberx_mean_aggregate_sinking_speed(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, ppao, prho)
  !FIXME: structure tra with intent(INOUT) may be causing issues / unexpected results
  !       hence, for now not passing complete tra, but individual variables with set intent IN or OUT
  !subroutine climberx_mean_aggregate_sinking_speed(tra, kbo, pddpo, ptho, psao, layer_depth, prho)
  subroutine climberx_mean_aggregate_sinking_speed(CC_det, CC_opal, CC_calc, CC_dust, &         ! IN
  &                                                kbo, pddpo, ptho, psao, layer_depth, prho, & ! IN
  &                                                wsagg, visco, xLmaxagg, xavdp, xavrhop, xdfagg, xbagg, xavdC, xavrhofV, xavporV, & ! OUT
  &                                                xstickagg, xstickfrust)                                                            ! OUT

    !-----------------------------------------------------------------------
    !>
    !! calculates the mass concentration-weighted mean sinking velocity of marine
    !! aggregates
    !!

    implicit none

    ! TYPE(carb_t), INTENT(inout)  :: tra
    ! integer, intent(in)  :: kpke                !< 3rd (vertical) real of model grid.
    real(wp), intent(in)  :: CC_det(kpke)         !< detritus concentration []
    real(wp), intent(in)  :: CC_opal(kpke)        !< opal concentration []
    real(wp), intent(in)  :: CC_calc(kpke)        !< calcite concentration []
    real(wp), intent(in)  :: CC_dust(kpke)        !< free dust concentration []
    integer,  intent(in)  :: kbo                      !< 3rd (vertical) real of model grid.
    real(wp), intent(in)  :: pddpo(kpke)          !< size of scalar grid cell (3rd dimension) [m]
    real(wp), intent(in)  :: ptho(kpke)           !< potential temperature [deg C]
    real(wp), intent(in)  :: psao(kpke)           !< salinity [psu.].
    real(wp), intent(in)  :: layer_depth(kpke)    !< layer depth [m] 
    real(wp), intent(in)  :: prho(kpke)           !< water density [g/cm3]

    real(wp), intent(out) :: wsagg(kpke)
    real(wp), intent(out) :: visco(kpke)
    real(wp), intent(out) :: xLmaxagg(kpke)
    real(wp), intent(out) :: xavdp(kpke)
    real(wp), intent(out) :: xavrhop(kpke)
    real(wp), intent(out) :: xdfagg(kpke)
    real(wp), intent(out) :: xbagg(kpke)
    real(wp), intent(out) :: xavdC(kpke)
    real(wp), intent(out) :: xavrhofV(kpke)
    real(wp), intent(out) :: xavporV(kpke)
    real(wp), intent(out) :: xstickagg(kpke)
    real(wp), intent(out) :: xstickfrust(kpke)

    integer  :: k
    type(agg_environment) :: agg_env
    type(aggregates)      :: aggs
    
    allocate(aggs%dp_pp(NPrimPartTypes))
    allocate(aggs%rho_pp(NPrimPartTypes))
    allocate(aggs%stickiness_pp(NPrimPartTypes))
    allocate(aggs%n_pp(NPrimPartTypes))
    allocate(aggs%A_pp(NPrimPartTypes))
    allocate(aggs%V_pp(NPrimPartTypes))


    !$OMP PARALLEL DO PRIVATE(k,aggs,agg_env)
    do k = 1,kbo
        visco(k)=0.
        wsagg(k)=0.


        ! ------ provide aggregates environment
        visco(k) = mol_dyn_vis(layer_depth(k),ptho(k),psao(k)) ! molecular dynamic viscosity [kg/m/s]
                                                               ! layer_depth [m] = pressure [dbar]
        agg_env%rho_aq = rho_aq
        agg_env%mu     = visco(k)

        ! ------ prepare primary particle information to calculate aggregate properties
        !call prepare_primary_particles(tra, k)
        !call prepare_primary_particles(CC_det(k), 3.0e-14_wp, 0._wp, 0._wp)
        call prepare_primary_particles(CC_det(k), CC_opal(k), CC_calc(k), CC_dust(k), aggs, agg_env)

        ! ------ calculate aggregate properties from individual primary particle information
        call aggregate_properties(aggs, agg_env)
 
        ! ======== calculate the mean sinking velocity of aggregates =======
        call ws_Re_approx(aggs, agg_env)
 
        ! Limit settling velocity wrt CFL:
        wsagg(k) = min(aggs%ws_aggregates*dtbgc, 0.99_wp*pddpo(k)) ! (m/s -> m/d)*dtb
                                                               ! dtbgc: time step length in seconds
                                                               ! ws_agg [m/timestep]

        !tra%visco(k) = mol_dyn_vis(layer_depth(k), ptho(k), psao(k)) ! visco(pressure [dbar], temp [C], salinity [psu]) from m4ago_core

        if (wsagg(k) .lt. 0.0 .OR. wsagg(k)/dtb .gt. 530.0) then
          WRITE(*,*)'+++ WEIRD: aggregates rising / sinking faster than 530 m/day; aggregate speed [m/day] = ', wsagg(k)/dtb,' at k = ', k
        endif

        xLmaxagg(k)    = aggs%dmax_agg       ! applied max. diameter
        xavdp(k)       = aggs%av_dp          ! mean primary particle diameter
        xavrhop(k)     = aggs%av_rho_p       ! mean primary particle density
        xdfagg(k)      = aggs%df_agg         ! aggregate fractal dim
        xbagg(k)       = aggs%b_agg          ! aggre number distr. slope
        xavdC(k)       = conc_weighted_mean_agg_diameter(aggs)  ! conc-weighted mean agg. diameter
        xavrhofV(k)    = volweighted_agg_density(aggs,agg_env)  ! volume-weighted aggregate density
        xavporV(k)     = volweighted_agg_porosity(aggs)         ! volume-weighted aggregate porosity
        xstickagg(k)   = aggs%stickiness_agg           ! aggre. stickiness
        xstickfrust(k) = aggs%stickiness_frustule      ! frustle stickiness


!         tra%wpoc(k)                = ws_agg      ! particle sinking speeds [m/timestep]
!         tra%wopal(k)               = ws_agg      !   [m/timestep]
!         tra%wcal(k)                = ws_agg      !   [m/timestep]
!         tra%wdust(k)               = ws_agg      !   [m/timestep]
!         tra%visco(k)               = visco ! dyn_vis(k)     ! dynamic molecular viscosity
!         tra%Lmax_agg(k)            = Lmax_agg       ! applied max. diameter
!         tra%av_dp(k)               = av_dp          ! mean primary particle diameter
!         tra%av_rho_p(k)            = av_rho_p       ! mean primary particle density
!         tra%df_agg(k)              = df_agg         ! aggregate fractal dim
!         tra%b_agg(k)               = b_agg          ! aggre number distr. slope
!         tra%av_d_C(k)              = conc_weighted_mean_agg_diameter()  ! conc-weighted mean agg. diameter
!         tra%av_rhof_V(k)           = volweighted_agg_density()          ! volume-weighted aggregate density
!         tra%av_por_V(k)            = volweighted_agg_porosity()         ! volume-weighted aggregate porosity
!         tra%stickiness_agg(k)      = stickiness_agg           ! aggre. stickiness
!         tra%stickiness_frustule(k) = stickiness_frustule      ! frustle stickiness

    enddo
  end subroutine climberx_mean_aggregate_sinking_speed

  !===================================================================================== aggregate_properties
  ! FIXME subroutine prepare_primary_particles(k,tra,aggs,agg_env)
  subroutine prepare_primary_particles(C_det, C_opal, C_calc, C_dust, aggs, agg_env)
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

    !TYPE(carb_t), INTENT(in)  :: tra
    ! integer, intent(in)  :: k                  !< 3rd (vertical) real of model grid.
    type(aggregates),intent(inout)   :: aggs
    type(agg_environment),intent(in) :: agg_env

    real(wp),intent(in) :: C_det,C_opal,C_calc,C_dust         ! Concentration of tracers
    real(wp) :: n_det,n_opal,n_calc,n_dust         ! total primary particle number (#)
    real(wp) :: A_dust,A_det,A_calc,A_opal,A_total ! total surface area of primary particles per unit volume (L^2/L^3)
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

    V_det   = 0. ! total volume of primary particles in a unit volume
    V_opal  = 0.
    V_calc  = 0.
    V_dust  = 0.
    V_solid = 0.

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

    ! number of opal frustules (/NUM_FAC)
    ! FIXME: What if C_opal=0? 
    !        Better to set n_opal to at least 1/NUM_FAC?
    n_opal = C_opal*opalwei/rho_V_frustule_opal
    ! maximum mass of detritus inside a frustule
    cell_pot_det_mass = n_opal*V_frustule_inner*agg_org_dens

    ! detritus mass inside frustules
    ! cell_det_mass = min(cell_pot_det_mass, C_det*det_mol2mass - EPS_ONE)
    ! cell_det_mass = min(cell_pot_det_mass, C_det*det_mol2mass + EPS_ONE) ! "-" before, mistake?
    cell_det_mass = max(0._wp,min(cell_pot_det_mass,C_det*det_mol2mass))

    ! volume of detritus component in cell
    ! CAREFUL / fixme: (0/0)/agg_org_dens for C_opal=0
    V_POM_cell = (cell_det_mass/n_opal)/agg_org_dens

    ! if not detritus is available, water is added
    V_aq = V_frustule_inner -  V_POM_cell

    ! density of the diatom frsutules incl. opal, detritus and water
    rho_frustule = (rho_V_frustule_opal + cell_det_mass/n_opal + V_aq*agg_env%rho_aq)/V_dp_opal

    ! mass of extra cellular detritus particles
    free_detritus = C_det*det_mol2mass  - cell_det_mass
    rho_diatom = (rho_frustule + cell_det_mass/cell_pot_det_mass*rho_TEP)                          &
                   /(1._wp + cell_det_mass/cell_pot_det_mass)

    ! number of primary particles
    ! FIXME: See above for opal / better to set n_* to at least 1/NUM_FAC?
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

    ! calc frustule stickiness
    stickiness_frustule = cell_det_mass/(cell_pot_det_mass +EPS_ONE)*stickiness_TEP                &
                               & + (1._wp - cell_det_mass/(cell_pot_det_mass + EPS_ONE))           &
                               &   *stickiness_opal


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


end module mo_m4ago_climberx
