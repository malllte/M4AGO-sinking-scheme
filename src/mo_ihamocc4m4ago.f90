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
!! @file mo_ihamocc4m4ago.F90
!! @brief Module for Marine Aggregates:
!!        The Microstructure, Multiscale, Mechanistic, Marine Aggregates
!!        in the Global Ocean (M4AGO) sinking scheme
!!
!! The mo_ihamocc4m4ago module contains routines to calculate:
!!      - primary particles from iHAMOCC tracers
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
!!
!! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
!!
!!

module mo_ihamocc4m4ago

  ! iHAMOCC ocean biogeochemistry model-specific routines
  use mo_vgrid,         only: dp_min
  use mo_control_bgc,   only: dtb, dtbgc,io_stdo_bgc
  use mo_param_bgc,     only: calcdens, claydens, opaldens, calcwei, opalwei, ropal
  use mo_carbch,        only: ocetra
  use mo_param1_bgc,    only: iopal, ifdust, icalc, idet

  ! M4AGO routines:
  use mo_m4ago_core,    only: rho_aq,ONE_SIXTH,PI,aggregates,                                      &
                            & ws_Re_approx,volweighted_agg_density,                                &
                            & volweighted_agg_porosity,conc_weighted_mean_agg_diameter,            &
                            & aggregate_properties, init_m4ago_core_parameters

  use mo_m4ago_physics, only: mol_dyn_vis

  implicit none

  private

  ! Public subroutines
  public :: ihamocc_mean_aggregate_sinking_speed, init_m4ago_nml_params, init_m4ago_params,        &
          & alloc_mem_m4ago, cleanup_mem_m4ago


  ! Public fields and parameters
  public :: ws_agg,                                                                                &
          & aggregate_diagnostics,kav_dp,kav_rho_p,kav_d_C,kws_agg,kdf_agg,kstickiness_agg,kb_agg, &
          & kstickiness_frustule,kLmax_agg,kdynvis,kav_rhof_V,kav_por_V


  ! ------------------------------------------------------------------------------------------------
  ! biogeochemistry model-specific parameters
  ! ------------------------------------------------------------------------------------------------
  ! primary particle diameter for POM & PIM species involved in parametrized aggregation (m)
  real, protected :: dp_dust ! primary particle diameter dust
  real, protected :: dp_det  ! primary particle diameter detritus
  real, protected :: dp_calc ! primary particle diameter calc
  real, protected :: dp_opal ! primary particle diameter opal

  ! Stickiness of primary particles
  real, protected :: stickiness_TEP  ! stickiness of TEP (related to opal frustules)
  real, protected :: stickiness_det  ! normal detritus stickiness
  real, protected :: stickiness_opal ! stickiness of opal (without TEP - just normal coating)
  real, protected :: stickiness_calc ! stickiness of calc particles (coated with organics)
  real, protected :: stickiness_dust ! stickiness of dust particles (coated with organics)

  real, protected :: agg_df_max      ! maximum fractal dimension of aggregates (~2.5)
  real, protected :: agg_df_min      ! minimum fractal dimension of aggregates (~1.2 - 1.6)
  real, protected :: rho_TEP         ! density of TEP particles
  real, protected :: agg_org_dens    ! organic detritus density (alternative to orgdens to avoid negative ws)
  real, protected :: agg_Re_crit     ! critical particle Reynolds number for fragmentation

  ! calculated model-specific parameters
  real, protected :: det_mol2mass ! mol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)
  real, protected :: V_dp_dust,V_dp_det,V_dp_calc,V_dp_opal   ! volumes of primary particles (L^3)
  real, protected :: A_dp_dust,A_dp_det,A_dp_calc,A_dp_opal   ! surface areas of primary particles (L^2)
  real, protected :: stickiness_min, stickiness_max           ! minimum and maximum stickiness of primary particles
  real, protected :: rho_V_dp_dust,rho_V_dp_det,rho_V_dp_calc ! rho_V_dp_opal ! mass of primary particles (M)
  real, protected :: Rm_SiP                                   ! molar mass ratio opal (SiO_2) to POM
  real, protected :: thick_shell                              ! diatom frustule shell thickness (L)
  real, protected :: d_frustule_inner                         ! diameter of hollow part in diatom frustule (L)
  real, protected :: V_frustule_inner                         ! volume of hollow part in diatom frustule (L^3)
  real, protected :: V_frustule_opal                          ! volume of opal shell material (L^3)
  real, protected :: rho_V_frustule_opal                      ! mass of frustule material (M)

  ! Parameter for M4AGO core
  integer, parameter :: NPrimPartTypes = 4 ! Number of primary particle types generated from the biogeochemistry model

  ! Fields
  real,allocatable :: ws_agg(:,:,:)       ! mass concentration-weighted aggregate mean sinking velocity
  real,allocatable :: dyn_vis(:,:,:)      ! molecular dynamic viscosity
  real,allocatable :: m4ago_ppo(:,:,:)    ! pressure

  ! Marine aggregate diagnostics
  real, dimension (:,:,:,:), allocatable, target :: aggregate_diagnostics    ! diagnostics for marine aggregates

  integer, parameter :: kav_dp               =  1, &
                        kav_rho_p            =  2, &
                        kav_d_C              =  3, &
                        kws_agg              =  4, &
                        kdf_agg              =  5, &
                        kstickiness_agg      =  6, &
                        kb_agg               =  7, &
                        kstickiness_frustule =  8, &
                        kLmax_agg            =  9, &
                        kdynvis              = 10, &
                        kav_rhof_V           = 11, &
                        kav_por_V            = 12, &
                        naggdiag             = 12

  ! Internally used parameters and values
  real, parameter :: NUM_FAC = 1.e9             ! factor to avoid numerical precision problems
  real, parameter :: EPS_ONE = EPSILON(1.)

contains

  !===================================================================================== m4ago_init_params
  subroutine init_m4ago_nml_params
    !>
    !! Initialization of namelist parameters
    !!

    implicit none

    ! Primary particle sizes
    dp_dust = 2.e-6      ! following the classical HAMOCC parametrization
    dp_det  = 4.e-6      ! not well defined
    dp_calc = 3.e-6      ! following Henderiks 2008, Henderiks & Pagani 2008
    dp_opal = 20.e-6     ! mean frustule diameter of diatoms

    ! Stickiness values - note that their relative values to each other matter!
    stickiness_TEP    = 0.19
    stickiness_det    = 0.1
    stickiness_opal   = 0.08
    stickiness_calc   = 0.09
    stickiness_dust   = 0.07

    ! Minimum and maximum aggregate fractal dimension
    agg_df_min        = 1.6
    agg_df_max        = 2.4

    ! Density of primary particles
    rho_TEP           = 800. ! 700.-840. kg/m^3 Azetsu-Scott & Passow 2004
    agg_org_dens      = 1100. ! detritus density - don't use orgdens to avoid negative ws

    ! Critical particle Reynolds number for limiting nr-distribution
    agg_Re_crit       = 20.

  end subroutine init_m4ago_nml_params


  subroutine init_m4ago_params
    !>
    !! Initilization of parameters
    !!

    implicit none
    det_mol2mass   = 3166.  ! mol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)

    ! Volume of an individual primary particle*NUMFAC
    V_dp_dust = ONE_SIXTH*PI*dp_dust**3.*NUM_FAC
    V_dp_det  = ONE_SIXTH*PI*dp_det**3. *NUM_FAC
    V_dp_calc = ONE_SIXTH*PI*dp_calc**3.*NUM_FAC
    V_dp_opal = ONE_SIXTH*PI*dp_opal**3.*NUM_FAC

    ! Surface area of an individual primary particle*NUMFAC
    A_dp_dust = PI*dp_dust**2.*NUM_FAC
    A_dp_det  = PI*dp_det**2. *NUM_FAC
    A_dp_calc = PI*dp_calc**2.*NUM_FAC
    A_dp_opal = PI*dp_opal**2.*NUM_FAC

    ! Mass of an individual primary particle*NUMFAC
    rho_V_dp_dust = V_dp_dust*claydens
    rho_V_dp_det  = V_dp_det*agg_org_dens
    rho_V_dp_calc = V_dp_calc*calcdens

    Rm_SiP              = ropal*opalwei/det_mol2mass
    ! shell thickness
    thick_shell         = 0.5*dp_opal*(1. - (opaldens/(Rm_SiP*agg_org_dens+opaldens))**(1./3.))
    d_frustule_inner    = dp_opal - 2.*thick_shell
    ! volume of hollow part of frustule
    V_frustule_inner    = ONE_SIXTH* PI*d_frustule_inner**3.*NUM_FAC
    ! volume of opal part of frustule
    V_frustule_opal     = ONE_SIXTH*PI*(dp_opal**3. - d_frustule_inner**3.)*NUM_FAC
    rho_V_frustule_opal = V_frustule_opal*opaldens

    ! Minimum and maximum reachable stickiness
    stickiness_min      = min(stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust)
    stickiness_max      = max(stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust)

    ! Init core M4AGO parameters
    call init_m4ago_core_parameters(agg_Re_crit,agg_df_min,agg_df_max,stickiness_min,stickiness_max)

  end subroutine init_m4ago_params

  !===================================================================================== mean_agg_ws
  subroutine ihamocc_mean_aggregate_sinking_speed(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, ppao, prho)
    !-----------------------------------------------------------------------
    !>
    !! calculates the mass concentration-weighted mean sinking velocity of marine
    !! aggregates
    !!

    implicit none

    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.
    integer, intent(in)  :: kbnd                  !< time step boundary
    real,    intent(in)  :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
    real,    intent(in)  :: omask(kpie,kpje)      !< ocean mask
    real,    intent(in)  :: ptho (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< potential temperature [deg C]
    real,    intent(in)  :: psao (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< salinity [psu.].
    real,    intent(in)  :: ppao (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)      !< pressure at sea level [Pa].
    real,    intent(in)  :: prho (kpie,kpje,kpke)                         !< water density [g/cm3]

    integer :: i,j,k
    type(aggregates) :: aggs
    allocate(aggs%dp_pp(NPrimPartTypes))
    allocate(aggs%rho_pp(NPrimPartTypes))
    allocate(aggs%stickiness_pp(NPrimPartTypes))
    allocate(aggs%n_pp(NPrimPartTypes))
    allocate(aggs%A_pp(NPrimPartTypes))
    allocate(aggs%V_pp(NPrimPartTypes))

    ! get pressure
    call calc_pressure(kpie, kpje, kpke,kbnd, pddpo, omask)

    ! molecular dynamic viscosity
    call dynvis(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, m4ago_ppo)

    !$OMP PARALLEL DO PRIVATE(i,j,k,aggs)
    do j = 1,kpje
      do i = 1,kpie
        do k = 1,kpke
          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then

            ! ------ prepare primary particle information to calculate aggregate properties
            call prepare_primary_particles(i, j, k,aggs)

            ! ------ calculate aggregate properties from individual primary particle information
            call aggregate_properties(aggs, dyn_vis(i,j,k))

            ! ======== calculate the mean sinking velocity of aggregates =======
            call ws_Re_approx(aggs,dyn_vis(i,j,k))

            ! Limit settling velocity wrt CFL:
            ws_agg(i,j,k) = min(aggs%ws_aggregates*dtbgc, 0.99*pddpo(i,j,k)) ! (m/s -> m/d)*dtb


            ! ============================== Write general diagnostics ============
            aggregate_diagnostics(i,j,k,kws_agg)    = ws_agg(i,j,k)/dtb  ! applied ws conversion  m/time_step  to  m/d for output
            aggregate_diagnostics(i,j,k,kdynvis)    = dyn_vis(i,j,k)     ! dynamic molecular viscosity
            aggregate_diagnostics(i,j,k,kLmax_agg)  = aggs%dmax_agg      ! applied max. diameter
            aggregate_diagnostics(i,j,k,kav_dp)     = aggs%av_dp         ! mean primary particle diameter
            aggregate_diagnostics(i,j,k,kav_rho_p)  = aggs%av_rho_p      ! mean primary particle density
            aggregate_diagnostics(i,j,k,kdf_agg)    = aggs%df_agg        ! aggregate fractal dim
            aggregate_diagnostics(i,j,k,kb_agg)     = aggs%b_agg         ! aggre number distr. slope
            aggregate_diagnostics(i,j,k,kav_d_C)    = conc_weighted_mean_agg_diameter(aggs) ! conc-weighted mean agg. diameter
            aggregate_diagnostics(i,j,k,kav_rhof_V) = volweighted_agg_density(aggs)         ! volume-weighted aggregate density
            aggregate_diagnostics(i,j,k,kav_por_V)  = volweighted_agg_porosity(aggs)        ! volume-weighted aggregate porosity
            aggregate_diagnostics(i,j,k,kstickiness_agg)      = aggs%stickiness_agg         ! aggre. stickiness
            aggregate_diagnostics(i,j,k,kstickiness_frustule) = aggs%stickiness_frustule    ! frustle stickiness
          endif
        enddo
      enddo
    enddo
  end subroutine ihamocc_mean_aggregate_sinking_speed

  !===================================================================================== aggregate_properties
  subroutine prepare_primary_particles(i, j, k,aggs)
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

    integer, intent(in)  :: i                  !< 1st real of model grid.
    integer, intent(in)  :: j                  !< 2nd real of model grid.
    integer, intent(in)  :: k                  !< 3rd (vertical) real of model grid.
    type(aggregates),intent(inout) :: aggs

    real :: C_det,C_opal,C_calc,C_dust         ! Concentration of tracers
    real :: n_det,n_opal,n_calc,n_dust         ! total primary particle number (#)
    real :: A_dust,A_det,A_calc,A_opal,A_total ! total surface area of primary particles per unit volume (L^2/L^3)
    real :: V_det,V_opal,V_calc,V_dust,V_solid ! total volume of primary particles in a unit volume (L^3/L^3)

    real :: cell_det_mass                      ! mass of detritus material in diatoms
    real :: cell_pot_det_mass                  ! potential (max) mass detritus material in diatoms
    real :: free_detritus                      ! freely available detritus mass outside the frustule
    real :: V_POM_cell                         ! volume of POM in frustule
    real :: V_aq                               ! volume of water space in frustule
    real :: rho_frustule                       ! density of diatom frustule incl. opal, detritus and water
    real :: rho_diatom                         ! density of either hollow frustule or with additions of detritus and water
    real :: stickiness_frustule                ! stickiness of the diatom frustile as primary particle

    C_det  = abs(ocetra(i,j,k,idet))
    C_opal = abs(ocetra(i,j,k,iopal))
    C_calc = abs(ocetra(i,j,k,icalc))
    C_dust = abs(ocetra(i,j,k,ifdust))

    n_det   = 0. ! number of primary particles in a unit volume
    n_opal  = 0.
    n_dust  = 0.
    n_calc  = 0.

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
    free_detritus     = 0.
    rho_diatom        = 0.
    cell_det_mass     = 0.
    cell_pot_det_mass = 0.
    V_POM_cell        = 0.
    V_aq              = 0.
    rho_frustule      = 0.

    ! number of opal frustules (/NUM_FAC)
    n_opal = C_opal*opalwei/rho_V_frustule_opal
    ! maximum mass of detritus inside a frustule
    cell_pot_det_mass = n_opal*V_frustule_inner*agg_org_dens

    ! detritus mass inside frustules
    cell_det_mass = min(cell_pot_det_mass, C_det*det_mol2mass - EPS_ONE)

    ! volume of detritus component in cell
    V_POM_cell = (cell_det_mass/n_opal)/agg_org_dens

    ! if not detritus is available, water is added
    V_aq = V_frustule_inner -  V_POM_cell

    ! density of the diatom frsutules incl. opal, detritus and water
    rho_frustule = (rho_V_frustule_opal + cell_det_mass/n_opal + V_aq*rho_aq)/V_dp_opal

    ! mass of extra cellular detritus particles
    free_detritus = C_det*det_mol2mass  - cell_det_mass
    rho_diatom = (rho_frustule + cell_det_mass/cell_pot_det_mass*rho_TEP)                          &
                   /(1. + cell_det_mass/cell_pot_det_mass)

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

    ! calc frustule stickiness
    stickiness_frustule = cell_det_mass/(cell_pot_det_mass +EPS_ONE)*stickiness_TEP           &
                               & + (1. - cell_det_mass/(cell_pot_det_mass + EPS_ONE))              &
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

  !===================================================================================== pressure
  subroutine alloc_mem_m4ago(kpie, kpje, kpke)
    !-----------------------------------------------------------------------
    !>
    !! Initialization/allocation fields
    !! Called in ini_bgc after read_namelist
    !!

    implicit none

    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.

    ! allocate memory space for:
    ! -> aggregate diagnostics (output)
    allocate(aggregate_diagnostics(kpie, kpje, kpke, naggdiag))

    ! -> mean sinking velocity
    allocate(ws_agg(kpie,kpje,kpke))

    ! -> molecular dynamic viscosity
    allocate(dyn_vis(kpie, kpje, kpke))

    ! -> pressure
    allocate(m4ago_ppo(kpie,kpje,kpke))

    ! Initialization
    aggregate_diagnostics = 0.
    m4ago_ppo             = 0.
    ws_agg                = 0.
    dyn_vis               = 0.

  end subroutine alloc_mem_m4ago

  !===================================================================================== pressure
  subroutine cleanup_mem_m4ago
    implicit none
    ! clean memory
    deallocate(aggregate_diagnostics)
    deallocate(ws_agg)
    deallocate(dyn_vis)
    deallocate(m4ago_ppo)
  end subroutine cleanup_mem_m4ago

  !===================================================================================== pressure
  subroutine calc_pressure(kpie, kpje, kpke,kbnd, pddpo,omask)

    use mo_vgrid, only: ptiestu

    implicit none

    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.
    integer, intent(in)  :: kbnd
    real, intent(in) :: pddpo(kpie,kpje,kpke)     !< size of scalar grid cell (3rd dimension) [m]
    real, intent(in) :: omask(kpie,kpje)          !< mask

    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = 1,kpke
      do j = 1,kpje
        do i = 1,kpie
          if(omask(i,j) > 0.5 .and. pddpo(i,j,k) .gt. dp_min) then
            m4ago_ppo(i,j,k) = 1e5 * ptiestu(i,j,k)*98060.*1.027e-6 ! pressure in unit Pa, 98060 = onem
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_pressure

  !===================================================================================== dynvis
  subroutine dynvis(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, ppo)
    !-----------------------------------------------------------------------
    !>
    !! dynvis calculates the molecular dynamic viscosity according to
    !! Richards 1998: The effect of temperature, pressure, and salinity
    !! on sound attenuation in turbid seawater. J. Acoust. Soc. Am. 103 (1),
    !! originally published by  Matthaeus, W. (1972): Die Viskositaet des
    !! Meerwassers. Beitraege zur Meereskunde, Heft 29 (in German).
    !!

    implicit none

    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.
    integer, intent(in)  :: kbnd

    real, intent(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
    real, intent(in) :: omask(kpie,kpje)      !< ocean mask
    real, intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< potential temperature [deg C]
    real, intent(in) :: psao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)  !< salinity [psu.].
    real, intent(in) :: ppo(kpie,kpje,kpke)  !< pressure [Pa].

    ! Local variables
    real    :: press_val  ! Pascal/rho -> dbar
    real    :: ptho_val,psao_val
    integer :: i,j,k,kch
    kch = 0
    !$OMP PARALLEL DO PRIVATE(i,j,k,press_val,ptho_val,psao_val,kch)
    do j = 1,kpje
      do i = 1,kpie
        do k = 1,kpke
          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
            kch = merge(k+1,k,k<kpke)
            if(pddpo(i,j,kch) > 0.5) then
              press_val    = 0.5*(ppo(i,j,k)  + ppo(i,j,kch))*1.e-5 ! Pascal -> dbar
              ptho_val     = 0.5*(ptho(i,j,k) + ptho(i,j,kch))
              psao_val     = 0.5*(psao(i,j,k) + ptho(i,j,kch))
            else
              press_val    = ppo(i,j,k)*1.e-5 ! Pascal -> dbar
              ptho_val     = ptho(i,j,k)
              psao_val     = psao(i,j,k)
            endif

            ! molecular dynamic viscosity
            dyn_vis(i,j,k) = mol_dyn_vis(press_val,ptho_val,psao_val)
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine dynvis


end module mo_ihamocc4m4ago
