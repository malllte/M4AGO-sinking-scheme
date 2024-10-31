program M4AGO_driver

use mo_m4ago_kind,    only: wp
use mo_m4ago_core,    only: rho_aq,ONE_SIXTH,PI,aggregates,agg_environment,                        &
                          & ws_Re_approx,volweighted_agg_density,                                  &
                          & volweighted_agg_porosity,conc_weighted_mean_agg_diameter,              &
                          & aggregate_properties, init_m4ago_core_parameters

!use mo_m4ago_physics, only: mol_dyn_vis

  ! ------------------------------------------------------------------------------------------------
  ! biogeochemistry model-specific parameters
  ! ------------------------------------------------------------------------------------------------
  ! primary particle diameter for POM & PIM species involved in parametrized aggregation (m)
  real(wp) :: dp_dust ! primary particle diameter dust
  real(wp) :: dp_det  ! primary particle diameter detritus
  real(wp) :: dp_calc ! primary particle diameter calc
  real(wp) :: dp_opal ! primary particle diameter opal

  ! Stickiness of primary particles
  real(wp) :: stickiness_TEP  ! stickiness of TEP (related to opal frustules)
  real(wp) :: stickiness_det  ! normal detritus stickiness
  real(wp) :: stickiness_opal ! stickiness of opal (without TEP - just normal coating)
  real(wp) :: stickiness_calc ! stickiness of calc particles (coated with organics)
  real(wp) :: stickiness_dust ! stickiness of dust particles (coated with organics)

  real(wp) :: agg_df_max      ! maximum fractal dimension of aggregates (~2.5)
  real(wp) :: agg_df_min      ! minimum fractal dimension of aggregates (~1.2 - 1.6)
  real(wp) :: rho_TEP         ! density of TEP particles
  real(wp) :: agg_org_dens    ! organic detritus density (alternative to orgdens to avoid negative ws)
  real(wp) :: agg_Re_crit     ! critical particle Reynolds number for fragmentation

  ! calculated model-specific parameters
  real(wp) :: det_mol2mass ! mol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)
  real(wp) :: V_dp_dust,V_dp_det,V_dp_calc,V_dp_opal   ! volumes of primary particles (L^3)
  real(wp) :: A_dp_dust,A_dp_det,A_dp_calc,A_dp_opal   ! surface areas of primary particles (L^2)
  real(wp) :: stickiness_min, stickiness_max           ! minimum and maximum stickiness of primary particles
  real(wp) :: rho_V_dp_dust,rho_V_dp_det,rho_V_dp_calc ! rho_V_dp_opal ! mass of primary particles (M)
  real(wp) :: Rm_SiP                                   ! molar mass ratio opal (SiO_2) to POM
  real(wp) :: thick_shell                              ! diatom frustule shell thickness (L)
  real(wp) :: d_frustule_inner                         ! diameter of hollow part in diatom frustule (L)
  real(wp) :: V_frustule_inner                         ! volume of hollow part in diatom frustule (L^3)
  real(wp) :: V_frustule_opal                          ! volume of opal shell material (L^3)
  real(wp) :: rho_V_frustule_opal                      ! mass of frustule material (M)

  ! Parameter for M4AGO core
  integer, parameter :: NPrimPartTypes = 4 ! Number of primary particle types generated from the biogeochemistry model


  real(wp) :: dynvis  = 0.001567 ! [kg/(m s)] dynamic molecular viscosity
  real(wp) :: calcdens = 2600.
  real(wp) :: claydens=2600.
  real(wp) :: NUM_FAC=1e9
  real(wp) :: opaldens = 2200.
  real(wp) :: opalwei = 60.
  real(wp) :: ropal=20.

  integer i,j,k

  call init_m4ago_nml_params
  call init_m4ago_params




contains

  subroutine init_m4ago_nml_params
    !>
    !! Initialization of namelist parameters
    !!

    implicit none

    ! Primary particle sizes
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


end program M4AGO_driver
