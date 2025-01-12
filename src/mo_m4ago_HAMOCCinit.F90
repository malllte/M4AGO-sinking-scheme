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
module mo_m4ago_HAMOCCinit

  use mo_m4ago_kind,    only: wp
  use mo_m4ago_params,  only: rho_aq,ONE_SIXTH,PI
  use mo_m4ago_core,    only: init_m4ago_core_parameters

  implicit none

  private

  public :: init_m4ago_nml_params, init_m4ago_derived_params !,reset_m4ago_nml_params
  public :: agg_df_max,agg_df_min,agg_Re_crit,NPrimPartTypes,NUM_FAC,                              &
          & stickiness_min,stickiness_max,                                                         &
          & stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust,     &
          & dp_det,dp_opal,dp_calc,dp_dust,                                                        &
          & A_dp_det,A_dp_opal,A_dp_calc,A_dp_dust,                                                &
          & V_dp_det,V_dp_opal,V_dp_calc,V_dp_dust,V_frustule_inner,                               &
          & rho_TEP,rho_V_frustule_opal,rho_V_dp_det,rho_V_dp_calc,rho_V_dp_dust,                  &
          & det_mol2mass,rho_det,rho_calc,calc_weight,rho_dust,opal_weight


  ! ------------------------------------------------------------------------------------------------
  ! Biogeochemistry model-specific parameters - HAMOCC
  ! ------------------------------------------------------------------------------------------------
  ! Primary particle diameter for POM & PIM species involved in parametrized aggregation
  real(wp), protected :: dp_dust ! (m) Primary particle diameter dust
  real(wp), protected :: dp_det  ! (m) Primary particle diameter detritus
  real(wp), protected :: dp_calc ! (m) Primary particle diameter calc
  real(wp), protected :: dp_opal ! (m) Primary particle diameter opal

  ! Stickiness of primary particles
  real(wp), protected :: stickiness_TEP  ! (-) Stickiness of TEP (related to opal frustules)
  real(wp), protected :: stickiness_det  ! (-) Normal detritus stickiness
  real(wp), protected :: stickiness_opal ! (-) Stickiness of opal (without TEP - just normal coating)
  real(wp), protected :: stickiness_calc ! (-) Stickiness of calc particles (coated with organics)
  real(wp), protected :: stickiness_dust ! (-) Stickiness of dust particles (coated with organics)

  ! Minimum and maximum fractal dimension
  real(wp), protected :: agg_df_max      ! (-) Maximum fractal dimension of aggregates (~2.5)
  real(wp), protected :: agg_df_min      ! (-) Minimum fractal dimension of aggregates (~1.2 - 1.6)

  ! Base densities and mol-weights for primary particles
  real(wp), protected :: rho_TEP         ! (kg/m3)   Density of TEP particles
  real(wp), protected :: rho_det         ! (kg/m3)   Organic detritus density (alternative to orgdens to avoid negative ws)
  real(wp), protected :: rho_calc        ! (kg/m3)   Density of CaCO3
  real(wp), protected :: rho_opal        ! (kg/m3)   Density of opal
  real(wp), protected :: rho_dust        ! (kg/m3)   Density of dust/clay
  real(wp), protected :: opal_weight     ! (kg/kmol) mol-weight of opal
  real(wp), protected :: calc_weight     ! (kg/kmol) mol-weight of CaCO3

  ! Citical particle Reynolds number for fragmentation
  real(wp)            :: agg_Re_crit     ! (-)       Critical particle Reynolds number for fragmentation

  ! Calculated model-specific parameters - some include NUM_FAC
  real(wp), protected :: det_mol2mass ! (kg/kmol detritus P) mol-weight of detritus (according to HAMOCC stoichiometry)
  real(wp), protected :: V_dp_dust,V_dp_det,V_dp_calc,V_dp_opal   ! (m^3)*NUM_FAC  Volumes of primary particles
  real(wp), protected :: A_dp_dust,A_dp_det,A_dp_calc,A_dp_opal   ! (m^2)*NUM_FAC  Surface areas of primary particles
  real(wp), protected :: stickiness_min, stickiness_max           ! (-)            Minimum and maximum stickiness of primary particles
  real(wp), protected :: rho_V_dp_dust,rho_V_dp_det,rho_V_dp_calc ! (kg)           Mass of primary particle types
  real(wp), protected :: Rm_SiP                                   ! (kg Si/kg det) Molar mass ratio opal (SiO_2) to POM
  real(wp), protected :: thick_shell                              ! (m)            Diatom frustule shell thickness
  real(wp), protected :: d_frustule_inner                         ! (m)            Diameter of hollow part in diatom frustule
  real(wp), protected :: V_frustule_inner                         ! (m^3)*NUM_FAC  Volume of hollow part in diatom frustule
  real(wp), protected :: V_frustule_opal                          ! (m^3)*NUM_FAC  Volume of opal shell material
  real(wp), protected :: rho_V_frustule_opal                      ! (kg)*NUM_FAC   Mass of frustule material

  real(wp), protected :: NUM_FAC  = 1e9_wp                        ! (-) Numerical factor to avoid numerical issues

  integer, protected  :: NPrimPartTypes = 4 ! Number of primary particle types generated from the biogeochemistry model

contains
  ! ------------------------------------------------------------------------------------------------
  subroutine init_m4ago_derived_params(ropal)
    !>
    !! Initilization of parameters
    !!

    implicit none

    real(wp), intent(in) :: ropal ! (mol Si/mol P) Opal to phosporus production ratio

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
    rho_V_dp_dust = V_dp_dust*rho_dust
    rho_V_dp_det  = V_dp_det*rho_det
    rho_V_dp_calc = V_dp_calc*rho_calc

    Rm_SiP              = ropal*opal_weight/det_mol2mass
    ! Shell thickness
    thick_shell         = 0.5_wp*dp_opal*(1._wp - (rho_opal/(Rm_SiP*rho_det+rho_opal))**(1._wp/3._wp))
    d_frustule_inner    = dp_opal - 2._wp*thick_shell
    ! Volume of hollow part of frustule * NUM_FAC
    V_frustule_inner    = ONE_SIXTH* PI*d_frustule_inner**3*NUM_FAC
    ! Volume of opal part of frustule * NUM_FAC
    V_frustule_opal     = ONE_SIXTH*PI*(dp_opal**3 - d_frustule_inner**3)*NUM_FAC
    ! Mass of opal part of frustule * NUM_FAC
    rho_V_frustule_opal = V_frustule_opal*rho_opal

    ! Minimum and maximum reachable stickiness
    stickiness_min      = min(stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust)
    stickiness_max      = max(stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust)

    ! Init core M4AGO parameters
    call init_m4ago_core_parameters(agg_Re_crit,agg_df_min,agg_df_max,stickiness_min,stickiness_max)

  end subroutine init_m4ago_derived_params

!  ! ------------------------------------------------------------------------------------------------
!  subroutine reset_m4ago_nml_params()
!  ! could be introduced and used to hard-reset parameters after init_m4ago_nml_params
!  ! to enable setting of nml parameters directly through the driving OBGC model
!  ! currently left blank
!  implicit none
!
!  end subroutine

  ! ------------------------------------------------------------------------------------------------
  subroutine init_m4ago_nml_params(claydens,calcdens,calcwei,opaldens,opalwei,nmlfile)
    !>
    !! Initialization of (namelist) parameters
    !!

    implicit none

    real(wp), intent(in) :: claydens  ! (kg/m3)   Density of dust/clay
    real(wp), intent(in) :: calcdens  ! (kg/m3)   Density of CaCO3
    real(wp), intent(in) :: calcwei   ! (kg/kmol) mol-weight of CaCO3
    real(wp), intent(in) :: opaldens  ! (kg/m3)   Density of opal
    real(wp), intent(in) :: opalwei   ! (kg/kmol) mol-weight of opal
    character(*),intent(in),optional :: nmlfile ! namelist file for M4AGO

    integer :: iounit
    namelist /m4agoparams/ dp_dust,dp_det,dp_calc,dp_opal,rho_calc,rho_opal,rho_dust,rho_det,      &
                           rho_TEP,agg_Re_crit,stickiness_det,stickiness_opal,stickiness_calc,     &
                           stickiness_dust,stickiness_TEP,agg_df_min,agg_df_max,calc_weight,       &
                           opal_weight

    ! Detritus mol-weight for HAMOCC stoichiometry:
    det_mol2mass   = 3166._wp  ! (kg/kmol) kmol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)

    ! Densities and mol-weights of primary particles provided through HAMOCC
    rho_calc    = calcdens
    rho_opal    = opaldens
    rho_dust    = claydens
    opal_weight = opalwei
    calc_weight = calcwei

    ! Density of primary particles defined via M4AGO
    rho_TEP     = 800._wp  ! (kg/m^3) 700.-840. kg/m^3 Azetsu-Scott & Passow 2004
    rho_det     = 1100._wp ! (kg/m^3) detritus density - don't use orgdens to avoid negative ws

    ! Primary particle diameters
    dp_dust = 2.e-6_wp      ! (m) following the classical HAMOCC parametrization
    dp_det  = 4.e-6_wp      ! (m) not well defined
    dp_calc = 3.e-6_wp      ! (m) following Henderiks 2008, Henderiks & Pagani 2008
    dp_opal = 20.e-6_wp     ! (m) mean frustule diameter of diatoms

    ! Stickiness values - note that their relative values to each other matter!
    stickiness_TEP    = 0.19_wp ! (-) range 0-1
    stickiness_det    = 0.1_wp  ! (-) range 0-1
    stickiness_opal   = 0.08_wp ! (-) range 0-1
    stickiness_calc   = 0.09_wp ! (-) range 0-1
    stickiness_dust   = 0.07_wp ! (-) range 0-1

    ! Minimum and maximum aggregate fractal dimension
    agg_df_min        = 1.6_wp ! (-)
    agg_df_max        = 2.4_wp ! (-)

    ! Critical particle Reynolds number (based on diameter) for limiting nr-distribution
    agg_Re_crit       = 20._wp ! (-)

    if (present(nmlfile)) then
      ! if a namelist file is provided for M4AGO, read it
      print*,'Reading namelist file ...',nmlfile
      open (newunit=iounit, file=nmlfile,status='old',action='read')
      read (unit=iounit,nml=M4AGOPARAMS)
      close(unit=iounit)
      print*,'... successful'
    else
      print*,'No namelist file for M4AGO provided - continue with default parameters for HAMOCC'
    endif

  end subroutine init_m4ago_nml_params

end module mo_m4ago_HAMOCCinit
