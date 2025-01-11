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
module mo_m4ago_HAMOCCPrimPart

  use mo_m4ago_kind,  only: wp
  use mo_m4ago_types, only: aggregates, agg_environment
  use mo_m4ago_HAMOCCinit, only: NUM_FAC,det_mol2mass,NPrimPartTypes,                              &
                                 rho_det,rho_calc,calc_weight,rho_dust,opal_weight,rho_TEP,   &
                                 dp_det,dp_opal,dp_calc,dp_dust,                                   &
                                 A_dp_det,A_dp_calc,A_dp_dust,A_dp_opal,                           &
                                 V_dp_det,V_dp_opal,V_dp_calc,V_dp_dust,                           &
                                 rho_V_frustule_opal,V_frustule_inner,rho_V_dp_det,                &
                                 rho_V_dp_calc,rho_V_dp_dust,                                      &
                                 stickiness_TEP, stickiness_det, stickiness_opal,                  &
                                 stickiness_calc, stickiness_dust

  implicit none

  private

  public :: prepare_primary_particles

contains

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

    real(wp), intent(in)  :: C_det                  ! (kmol P/m3)  Detritus concentration
    real(wp), intent(in)  :: C_opal                 ! (kmol Si/m3) Opal concentration
    real(wp), intent(in)  :: C_calc                 ! (kmol C/m3)  CaCO3 concentration
    real(wp), intent(in)  :: C_dust                 ! (kg/m3)      Dust concentration
    type(aggregates),intent(inout)   :: aggs
    type(agg_environment),intent(in) :: agg_env

    real(wp) :: n_det,n_opal,n_calc,n_dust         ! (#/m3)  Total primary particle number of each type
    real(wp) :: A_dust,A_det,A_calc,A_opal         ! (m2/m3) Total surface area of primary particles per unit volume
    real(wp) :: V_det,V_opal,V_calc,V_dust,V_solid ! (m3/m3) Total volume of primary particles in a unit volume

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
    n_opal = C_opal*opal_weight/rho_V_frustule_opal

    ! maximum mass of detritus inside a frustule
    cell_pot_det_mass = n_opal*V_frustule_inner*rho_det

    ! detritus mass inside frustules
    cell_det_mass = max(0._wp,min(cell_pot_det_mass,C_det*det_mol2mass))

    if (n_opal > 0._wp) then
      ! volume of detritus component in diatom cell
      V_POM_cell = (cell_det_mass/n_opal)/rho_det

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
    n_calc = C_calc*calc_weight/rho_V_dp_calc
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
    aggs%dp_pp          = (/ dp_dust,  dp_calc,  dp_det,   dp_opal /)
    aggs%rho_pp         = (/ rho_dust, rho_calc, rho_det,  rho_diatom /)
    aggs%n_pp           = (/ n_dust,   n_calc,   n_det,    n_opal  /)
    aggs%A_pp           = (/ A_dust,   A_calc,   A_det,    A_opal /)
    aggs%V_pp           = (/ V_dust,   V_calc,   V_det,    V_opal /)
    aggs%stickiness_pp  = (/ stickiness_dust, stickiness_calc, stickiness_det, stickiness_frustule /)
    aggs%stickiness_frustule = stickiness_frustule

  end subroutine prepare_primary_particles

end module mo_m4ago_HAMOCCPrimPart
