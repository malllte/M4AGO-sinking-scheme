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
use mo_driver_routines,only: agg_df_max,agg_df_min,agg_Re_crit,dynvis,NPrimPartTypes,stickiness_max,&
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
