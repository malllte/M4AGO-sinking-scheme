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
module mo_driver_routines
  use mo_m4ago_kind,    only: wp
  use mo_m4ago_types,   only: aggregates,agg_environment
  use mo_m4ago_core,    only: volweighted_agg_density,                                             &
                            & volweighted_agg_porosity,conc_weighted_mean_agg_diameter
  use mo_m4ago_HAMOCCinit, only: NUM_FAC

  implicit none

  private

  public :: print_information

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

end module mo_driver_routines

