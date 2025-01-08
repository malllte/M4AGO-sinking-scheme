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
!! @file mo_m4ago_types.F90
!! @brief Module for Marine Aggregates:
!!        The Microstructure, Multiscale, Mechanistic, Marine Aggregates
!!        in the Global Ocean (M4AGO) sinking scheme
!!
!! The mo_m4ago_types module contains type definitions:
!!      - aggregate properties
!!      - aggregate environment
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

module mo_m4ago_types

  use mo_m4ago_kind, only: wp

  implicit none

  private

  ! type aggregates holds required information on aggregates, their composition, their size distribution, etc.
  type, public :: aggregates
    integer  :: NPrimPartTypes                       ! Number of primary particle types
    real(wp) :: av_dp                                ! mean primary particle diameter (m)
    real(wp) :: av_rho_p                             ! mean primary particle density (kg/m3)
    real(wp) :: df_agg                               ! aggregate fractal dimension - range: agg_df_min-agg_df_max (-)
    real(wp) :: b_agg                                ! aggregate number distribution slope (-)
    real(wp) :: dmax_agg                             ! maximum aggregate diameter (m)
    real(wp) :: stickiness_agg                       ! aggregate stickiness - range : 0-1 (-)
    real(wp) :: stickiness_frustule                  ! opal frustule stickiness
    real(wp) :: Re_crit_agg                          ! critical diameter-based particle Reynolds number for fragmentation
    real(wp) :: ws_aggregates                        ! mean aggregate sinking velocity (m/s)
    real(wp) :: n_pptotal                            ! total number of primary particles (#/L^3)
    real(wp),dimension(:), allocatable :: dp_pp         ! primary particle diameter of each primary particle type (L)
    real(wp),dimension(:), allocatable :: rho_pp        ! primary particle density of each primary particle type (M/L^3)
    real(wp),dimension(:), allocatable :: stickiness_pp ! stickiness of each primary particle type (-)
    real(wp),dimension(:), allocatable :: n_pp          ! total number of each primary particle type (#/L^3)
    real(wp),dimension(:), allocatable :: A_pp          ! surface area of each primary particle type (L^2/L^3)
    real(wp),dimension(:), allocatable :: V_pp          ! total volume of each primary particle type (L^3/L^3)
  end type aggregates

  ! type agg_environment holds information on the local environment of the aggregtes
  type, public ::agg_environment
    real(wp) :: mu                            ! molecular dynamic viscosity
    real(wp) :: rho_aq                        ! density of surrounding water
  end type agg_environment

end module mo_m4ago_types
