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
!! @file mo_m4ago_physics.f90
!! @brief Module for Marine Aggregates:
!!        The Microstructure, Multiscale, Mechanistic, Marine Aggregates
!!        in the Global Ocean (M4AGO) sinking scheme
!!
!! The mo_m4ao_physics module contains routines to calculate:
!!      - molecular dynamic viscosity
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
!! 2024 packaged as individual module by joeran maerz, UiB, Bergen
!!
!! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
!!
!!

module mo_m4ago_physics

  implicit none

  private

  public :: mol_dyn_vis

contains

  real function mol_dyn_vis(press_val,ptho_val,psao_val)
    !-----------------------------------------------------------------------
    !>
    !! mol_dyn_vis calculates the molecular dynamic viscosity according to
    !! Richards 1998: The effect of temperature, pressure, and salinity
    !! on sound attenuation in turbid seawater. J. Acoust. Soc. Am. 103 (1),
    !! originally published by  Matthaeus, W. (1972): Die Viskositaet des
    !! Meerwassers. Beitraege zur Meereskunde, Heft 29 (in German).
    !!
    !! returns: molecular dynamic viscosity [kg/(m*s)]

    real, intent(in) :: press_val ! pressure [dbar]
    real, intent(in) :: ptho_val  ! temperature [deg C]
    real, intent(in) :: psao_val  ! salinity [psu]

    mol_dyn_vis = 0.1 & ! Unit: g / (cm*s) -> kg / (m*s)
                      &     *(1.79e-2                                                         &
                      &     - 6.1299e-4*ptho_val + 1.4467e-5*ptho_val**2.                     &
                      &     - 1.6826e-7*ptho_val**3.                                          &
                      &     - 1.8266e-7*press_val  + 9.8972e-12*press_val**2.                 &
                      &     + 2.4727e-5*psao_val                                              &
                      &     + psao_val*(4.8429e-7*ptho_val - 4.7172e-8*ptho_val**2.           &
                      &     + 7.5986e-10*ptho_val**3.)                                        &
                      &     + press_val*(1.3817e-8*ptho_val - 2.6363e-10*ptho_val**2.)        &
                      &     - press_val**2.*(6.3255e-13*ptho_val - 1.2116e-14*ptho_val**2.))



  end function mol_dyn_vis

end module mo_m4ago_physics
