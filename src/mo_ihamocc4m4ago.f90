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
  use mo_carbch,        only: ocetra
  use mo_param1_bgc,    only: iopal, ifdust, icalc, idet

  ! M4AGO routines:
  use mo_m4ago_kind,    only: wp
  use mo_m4ago_types,   only: aggregates,agg_environment
  use mo_m4ago_params,  only: rho_aq
  use mo_m4ago_core,    only: mean_aggregate_sinking_speed,volweighted_agg_density,                &
                            & volweighted_agg_porosity,conc_weighted_mean_agg_diameter,            &
                            & aggregate_properties, init_m4ago_core_parameters

  use mo_m4ago_physics, only: mol_dyn_vis
  use mo_m4ago_HAMOCCinit,     only: NPrimPartTypes
  use mo_m4ago_HAMOCCPrimPart, only: prepare_primary_particles

  implicit none

  private

  ! Public subroutines
  public :: ihamocc_mean_aggregate_sinking_speed,alloc_mem_m4ago, cleanup_mem_m4ago

  ! Public fields and parameters
  public :: ws_agg,                                                                                &
          & aggregate_diagnostics,kav_dp,kav_rho_p,kav_d_C,kws_agg,kdf_agg,kstickiness_agg,kb_agg, &
          & kstickiness_frustule,kLmax_agg,kdynvis,kav_rhof_V,kav_por_V

  ! Fields
  real(wp),allocatable :: ws_agg(:,:,:)       ! mass concentration-weighted aggregate mean sinking velocity
  real(wp),allocatable :: dyn_vis(:,:,:)      ! molecular dynamic viscosity
  real(wp),allocatable :: m4ago_ppo(:,:,:)    ! pressure

  ! Marine aggregate diagnostics
  real(wp), dimension (:,:,:,:), allocatable, target :: aggregate_diagnostics    ! diagnostics for marine aggregates

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

contains

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
    type(agg_environment) :: agg_env
    type(aggregates)      :: aggs

    ! Allocate memory for primary particle types information
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

    !$OMP PARALLEL DO PRIVATE(i,j,k,aggs,agg_env)
    do j = 1,kpje
      do i = 1,kpie
        do k = 1,kpke
          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5_wp) then

            ! Provide aggregates environment
            agg_env%rho_aq = rho_aq
            agg_env%mu     = dyn_vis(i,j,k)

            ! ------ prepare primary particle information to calculate aggregate properties
            call prepare_primary_particles(abs(ocetra(i,j,k,idet)),abs(ocetra(i,j,k,iopal)),       &
                                           abs(ocetra(i,j,k,icalc)),abs(ocetra(i,j,k,ifdust)),     &
                                           aggs,agg_env)

            ! ------ calculate aggregate properties from individual primary particle information
            call aggregate_properties(aggs, agg_env)

            ! ======== calculate the mean sinking velocity of aggregates =======
            call mean_aggregate_sinking_speed(aggs, agg_env)

            ! Limit settling velocity wrt CFL:
            ws_agg(i,j,k) = min(aggs%ws_aggregates*dtbgc, 0.99_wp*pddpo(i,j,k)) ! (m/s -> m/d)*dtb


            ! ============================== Write general diagnostics ============
            aggregate_diagnostics(i,j,k,kws_agg)    = ws_agg(i,j,k)/dtb  ! applied ws conversion  m/time_step  to  m/d for output
            aggregate_diagnostics(i,j,k,kdynvis)    = agg_env%mu     ! dynamic molecular viscosity
            aggregate_diagnostics(i,j,k,kLmax_agg)  = aggs%dmax_agg      ! applied max. diameter
            aggregate_diagnostics(i,j,k,kav_dp)     = aggs%av_dp         ! mean primary particle diameter
            aggregate_diagnostics(i,j,k,kav_rho_p)  = aggs%av_rho_p      ! mean primary particle density
            aggregate_diagnostics(i,j,k,kdf_agg)    = aggs%df_agg        ! aggregate fractal dim
            aggregate_diagnostics(i,j,k,kb_agg)     = aggs%b_agg         ! aggre number distr. slope
            aggregate_diagnostics(i,j,k,kav_d_C)    = conc_weighted_mean_agg_diameter(aggs) ! conc-weighted mean agg. diameter
            aggregate_diagnostics(i,j,k,kav_rhof_V) = volweighted_agg_density(aggs,agg_env) ! volume-weighted aggregate density
            aggregate_diagnostics(i,j,k,kav_por_V)  = volweighted_agg_porosity(aggs)        ! volume-weighted aggregate porosity
            aggregate_diagnostics(i,j,k,kstickiness_agg)      = aggs%stickiness_agg         ! aggre. stickiness
            aggregate_diagnostics(i,j,k,kstickiness_frustule) = aggs%stickiness_frustule    ! frustle stickiness
          endif
        enddo
      enddo
    enddo
  end subroutine ihamocc_mean_aggregate_sinking_speed

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
    aggregate_diagnostics = 0._wp
    m4ago_ppo             = 0._wp
    ws_agg                = 0._wp
    dyn_vis               = 0._wp

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
          if(omask(i,j) > 0.5_wp .and. pddpo(i,j,k) > dp_min) then
            m4ago_ppo(i,j,k) = 1e5_wp * ptiestu(i,j,k)*98060._wp*1.027e-6_wp ! pressure in unit Pa, 98060 = onem
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
    real, intent(in) :: psao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< salinity [psu.].
    real, intent(in) :: ppo(kpie,kpje,kpke)  !< pressure [Pa].

    ! Local variables
    real(wp) :: press_val  ! Pascal/rho -> dbar
    real(wp) :: ptho_val,psao_val
    integer  :: i,j,k,kch

    kch = 0
    !$OMP PARALLEL DO PRIVATE(i,j,k,press_val,ptho_val,psao_val,kch)
    do j = 1,kpje
      do i = 1,kpie
        do k = 1,kpke
          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5_wp) then
            kch = merge(k+1,k,k<kpke)
            if(pddpo(i,j,kch) > 0.5_wp) then
              press_val    = 0.5_wp*(ppo(i,j,k)  + ppo(i,j,kch))*1.e-5_wp ! Pascal -> dbar
              ptho_val     = 0.5_wp*(ptho(i,j,k) + ptho(i,j,kch))
              psao_val     = 0.5_wp*(psao(i,j,k) + ptho(i,j,kch))
            else
              press_val    = ppo(i,j,k)*1.e-5_wp ! Pascal -> dbar
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
