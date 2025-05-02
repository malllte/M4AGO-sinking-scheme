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

  ! CLIMBER-X ocean biogeochemistry model-specific routines:
  use control,        only: wp
  use bgc_grid,         only: kpke => kt
  use bgc_params,       only: dtb, dtbgc
  use bgc_params,       only: calcdens, claydens, opaldens, calcwei, opalwei, ropal
  !use bgc_def,          only: carb_t
  !use bgc_params,       only: iopal, ifdust, icalc, idet

  ! M4AGO routines:
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

  ! Public subroutines / called from subroutine ocprod in CLIMBER-X
  public :: climberx_mean_aggregate_sinking_speed !, init_m4ago_nml_params, init_m4ago_params
                                                  ! TODO:alloc/cleanup_mem_m4ago?


  ! Public fields and parameters
  !public :: ws_agg

  ! Fields
  !real(wp),allocatable :: ws_agg(:,:,:)       ! mass concentration-weighted aggregate mean sinking velocity
  !real(wp),allocatable :: dyn_vis(:,:,:)      ! molecular dynamic viscosity
  !real(wp),allocatable :: m4ago_ppo(:,:,:)    ! pressure


contains


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
        call mean_aggregate_sinking_speed(aggs, agg_env)
 
        ! Limit settling velocity wrt CFL:
        wsagg(k) = min(aggs%ws_aggregates*dtbgc, 0.99_wp*pddpo(k)) ! (m/s -> m/d)*dtb
                                                               ! dtbgc: time step length in seconds
                                                               ! wsagg [m/timestep]

        !tra%visco(k) = mol_dyn_vis(layer_depth(k), ptho(k), psao(k)) ! visco(pressure [dbar], temp [C], salinity [psu]) from m4ago_core

        if (wsagg(k) .lt. 0.0 .OR. wsagg(k)/dtb .gt. 530.0) then
          WRITE(*,*)'+++ WEIRD: aggregates rising / sinking faster than 530 m/day; aggregate speed [m/day] = ', wsagg(k)/dtb,' at k = ', k
          WRITE(*,*) 'CC_det(k), CC_opal(k), CC_calc(k), CC_dust(k) = ', CC_det(k), CC_opal(k), CC_calc(k), CC_dust(k) 
          WRITE(*,*) 'dtb, dtbgc = ', dtb, dtbgc
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


end module mo_m4ago_climberx
