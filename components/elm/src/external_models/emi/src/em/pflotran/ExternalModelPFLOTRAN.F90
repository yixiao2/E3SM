module ExternalModelPFLOTRANMod
#ifdef USE_PETSC_LIB

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscviewer.h"

  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use EMI_DataMod                  , only : emi_data_list, emi_data
  use elm_varctl                   , only : iulog
  use ExternalModelBaseType        , only : em_base_type
  use MultiPhysicsProbVSFM         , only : mpp_vsfm_type
  use decompMod                    , only : bounds_type
  use pflotran_model_module
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  use elm_pflotran_interface_data
  !
  implicit none
  !

  type, public, extends(em_base_type) :: em_pflotran_type
     ! ----------------------------------------------------------------------
     ! Indicies required during the initialization
     ! ----------------------------------------------------------------------
     integer :: index_l2e_init_col_active
     integer :: index_l2e_init_col_type
     integer :: index_l2e_init_col_landunit_index
     integer :: index_l2e_init_col_gridcell_index
     integer :: index_l2e_init_col_zi
     integer :: index_l2e_init_col_dz
     integer :: index_l2e_init_col_z
     integer :: index_l2e_init_col_area

     integer :: index_l2e_init_state_wtd
     integer :: index_l2e_init_state_soilp

     integer :: index_l2e_init_h2osoi_liq
     integer :: index_l2e_init_h2osoi_ice

     integer :: index_e2l_init_state_h2osoi_liq
     integer :: index_e2l_init_state_h2osoi_ice
     ! integer :: index_e2l_init_state_h2osoi_vol
     integer :: index_e2l_init_state_wtd
     ! integer :: index_e2l_init_parameter_watsatc
     ! integer :: index_e2l_init_parameter_hksatc
     ! integer :: index_e2l_init_parameter_bswc
     ! integer :: index_e2l_init_parameter_sucsatc

     integer :: index_e2l_init_flux_mflx_snowlyr_col
     integer :: index_l2e_init_flux_mflx_snowlyr_col

     integer :: index_l2e_init_landunit_type
     integer :: index_l2e_init_landunit_lakepoint
     integer :: index_l2e_init_landunit_urbanpoint

     integer :: index_l2e_init_parameter_watsatc
     integer :: index_l2e_init_parameter_hksatc
     integer :: index_l2e_init_parameter_bswc
     integer :: index_l2e_init_parameter_sucsatc
     integer :: index_l2e_init_parameter_effporosityc

     ! ----------------------------------------------------------------------
     ! Indicies required during timestepping
     ! ----------------------------------------------------------------------

     integer :: index_l2e_state_tsoil
     integer :: index_l2e_state_h2osoi_liq
     integer :: index_l2e_state_h2osoi_ice

     integer :: index_e2l_state_h2osoi_liq
     integer :: index_e2l_state_h2osoi_ice
     integer :: index_e2l_state_wtd
     integer :: index_e2l_state_soilp

     integer :: index_l2e_flux_infil
     integer :: index_l2e_flux_et
     integer :: index_l2e_flux_dew
     integer :: index_l2e_flux_snow_sub
     integer :: index_l2e_flux_snowlyr
     integer :: index_l2e_flux_drainage

     integer :: index_e2l_flux_qrecharge
     ! integer :: index_e2l_flux_drain_perched
     ! integer :: index_e2l_flux_drain
     ! integer :: index_e2l_flux_qrgwl
     ! integer :: index_e2l_flux_rsub_sat

     integer :: index_l2e_filter_hydrologyc
     integer :: index_l2e_filter_num_hydrologyc

     integer :: index_l2e_column_zi
     integer :: index_l2e_col_active
     integer :: index_l2e_col_gridcell
     integer :: index_l2e_col_dz

     type(pflotran_model_type), pointer :: pflotran_m

   contains

     procedure, public :: Populate_L2E_Init_List  => EM_PFLOTRAN_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EM_PFLOTRAN_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EM_PFLOTRAN_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_PFLOTRAN_Populate_E2L_List
     procedure, public :: PreInit                 => EM_PFLOTRAN_PreInit
     procedure, public :: Init                    => EM_PFLOTRAN_Init
     procedure, public :: Solve                   => EM_PFLOTRAN_Solve
     procedure, public :: Finalize                => EM_PFLOTRAN_Finalize

  end type em_pflotran_type

contains

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PFLOTRAN from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)             :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                         = L2E_STATE_WTD
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_wtd              = index

    id                                         = L2E_STATE_VSFM_PROGNOSTIC_SOILP
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_soilp            = index

    id                                         = L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_flux_mflx_snowlyr_col  = index

    id                                         = L2E_COLUMN_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_active             = index

    id                                         = L2E_COLUMN_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_type               = index

    id                                         = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_landunit_index     = index

    id                                         = L2E_COLUMN_GRIDCELL_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_gridcell_index     = index

    id                                         = L2E_COLUMN_ZI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_zi                 = index

    id                                         = L2E_COLUMN_DZ
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_dz                 = index

    id                                         = L2E_COLUMN_Z
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_z                  = index

    id                                         = L2E_COLUMN_AREA
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_area               = index

    id                                         = L2E_LANDUNIT_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_type          = index

    id                                         = L2E_LANDUNIT_LAKEPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_lakepoint     = index

    id                                         = L2E_LANDUNIT_URBANPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_urbanpoint    = index

    id                                         = L2E_PARAMETER_WATSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_watsatc      = index

    id                                         = L2E_PARAMETER_HKSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_hksatc       = index

    id                                         = L2E_PARAMETER_BSWC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_bswc         = index

    id                                         = L2E_PARAMETER_SUCSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_sucsatc      = index

    id                                         = L2E_PARAMETER_EFFPOROSITYC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_effporosityc = index

    id                                        = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_liq            = index

    id                                        = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_ice            = index

    deallocate(em_stages)

  end subroutine EM_PFLOTRAN_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Populate_E2L_Init_List(this, e2l_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PFLOTRAN from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)             :: this
    class(emi_data_list), intent(inout) :: e2l_init_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                        = E2L_STATE_H2OSOI_LIQ
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_liq      = index

    id                                        = E2L_STATE_H2OSOI_ICE
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_ice      = index

   !  id                                        = E2L_STATE_H2OSOI_VOL_NLEVGRND
   !  call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
   !  this%index_e2l_init_state_h2osoi_vol      = index

    id                                        = E2L_STATE_WTD
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_wtd             = index

    id                                        = E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_flux_mflx_snowlyr_col = index

   !  id                                         = E2L_PARAMETER_WATSATC
   !  call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
   !  this%index_e2l_init_parameter_watsatc      = index

   !  id                                         = E2L_PARAMETER_HKSATC
   !  call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
   !  this%index_e2l_init_parameter_hksatc       = index

   !  id                                         = E2L_PARAMETER_BSWC
   !  call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
   !  this%index_e2l_init_parameter_bswc         = index

   !  id                                         = E2L_PARAMETER_SUCSATC
   !  call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
   !  this%index_e2l_init_parameter_sucsatc      = index

    deallocate(em_stages)

  end subroutine EM_PFLOTRAN_Populate_E2L_Init_List

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PFLOTRAN from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)             :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    !
    ! !LOCAL VARIABLES:
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_PFLOTRAN_SOIL_HYDRO_STAGE

    id                                   = L2E_STATE_TSOIL_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_tsoil           = index

    id                                   = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq      = index

    id                                   = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice      = index

    id                                   = L2E_FLUX_INFIL_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_infil            = index

    id                                   = L2E_FLUX_VERTICAL_ET_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_et               = index

    id                                   = L2E_FLUX_DEW_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_dew              = index

    id                                   = L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snow_sub         = index

    id                                   = L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snowlyr          = index

    id                                   = L2E_FLUX_DRAINAGE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_drainage         = index

    id                                   = L2E_FILTER_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_hydrologyc     = index

    id                                   = L2E_FILTER_NUM_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num_hydrologyc = index

    id                                   = L2E_COLUMN_ZI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_column_zi             = index

    id                                   = L2E_COLUMN_ACTIVE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_active            = index

    id                                   = L2E_COLUMN_GRIDCELL_INDEX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_gridcell          = index

    id                                    = L2E_COLUMN_DZ
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_dz                 = index

    deallocate(em_stages)

  end subroutine EM_PFLOTRAN_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Populate_E2L_List(this, e2l_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by PFLOTRAN from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)             :: this
    class(emi_data_list), intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_PFLOTRAN_SOIL_HYDRO_STAGE

    id                              = E2L_STATE_H2OSOI_LIQ
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_liq = index

    id                              = E2L_STATE_H2OSOI_ICE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_ice = index

    id                              = E2L_STATE_WTD
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_wtd        = index

    id                              = E2L_STATE_VSFM_PROGNOSTIC_SOILP
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_soilp      = index

    id                              = E2L_FLUX_AQUIFER_RECHARGE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_qrecharge   = index

   !  id                              = E2L_FLUX_DRAIN_PERCHED
   !  call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
   !  this%index_e2l_flux_drain_perched   = index

   !  id                              = E2L_FLUX_DRAIN
   !  call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
   !  this%index_e2l_flux_drain       = index

   !  id                              = E2L_FLUX_QRGWL
   !  call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
   !  this%index_e2l_flux_qrgwl   = index

   !  id                              = E2L_FLUX_RSUB_SAT
   !  call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
   !  this%index_e2l_flux_rsub_sat   = index

    deallocate(em_stages)

  end subroutine EM_PFLOTRAN_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_PreInit(this, prefix, restart_stamp)
    !
    use spmdMod     , only : mpicom
    !
    implicit none
    !
    class(em_pflotran_type)      :: this
    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: restart_stamp
    
    this%pflotran_m => pflotranModelCreate(mpicom, prefix)
    call pflotranModelSetupRestart(this%pflotran_m, restart_stamp)

  end subroutine EM_PFLOTRAN_PreInit

  !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)

    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)              :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! LOCAL VARAIBLES:
    integer :: elm_npts
    integer :: elm_surf_npts

    ! Create ELM-PFLOTRAN mapping files
    call CreateELMPFLOTRANInterfaceDate(this, bounds_clump, elm_npts, elm_surf_npts)
#ifdef DEBUG_ELMPFEH
     write(*,*) '[YX DEBUG][ExternalModelPFLOTRAN::EMPFLOTRAN_Init] CreateELMPFLOTRANInterfaceDate done'
#endif
    ! Create ELM-PFLOTRAN mapping files
    call CreateELMPFLOTRANMaps(this, bounds_clump, elm_npts, elm_surf_npts)
#ifdef DEBUG_ELMPFEH
     write(*,*) '[YX DEBUG][ExternalModelPFLOTRAN::EMPFLOTRAN_Init] CreateELMPFLOTRANMaps done'
#endif
    ! Initialize PFLOTRAN states
    call pflotranModelStepperRunInit(this%pflotran_m)
#ifdef DEBUG_ELMPFEH
     write(*,*) '[YX DEBUG][ExternalModelPFLOTRAN::EMPFLOTRAN_Init] pass pflotranModelStepperRunInit done'
#endif
    ! Get top surface area
    call pflotranModelGetTopFaceArea(this%pflotran_m)
#ifdef DEBUG_ELMPFEH
     write(*,*) '[YX DEBUG][ExternalModelPFLOTRAN::EMPFLOTRAN_Init] pass pflotranModelGetTopFaceArea done'
#endif
    ! Get PFLOTRAN states
    call pflotranModelGetUpdatedData(this%pflotran_m)
#ifdef DEBUG_ELMPFEH
     write(*,*) '[YX DEBUG][ExternalModelPFLOTRAN::EMPFLOTRAN_Init] pass pflotranModelGetUpdatedData done'
#endif
    ! Save the data need by ELM
    call extract_data_for_elm(this, l2e_init_list, e2l_init_list, bounds_clump)
#ifdef DEBUG_ELMPFEH
     write(*,*) '[YX DEBUG][ExternalModelPFLOTRAN::EMPFLOTRAN_Init] pass extract_data_for_elm done'
#endif
  end subroutine EM_PFLOTRAN_Init

  !-----------------------------------------------------------------------
  subroutine CreateELMPFLOTRANInterfaceDate(this, bounds, elm_npts, elm_surf_npts)
    !
    ! !DESCRIPTION:
    ! Allocates memory for ELM-PFLOTRAN interface
    !
    ! !USES:
    use decompMod                     , only : bounds_type
    use elm_varpar                    , only : nlevsoi, nlevgrnd
    use shr_kind_mod                  , only: r8 => shr_kind_r8
    ! pflotran
    use Option_module                 , only : printErrMsg
    use Simulation_Base_class         , only : simulation_base_type
    use Simulation_Subsurface_class   , only : simulation_subsurface_type
    ! use Simulation_Surface_class      , only : simulation_surface_type
    ! use Simulation_Surf_Subsurf_class , only : simulation_surfsubsurface_type
    use Realization_Subsurface_class  , only : realization_subsurface_type
    ! use Realization_Surface_class     , only : realization_surface_type
    use PFLOTRAN_Constants_module
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(inout) :: elm_npts
    integer                 , intent(inout) :: elm_surf_npts
    !
    ! LOCAL VARAIBLES:
    integer                                      :: nlevmapped
    character(len= 128)                          :: subname = 'CreateELMPFLOTRANInterfaceDate' ! subroutine name
    class(realization_subsurface_type) , pointer :: realization
    ! class(realization_surface_type)    , pointer :: surf_realization

    ! Initialize PETSc vector for data transfer between ELM and PFLOTRAN
    call ELMPFLOTRANIDataInit()

    nullify(realization)
    ! nullify(surf_realization)

    select type (simulation => this%pflotran_m%simulation)
    class is (simulation_subsurface_type)
       realization => simulation%realization
       ! nullify(surf_realization)
   !  class is (simulation_surface_type)
   !     nullify(realization)
   !     surf_realization => simulation%surf_realization
   !  class is (simulation_surfsubsurface_type)
   !     realization => simulation%realization
   !     surf_realization => simulation%surf_realization
    class default
       this%pflotran_m%option%io_buffer = "elm-pflotran only works with surface and subsurface simulations."
       write(*, '(/A/)') this%pflotran_m%option%io_buffer
       call printErrMsg(this%pflotran_m%option)
    end select

    ! Compute number of cells in ELM domain.
    ! Assumption-1: One column per ELM grid cell.

    ! Check if the number of ELM vertical soil layers defined in the mapping
    ! file read by PFLOTRAN matches either nlevsoi or nlevgrnd
    elm_pf_idata%nzelm_mapped = this%pflotran_m%map_elm_sub_to_pf_sub%elm_nlevsoi
    nlevmapped                = elm_pf_idata%nzelm_mapped
    if ( (nlevmapped /= nlevsoi) .and. (nlevmapped /= nlevgrnd) ) then
       call endrun(trim(subname)//' ERROR: Number of layers PFLOTRAN thinks ELM should '// &
            'have do not match either nlevsoi or nlevgrnd. Abortting' )
    end if

    elm_npts = (bounds%endg - bounds%begg + 1)*nlevmapped
    elm_surf_npts = (bounds%endg - bounds%begg + 1)

    ! ELM: Subsurface domain (local and ghosted cells)
    elm_pf_idata%nlelm_sub = elm_npts
    elm_pf_idata%ngelm_sub = elm_npts

    ! ELM: Surface of subsurface domain (local and ghosted cells)
    elm_pf_idata%nlelm_2dsub = (bounds%endg - bounds%begg + 1)
    elm_pf_idata%ngelm_2dsub = (bounds%endg - bounds%begg + 1)
    ! For ELM: Same as surface of subsurface domain
    elm_pf_idata%nlelm_srf = elm_surf_npts
    elm_pf_idata%ngelm_srf = elm_surf_npts

    ! PFLOTRAN: Subsurface domain (local and ghosted cells)
    elm_pf_idata%nlpf_sub = realization%patch%grid%nlmax
    elm_pf_idata%ngpf_sub = realization%patch%grid%ngmax

    ! PFLOTRAN: Surface of subsurface domain (local and ghosted cells)
    if (this%pflotran_m%option%iflowmode == TH_MODE) then
       elm_pf_idata%nlpf_2dsub = pflotranModelNSurfCells3DDomain(this%pflotran_m)
       elm_pf_idata%ngpf_2dsub = pflotranModelNSurfCells3DDomain(this%pflotran_m)
    else
       elm_pf_idata%nlpf_2dsub = 0
       elm_pf_idata%ngpf_2dsub = 0
    endif

    ! PFLOTRAN: Surface domain (local and ghosted cells)
   !  if (associated(surf_realization) .and. this%pflotran_m%option%nsurfflowdof > 0) then
   !     elm_pf_idata%nlpf_srf = surf_realization%patch%grid%nlmax
   !     elm_pf_idata%ngpf_srf = surf_realization%patch%grid%ngmax
   !  else
       elm_pf_idata%nlpf_srf = 0
       elm_pf_idata%ngpf_srf = 0
   !  endif

    ! Allocate vectors for data transfer between ELM and PFLOTRAN.
    call ELMPFLOTRANIDataCreateVec(MPI_COMM_WORLD)

  end subroutine CreateELMPFLOTRANInterfaceDate

  !-----------------------------------------------------------------------
  subroutine CreateELMPFLOTRANMaps(this, bounds, elm_npts, elm_surf_npts)
    !
    ! !DESCRIPTION:
    ! Creates maps to transfer data between ELM and PFLOTRAN
    !
    ! !USES:
    use elm_varctl      , only : pflotran_th_mode, pflotran_th_freezing
    use decompMod       , only : bounds_type, ldecomp
    use PFLOTRAN_Constants_module
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(inout) :: elm_npts
    integer                 , intent(inout) :: elm_surf_npts
    !
    ! LOCAL VARAIBLES:
    integer             :: g,j
    integer             :: nlevmapped
    integer, pointer    :: elm_cell_ids_nindex(:)
    integer, pointer    :: elm_surf_cell_ids_nindex(:)

    ! Save cell IDs of ELM grid
    allocate(elm_cell_ids_nindex(     1:elm_npts     ))
    allocate(elm_surf_cell_ids_nindex(1:elm_surf_npts))

    nlevmapped     = elm_pf_idata%nzelm_mapped
    elm_npts       = 0
    elm_surf_npts  = 0
    do g = bounds%begg, bounds%endg
       do j = 1,nlevmapped
          elm_npts = elm_npts + 1
          elm_cell_ids_nindex(elm_npts) = (ldecomp%gdc2glo(g)-1)*nlevmapped + j - 1
       enddo
       elm_surf_npts = elm_surf_npts + 1
       elm_surf_cell_ids_nindex(elm_surf_npts) = (ldecomp%gdc2glo(g)-1)*nlevmapped
    enddo

    ! Initialize maps for transferring data between ELM and PFLOTRAN.
    call pflotranModelInitMapping(this%pflotran_m, elm_cell_ids_nindex, &
                                  elm_npts, ELM_SUB_TO_PF_SUB)
    call pflotranModelInitMapping(this%pflotran_m, elm_cell_ids_nindex, &
                                  elm_npts, ELM_SUB_TO_PF_EXTENDED_SUB)
    call pflotranModelInitMapping(this%pflotran_m, elm_cell_ids_nindex, &
                                  elm_npts, PF_SUB_TO_ELM_SUB)

    if (this%pflotran_m%option%iflowmode == TH_MODE) then
      pflotran_th_mode = .true.
      if (this%pflotran_m%option%flow%th_freezing) pflotran_th_freezing  = .true.
    endif

   !  if (this%pflotran_m%option%nsurfflowdof > 0) then
   !    pflotran_surfaceflow = .true.
   !    call pflotranModelInitMapping(this%pflotran_m, elm_surf_cell_ids_nindex, &
   !                                  elm_surf_npts, PF_SRF_TO_ELM_SRF)
   !    call pflotranModelInitMapping(this%pflotran_m, elm_surf_cell_ids_nindex, &
   !                                  elm_surf_npts, ELM_SRF_TO_PF_SRF)
   !  else
      if (this%pflotran_m%option%iflowmode == TH_MODE) then
        call pflotranModelInitMapping(this%pflotran_m, elm_surf_cell_ids_nindex, &
                                      elm_surf_npts, ELM_SRF_TO_PF_2DSUB)
      endif
   !  endif

    deallocate(elm_cell_ids_nindex)
    deallocate(elm_surf_cell_ids_nindex)

  end subroutine CreateELMPFLOTRANMaps


  !-----------------------------------------------------------------------
  subroutine extract_data_for_elm(this, l2e_init_list, e2l_init_list, bounds_clump)
    !
    !DESCRIPTION
    !  Saves
    !
    use landunit_varcon           , only : istsoil
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use elm_varcon                , only : denice, denh2o
    use elm_varpar                , only : nlevgrnd
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_pflotran_type)              :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! !LOCAL VARIABLES:
    integer               :: c,j,g,l,pf_j  ! do loop indices
    integer               :: nlevmapped
    integer               :: gcount
    integer               :: bounds_proc_begc, bounds_proc_endc

    integer     , pointer :: col_active(:)
    integer     , pointer :: col_landunit(:)
    integer     , pointer :: col_gridcell(:)
    integer     , pointer :: lun_type(:)

    real(r8)    , pointer :: l2e_h2osoi_ice(:,:)

    real(r8)    , pointer :: e2l_h2osoi_liq(:,:)
    real(r8)    , pointer :: e2l_h2osoi_ice(:,:)
    ! real(r8)    , pointer :: e2l_h2osoi_vol(:,:)
    real(r8)    , pointer :: e2l_zwt(:)
    real(r8)    , pointer :: e2l_mflx_snowlyr_col(:)
    real(r8)    , pointer :: e2l_watsatc(:,:)
    ! real(r8)    , pointer :: e2l_hksatc(:,:)
    ! real(r8)    , pointer :: e2l_bswc(:,:)
    ! real(r8)    , pointer :: e2l_sucsatc(:,:)

    real(r8)    , pointer :: dz(:,:)

    PetscScalar , pointer :: sat_elm_loc(:)
    PetscScalar , pointer :: watsat_elm_loc(:)
    PetscScalar , pointer :: hksat_elm_loc(:)
    PetscScalar , pointer :: bsw_elm_loc(:)
    PetscScalar , pointer :: sucsat_elm_loc(:)
    PetscErrorCode        :: ierr

    character(len= 128)   :: subname = 'extract_data_for_elm'
    !-----------------------------------------------------------------------

    bounds_proc_begc = bounds_clump%begc
    bounds_proc_endc = bounds_clump%endc
    nlevmapped       = elm_pf_idata%nzelm_mapped

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active             , col_active   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index     , col_landunit )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_gridcell_index     , col_gridcell )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_dz                , dz           )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_type          , lun_type     )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_h2osoi_ice            , l2e_h2osoi_ice       )

    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_state_wtd             , e2l_zwt              )
    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_flux_mflx_snowlyr_col , e2l_mflx_snowlyr_col )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_liq      , e2l_h2osoi_liq       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_ice      , e2l_h2osoi_ice       )
    ! call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_vol      , e2l_h2osoi_vol       )

    ! call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_watsatc     , e2l_watsatc )
    ! call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_hksatc      , e2l_hksatc  )
    ! call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_bswc        , e2l_bswc    )
    ! call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_sucsatc     , e2l_sucsatc )

    call pflotranModelGetSoilProp(this%pflotran_m)

    ! Set initial value of for ELM
    e2l_mflx_snowlyr_col(:) = 0._r8
    e2l_zwt(:)              = 0._r8

    ! Initialize soil moisture
    call VecGetArrayF90(elm_pf_idata%sat_elms      , sat_elm_loc    , ierr)
    call VecGetArrayF90(elm_pf_idata%watsat2_elm  , watsat_elm_loc , ierr)
    call VecGetArrayF90(elm_pf_idata%hksat_x2_elm , hksat_elm_loc  , ierr)
    call VecGetArrayF90(elm_pf_idata%bsw2_elm     , bsw_elm_loc    , ierr)
    call VecGetArrayF90(elm_pf_idata%sucsat2_elm  , sucsat_elm_loc , ierr)

    do c = bounds_proc_begc, bounds_proc_endc
       if (col_active(c) == 1) then
          l = col_landunit(c)
          if (lun_type(l) == istsoil) then
             g = col_gridcell(c)
             gcount = g - bounds_clump%begg
             do j = 1, nlevgrnd
                pf_j = gcount*nlevmapped + j

                if (j <= nlevmapped) then
                   e2l_h2osoi_liq(c,j) = sat_elm_loc(pf_j)*watsat_elm_loc(pf_j)*dz(c,j)*1.e3_r8

                   ! e2l_h2osoi_vol(c,j) = e2l_h2osoi_liq(c,j)/dz(c,j)/denh2o + &
                   !     l2e_h2osoi_ice(c,j)/dz(c,j)/denice
                   ! e2l_h2osoi_vol(c,j) = min(e2l_h2osoi_vol(c,j),watsat_elm_loc(pf_j))
                   e2l_h2osoi_ice(c,j) = 0._r8

                   ! e2l_watsatc(c,j) = watsat_elm_loc(pf_j)
                   ! e2l_hksatc(c,j)  = hksat_elm_loc(pf_j)
                   ! e2l_bswc(c,j)    = bsw_elm_loc(pf_j)
                   ! e2l_sucsatc(c,j) = sucsat_elm_loc(pf_j)
                else
                   e2l_h2osoi_liq(c,j) = e2l_h2osoi_liq(c,nlevmapped)
                   ! e2l_h2osoi_vol(c,j) = e2l_h2osoi_vol(c,nlevmapped)
                   e2l_h2osoi_ice(c,j) = 0._r8
                   ! e2l_watsatc(c,j)    = e2l_watsatc(c,nlevmapped)
                   ! e2l_hksatc(c,j)     = e2l_hksatc(c,nlevmapped)
                   ! e2l_bswc(c,j)       = e2l_bswc(c,nlevmapped)
                   ! e2l_sucsatc(c,j)    = e2l_sucsatc(c,nlevmapped)
                end if
! #ifdef DEBUG_ELMPFEH
!    write(*,*) '[YX DEBUG][ExternalModelPFLOTRAN::extract_data_for_elm] nlevmapped=',nlevmapped
!    write(*,*) ' c=',c
!    write(*,*) ' j=',j
!    write(*,*) ' e2l_h2osoi_liq(c,j)=',e2l_h2osoi_liq(c,j)
! #endif
             enddo
          else
             write(iulog,*)'WARNING: Land Unit type other than soil type is present within the domain'
             call endrun( trim(subname)//' ERROR: Land Unit type not supported' )             
          endif
       endif
    enddo

    call VecRestoreArrayF90(elm_pf_idata%sat_elms      , sat_elm_loc    , ierr)
    call VecRestoreArrayF90(elm_pf_idata%watsat2_elm  , watsat_elm_loc , ierr)
    call VecRestoreArrayF90(elm_pf_idata%hksat_x2_elm , hksat_elm_loc  , ierr)
    call VecRestoreArrayF90(elm_pf_idata%bsw2_elm     , bsw_elm_loc    , ierr)
    call VecRestoreArrayF90(elm_pf_idata%sucsat2_elm  , sucsat_elm_loc , ierr)

   end subroutine extract_data_for_elm

   !------------------------------------------------------------------------
   subroutine EM_PFLOTRAN_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)
    !
    ! !DESCRIPTION:
    ! The VSFM dirver subroutine
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)              :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    select case (em_stage)
    case (EM_PFLOTRAN_SOIL_HYDRO_STAGE)
       call EM_PFLOTRAN_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, &
            bounds_clump)

    case default
       write(iulog,*)'EM_PFLOTRAN_Solve: Unknown em_stage.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine EM_PFLOTRAN_Solve


    !------------------------------------------------------------------------
  subroutine EM_PFLOTRAN_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    ! Solve the Variably Saturated Flow Model (VSFM) in soil columns.
    !
#include <petsc/finclude/petsc.h>
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_FRAC_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use MultiPhysicsProbConstants , only : VAR_LATERAL_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : VAR_BC_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use mpp_varpar                , only : nlevgrnd
    use elm_varcon                , only : denice, denh2o
    use petscsnes
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_pflotran_type)              :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer                              :: p,c,fc,j,g                                                       ! do loop indices
    integer                              :: pi                                                               ! pft index
    real(r8)                             :: dzsum                                                            ! summation of dzmm of layers below water table (mm)
    real(r8)                             :: dtime

    real(r8)  , pointer                  :: mflx_et_col_1d         (:)
    real(r8)  , pointer                  :: mflx_infl_col_1d       (:)
    real(r8)  , pointer                  :: mflx_dew_col_1d        (:)
    real(r8)  , pointer                  :: mflx_drain_col_1d      (:)
    real(r8)  , pointer                  :: mflx_sub_snow_col_1d   (:)
    real(r8)  , pointer                  :: mflx_snowlyr_col_1d    (:)
    real(r8)  , pointer                  :: t_soil_col_1d          (:)
    integer    , pointer                 :: col_active             (:)
    integer    , pointer                 :: col_gridcell           (:)

    real(r8)  , pointer                  :: fliq_col_1d       (:)
    real(r8)  , pointer                  :: mass_col_1d       (:)
    real(r8)  , pointer                  :: smpl_col_1d       (:)
    real(r8)  , pointer                  :: soilp_col_1d      (:)
    real(r8)  , pointer                  :: sat_col_1d        (:)

    real(r8)  , pointer                  :: frac_ice                    (:,:) ! fraction of ice
    real(r8)  , pointer                  :: total_mass_flux_col         (:)            ! Sum of all source-sinks conditions for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_et_col      (:)            ! ET sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_infl_col    (:)            ! Infiltration source for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_dew_col     (:)            ! Dew source for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_drain_col   (:)            ! Drainage sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_snowlyr_col (:)            ! Flux due to disappearance of snow for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_sub_col     (:)            ! Sublimation sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_lateral_col (:)            ! Lateral flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_seepage_col (:)            ! Seepage flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: qflx_seepage                (:)            ! Seepage flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: mass_prev_col          (:,:) ! Mass of water before a VSFM solve
    real(r8)  , pointer                  :: dmass_col              (:)            ! Change in mass of water after a VSFM solve
    real(r8)  , pointer                  :: mass_beg_col                (:)            ! Total mass before a VSFM solve
    real(r8)  , pointer                  :: mass_end_col                (:)            ! Total mass after a VSFM solve
    real(r8)  , pointer                  :: mass_bal_error_col          (:)            ! Mass balance error for a VSFM solve
    integer                              :: ier                                                              ! error status

    integer                              :: begc, endc
    integer                              :: g_idx, c_idx

    PetscInt                             :: soe_auxvar_id                                                    ! Index of system-of-equation's (SoE's) auxvar
    PetscErrorCode                       :: ierr                                                             ! PETSc return error code

    PetscBool                            :: converged                                                        ! Did VSFM solver converge to a solution with given PETSc SNES tolerances
    PetscInt                             :: converged_reason                                                 ! SNES converged due to which criteria
    PetscReal                            :: atol_default                                                     ! Default SNES absolute convergance tolerance
    PetscReal                            :: rtol_default                                                     ! Default SNES relative convergance tolerance
    PetscReal                            :: stol_default                                                     ! Default SNES solution convergance tolerance
    PetscInt                             :: max_it_default                                                   ! Default SNES maximum number of iteration
    PetscInt                             :: max_f_default                                                    ! Default SNES maximum number of function evaluation
    PetscReal                            :: stol                                                             ! solution convergance tolerance
    PetscReal                            :: rtol                                                             ! relative convergance tolerance
    PetscReal,parameter                  :: stol_alternate = 1.d-10                                          ! Alternate solution convergance tolerance

    PetscReal                            :: mass_beg                                                         ! Sum of mass of water for all active soil columns before VSFM is called
    PetscReal                            :: mass_end                                                         ! Sum of mass of water for all active soil columns after VSFM is called
    PetscReal                            :: total_mass_flux_et                                               ! Sum of mass ET mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_infl                                             ! Sum of mass infiltration mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_dew                                              ! Sum of mass dew mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_drain                                            ! Sum of mass drainage mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_snowlyr                                          ! Sum of mass snow layer disappearance mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_sub                                              ! Sum of mass sublimation mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_lateral                                          ! Sum of lateral mass flux for all active soil columns
    PetscReal                            :: total_mass_flux                                                  ! Sum of mass ALL mass flux of water for all active soil columns
    PetscInt                             :: iter_count                                                       ! How many times VSFM solver is called

    PetscInt, parameter                  :: max_iter_count = 10                                              ! Maximum number of times VSFM can be called
    PetscInt                             :: diverged_count                                                   ! Number of time VSFM solver diverged
    PetscInt                             :: mass_bal_err_count                                               ! Number of time VSFM solver returns a solution that isn't within acceptable mass balance error threshold
    PetscReal                            :: abs_mass_error_col                                               ! Maximum absolute error for any active soil column
    PetscReal, parameter                 :: max_abs_mass_error_col  = 1.e-5                                  ! Acceptable mass balance error
    PetscReal                            :: total_mass_bal_error                                             ! Sum of mass balance error for all active soil columns
    PetscBool                            :: successful_step                                                  ! Is the solution return by VSFM acceptable
    PetscReal , pointer                  :: soilp_col_ghosted_1d(:)
    PetscReal , pointer                  :: fliq_col_ghosted_1d(:)
    PetscReal , pointer                  :: mflx_lateral_col_1d(:)
    PetscReal , pointer                  :: lat_mass_exc_col_1d(:)
    PetscReal , pointer                  :: seepage_mass_exc_col_1d(:)
    PetscReal , pointer                  :: seepage_press_1d(:)

    integer                              :: jwt
    real(r8)                             :: z_dn, z_up

    real(r8)  , pointer                  :: l2e_mflux_infil(:)
    real(r8)  , pointer                  :: l2e_mflux_dew(:)
    real(r8)  , pointer                  :: l2e_mflux_sub_snow(:)
    real(r8)  , pointer                  :: l2e_mflux_snowlyr(:)
    real(r8)  , pointer                  :: l2e_mflux_et(:,:)
    real(r8)  , pointer                  :: l2e_mflux_drain(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: l2e_zi(:,:)
    real(r8)  , pointer                  :: col_dz(:,:)
    integer   , pointer                  :: l2e_filter_hydrologyc(:)
    integer                              :: l2e_num_hydrologyc
 
    real(r8)  , pointer                  :: e2l_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: e2l_smp(:,:)
    real(r8)  , pointer                  :: e2l_wtd(:)
    real(r8)  , pointer                  :: e2l_soilp(:,:)
    real(r8)  , pointer                  :: e2l_qrecharge(:)

    PetscScalar, pointer :: qflx_elm_loc(:)
    PetscScalar, pointer :: area_elm_loc(:)
    PetscScalar, pointer :: thetares2_elm_loc(:)
    PetscScalar, pointer :: watsat_elm_loc(:)
    PetscScalar, pointer :: sat_elm_loc(:)
    PetscScalar, pointer :: mass_elm_loc(:)
    ! PetscScalar, pointer :: e2l_drain_perched(:)
    ! PetscScalar, pointer :: e2l_drain(:)
    ! PetscScalar, pointer :: e2l_qrgwl(:)
    ! PetscScalar, pointer :: e2l_rsub_sat(:)

    PetscViewer :: viewer

    integer :: bounds_proc_begc, bounds_proc_endc
    integer :: nlevmapped
    real(r8) :: col_wtgcell
    !-----------------------------------------------------------------------

    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    ! Get time step

    dtime = dt

    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_infil       , l2e_mflux_infil       )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_dew         , l2e_mflux_dew         )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_snow_sub    , l2e_mflux_sub_snow    )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_snowlyr     , l2e_mflux_snowlyr     )

    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_et          , l2e_mflux_et          )
    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_drainage    , l2e_mflux_drain       )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq , l2e_h2osoi_liq        )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice , l2e_h2osoi_ice        )

    call l2e_list%GetPointerToInt1D(this%index_l2e_filter_hydrologyc , l2e_filter_hydrologyc )
    call l2e_list%GetIntValue(this%index_l2e_filter_num_hydrologyc   , l2e_num_hydrologyc    )

    call l2e_list%GetPointerToReal2D(this%index_l2e_column_zi        , l2e_zi                )
    call l2e_list%GetPointerToInt1D(this%index_l2e_col_active        , col_active            )
    call l2e_list%GetPointerToInt1D(this%index_l2e_col_gridcell      , col_gridcell          )
    call l2e_list%GetPointerToReal2D(this%index_l2e_init_col_dz      , col_dz                )

    call e2l_list%GetPointerToReal1D(this%index_e2l_state_wtd        , e2l_wtd               )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_liq , e2l_h2osoi_liq        )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_ice , e2l_h2osoi_ice        )
    !call e2l_list%GetPointerToReal2D(this%index_e2l_state_smp        , e2l_smp               )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_soilp      , e2l_soilp             )

    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_qrecharge    , e2l_qrecharge        )
    ! call e2l_list%GetPointerToReal1D(this%index_e2l_flux_drain_perched, e2l_drain_perched    )
    ! call e2l_list%GetPointerToReal1D(this%index_e2l_flux_drain        , e2l_drain            )
    ! call e2l_list%GetPointerToReal1D(this%index_e2l_flux_qrgwl        , e2l_qrgwl            )
    ! call e2l_list%GetPointerToReal1D(this%index_e2l_flux_rsub_sat     , e2l_rsub_sat         )

    begc = bounds_proc_begc
    endc = bounds_proc_endc

    allocate(frac_ice                    (begc:endc,1:nlevgrnd))
    allocate(total_mass_flux_col         (begc:endc))
    allocate(total_mass_flux_et_col      (begc:endc))
    allocate(total_mass_flux_infl_col    (begc:endc))
    allocate(total_mass_flux_dew_col     (begc:endc))
    allocate(total_mass_flux_drain_col   (begc:endc))
    allocate(total_mass_flux_snowlyr_col (begc:endc))
    allocate(total_mass_flux_sub_col     (begc:endc))
    allocate(total_mass_flux_lateral_col (begc:endc))
    allocate(total_mass_flux_seepage_col (begc:endc))
    allocate(qflx_seepage                (begc:endc))
    allocate(mass_prev_col          (begc:endc,1:nlevgrnd))
    allocate(dmass_col              (begc:endc))
    allocate(mass_beg_col                (begc:endc))
    allocate(mass_end_col                (begc:endc))
    allocate(mass_bal_error_col          (begc:endc))

    allocate(mflx_et_col_1d              ((endc-begc+1)*nlevgrnd))
    allocate(mflx_drain_col_1d           ((endc-begc+1)*nlevgrnd))
    allocate(mflx_infl_col_1d            (endc-begc+1))
    allocate(mflx_dew_col_1d             (endc-begc+1))
    allocate(mflx_sub_snow_col_1d        (endc-begc+1))
    allocate(mflx_snowlyr_col_1d         (endc-begc+1))
    allocate(t_soil_col_1d               ((endc-begc+1)*nlevgrnd))

    allocate(mass_col_1d            ((endc-begc+1)*nlevgrnd))
    allocate(fliq_col_1d            ((endc-begc+1)*nlevgrnd))
    allocate(smpl_col_1d            ((endc-begc+1)*nlevgrnd))
    allocate(soilp_col_1d           ((endc-begc+1)*nlevgrnd))
    allocate(sat_col_1d             ((endc-begc+1)*nlevgrnd))

    ! initialize
    mflx_et_col_1d(:)                = 0.d0
    mflx_infl_col_1d(:)              = 0.d0
    mflx_dew_col_1d(:)               = 0.d0
    mflx_drain_col_1d(:)             = 0.d0
    mflx_sub_snow_col_1d(:)          = 0.d0
    mflx_snowlyr_col_1d(:)           = 0.d0
    t_soil_col_1d(:)                 = 298.15d0

    mass_beg                         = 0.d0
    mass_end                         = 0.d0
    total_mass_flux                  = 0.d0
    total_mass_flux_et               = 0.d0
    total_mass_flux_infl             = 0.d0
    total_mass_flux_dew              = 0.d0
    total_mass_flux_drain            = 0.d0
    total_mass_flux_snowlyr          = 0.d0
    total_mass_flux_sub              = 0.d0
    total_mass_flux_lateral          = 0.d0

    mass_beg_col(:)                  = 0.d0
    mass_end_col(:)                  = 0.d0
    mass_bal_error_col(:)            = 0.d0
    total_mass_flux_col(:)           = 0.d0
    total_mass_flux_et_col(:)        = 0.d0
    total_mass_flux_infl_col(:)      = 0.d0
    total_mass_flux_dew_col(:)       = 0.d0
    total_mass_flux_drain_col(:)     = 0.d0
    total_mass_flux_snowlyr_col(:)   = 0.d0
    total_mass_flux_sub_col(:)       = 0.d0
    total_mass_flux_lateral_col(:)   = 0.d0
    total_mass_bal_error             = 0.d0

    mass_prev_col(:,:)          = 0.d0
    dmass_col(:)                = 0.d0

    nlevmapped = elm_pf_idata%nzelm_mapped

    ! Get total mass
    call pflotranModelGetUpdatedData( this%pflotran_m )
#ifdef DEBUG_ELMPFEH
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] pass 1st pflotranModelGetUpdatedData'
#endif

    call VecGetArrayF90(elm_pf_idata%mass_elms  , mass_elm_loc  , ierr); CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%area_top_face_elms, area_elm_loc, ierr); CHKERRQ(ierr)
#ifdef DEBUG_ELMPFEH
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] check elm_pf_idata update-1'
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] mass_elm_loc = ', mass_elm_loc
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] area_elm_loc = ', area_elm_loc
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] l2e_filter_hydrologyc = ', l2e_filter_hydrologyc
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] l2e_num_hydrologyc = ', l2e_num_hydrologyc
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] nlevmapped = ', nlevmapped
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] l2e_mflux_infil = ', l2e_mflux_infil
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] l2e_mflux_dew = ', l2e_mflux_dew
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] l2e_mflux_sub_snow = ', l2e_mflux_sub_snow
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] l2e_mflux_snowlyr = ', l2e_mflux_snowlyr
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] begc = ', begc
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] endc = ', endc
    !stop
#endif
    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)
       g = col_gridcell(c)

       do j = 1, nlevmapped

          c_idx = (c - begc)*nlevgrnd + j
          g_idx = (g - bounds_clump%begg)*nlevmapped + j

          mflx_et_col_1d(c_idx)          = l2e_mflux_et(c,j)
          mflx_drain_col_1d(c_idx)       = l2e_mflux_drain(c,j)

          total_mass_flux_et           = total_mass_flux_et           + mflx_et_col_1d(c_idx)
          total_mass_flux_et_col(c)    = total_mass_flux_et_col(c)    + mflx_et_col_1d(c_idx)

          total_mass_flux_drain        = total_mass_flux_drain        + mflx_drain_col_1d(c_idx)
          total_mass_flux_drain_col(c) = total_mass_flux_drain_col(c) + mflx_drain_col_1d(c_idx)

          mass_beg                     = mass_beg                     + mass_elm_loc(g_idx)/area_elm_loc(g_idx)
          mass_beg_col(c)              = mass_beg_col(c)              + mass_elm_loc(g_idx)/area_elm_loc(g_idx)
          mass_prev_col(c,j)      = mass_col_1d(c_idx)
       end do

       c_idx = c - begc+1

       mflx_dew_col_1d(c_idx)        = l2e_mflux_dew(c)
       mflx_infl_col_1d(c_idx)       = l2e_mflux_infil(c)
       mflx_snowlyr_col_1d(c_idx)    = l2e_mflux_snowlyr(c)
       mflx_sub_snow_col_1d(c_idx)   = l2e_mflux_sub_snow(c)

       total_mass_flux_dew            = total_mass_flux_dew            + mflx_dew_col_1d(c_idx)
       total_mass_flux_dew_col(c)     = total_mass_flux_dew_col(c)     + mflx_dew_col_1d(c_idx)

       total_mass_flux_infl           = total_mass_flux_infl           + mflx_infl_col_1d(c_idx)
       total_mass_flux_infl_col(c)    = total_mass_flux_infl_col(c)    + mflx_infl_col_1d(c_idx)

       total_mass_flux_snowlyr        = total_mass_flux_snowlyr        + mflx_snowlyr_col_1d(c_idx)
       total_mass_flux_snowlyr_col(c) = total_mass_flux_snowlyr_col(c) + mflx_snowlyr_col_1d(c_idx)

       total_mass_flux_sub            = total_mass_flux_sub            + mflx_sub_snow_col_1d(c_idx)
       total_mass_flux_sub_col(c)     = total_mass_flux_sub_col(c)     + mflx_sub_snow_col_1d(c_idx)

       total_mass_flux_col(c) = total_mass_flux_et_col(c)      + &
            total_mass_flux_infl_col(c)    + &
            total_mass_flux_dew_col(c)     + &
            total_mass_flux_drain_col(c)   + &
            total_mass_flux_snowlyr_col(c) + &
            total_mass_flux_sub_col(c)     + &
            total_mass_flux_lateral_col(c)
    end do
    total_mass_flux        = &
         total_mass_flux_et        + &
         total_mass_flux_infl      + &
         total_mass_flux_dew       + &
         total_mass_flux_drain     + &
         total_mass_flux_snowlyr   + &
         total_mass_flux_sub       + &
         total_mass_flux_lateral
    call VecRestoreArrayF90(elm_pf_idata%mass_elms  , mass_elm_loc  , ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%area_top_face_elms, area_elm_loc, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(elm_pf_idata%qflx_elm, qflx_elm_loc, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%area_top_face_elms, area_elm_loc, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%thetares2_elm, thetares2_elm_loc, ierr); CHKERRQ(ierr)
#ifdef DEBUG_ELMPFEH
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] check elm_pf_idata update0'
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] qflx_elm_loc = ', qflx_elm_loc
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] thetares2_elm_loc = ', thetares2_elm_loc
    !stop
#endif
    frac_ice(:,:)       = 0.d0
    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)
       do j = 1, nlevmapped
          frac_ice(c,j) = l2e_h2osoi_ice(c,j)/(l2e_h2osoi_liq(c,j) + l2e_h2osoi_ice(c,j))
       end do
    end do

    ! Initialize ET sink to ZERO
    do g = bounds_clump%begg, bounds_clump%endg
       do j = 1,nlevmapped
          g_idx = (g - bounds_clump%begg)*nlevmapped + j
          qflx_elm_loc(g_idx) = 0.0_r8
       end do
    end do
#ifdef DEBUG_ELMPFEH
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] check elm_pf_idata update0.1'
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] qflx_elm_loc = ', qflx_elm_loc
   !  !stop
#endif
    ! Account for following fluxes in the top soil layer:
    ! - infiltration
    ! - dew
    ! - disapperance of snow layer
    ! - sublimation of snow
    j = 1
    col_wtgcell = 1._r8
    do c = bounds_proc_begc, bounds_proc_endc
       if (col_active(c) == 1) then
          ! Set gridcell indices
          g = col_gridcell(c)
          g_idx = (g - bounds_clump%begg)*nlevmapped + j
          c_idx = c - begc+1
          qflx_elm_loc(g_idx) = qflx_elm_loc(g_idx) + &
               (&
               mflx_infl_col_1d(c_idx)    + &
               mflx_dew_col_1d(c_idx)     + &
               mflx_snowlyr_col_1d(c_idx) + &
               mflx_sub_snow_col_1d(c_idx)  &
               )*col_wtgcell*area_elm_loc(g_idx)
       end if
       total_mass_flux_col(c) = 0.d0
    enddo
#ifdef DEBUG_ELMPFEH
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] check elm_pf_idata update0.2'
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] qflx_elm_loc = ', qflx_elm_loc
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] mflx_infl_col_1d = ', mflx_infl_col_1d
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] mflx_dew_col_1d = ', mflx_dew_col_1d
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] mflx_snowlyr_col_1d = ', mflx_snowlyr_col_1d
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] mflx_sub_snow_col_1d = ', mflx_sub_snow_col_1d
   !  !stop
#endif
    ! Account for evapotranspiration flux
    do c = bounds_proc_begc, bounds_proc_endc
       do j = 1,nlevmapped
          g = col_gridcell(c)
          g_idx = (g - bounds_clump%begg)*nlevmapped + j
          c_idx = (c - bounds_proc_begc)*nlevgrnd+j
          if (col_active(c) == 1) then
             qflx_elm_loc(g_idx) = qflx_elm_loc(g_idx) + &
                  mflx_et_col_1d(c_idx)*area_elm_loc(g_idx)
             total_mass_flux_col(c) = total_mass_flux_col(c) + qflx_elm_loc(g_idx)/area_elm_loc(g_idx)
          end if
       end do
    end do
#ifdef DEBUG_ELMPFEH
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] check elm_pf_idata update1'
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] qflx_elm_loc = ', qflx_elm_loc
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] mflx_et_col_1d = ', mflx_et_col_1d
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] area_elm_loc = ', area_elm_loc
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] thetares2_elm_loc = ', thetares2_elm_loc
    !stop
#endif
    call VecRestoreArrayF90(elm_pf_idata%qflx_elm, qflx_elm_loc, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%thetares2_elm, thetares2_elm_loc, ierr); CHKERRQ(ierr)

    call pflotranModelUpdateFlowConds( this%pflotran_m )
#ifdef DEBUG_ELMPFEH
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] before pflotranModelStepperRunTillPauseTime'
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] nstep+1.0d0*dtime = ', (nstep+1.0d0)*dtime
    !stop
#endif
    call pflotranModelStepperRunTillPauseTime( this%pflotran_m, (nstep+1.0d0)*dtime )
    call pflotranModelGetUpdatedData( this%pflotran_m )

#ifdef DEBUG_ELMPFEH
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] pass pflotranModelStepperRunTillPauseTime'
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] nstep+1.0d0*dtime = ', (nstep+1.0d0)*dtime
    !stop
#endif
    call VecGetArrayF90(elm_pf_idata%sat_elms   , sat_elm_loc   , ierr); CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%mass_elms  , mass_elm_loc  , ierr); CHKERRQ(ierr)
    call VecGetArrayF90(elm_pf_idata%watsat_elm, watsat_elm_loc, ierr); CHKERRQ(ierr)
#ifdef DEBUG_ELMPFEH
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] check elm_pf_idata update2'
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] sat_elm_loc = ', sat_elm_loc
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] mass_elm_loc = ', mass_elm_loc
   !  write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] watsat_elm_loc = ', watsat_elm_loc
    !stop
#endif
    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)
       g = col_gridcell(c)

       ! initialization
       jwt = -1

       ! Loops in decreasing j so WTD can be computed in the same loop
       e2l_h2osoi_liq(c,:) = 0._r8
       e2l_h2osoi_ice(c,:) = 0._r8
       do j = nlevmapped, 1, -1
          g_idx = (g - bounds_clump%begg)*nlevmapped + j

          e2l_h2osoi_liq(c,j) =  (1.d0 - frac_ice(c,j))*mass_elm_loc(g_idx)/area_elm_loc(g_idx)
          e2l_h2osoi_ice(c,j) =  frac_ice(c,j)         *mass_elm_loc(g_idx)/area_elm_loc(g_idx)

          mass_end        = mass_end        + mass_elm_loc(g_idx)/area_elm_loc(g_idx)
          mass_end_col(c) = mass_end_col(c) + mass_elm_loc(g_idx)/area_elm_loc(g_idx)
#ifdef DEBUG_ELMPFEH
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] checkpoint e2l_h2osoi_liq, in do-loop'
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- c=', c
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- g=', g
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- j=', j
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- g_idx=', g_idx
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- e2l_h2osoi_liq(c,j)=', e2l_h2osoi_liq(c,j)
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro]     |- frac_ice(c,j)(c,j)=', frac_ice(c,j)
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro]     |- mass_elm_loc(g_idx)=', mass_elm_loc(g_idx)
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro]     |- area_elm_loc(g_idx)=', area_elm_loc(g_idx)
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro]     |- ratio=', mass_elm_loc(g_idx)/area_elm_loc(g_idx)
      ! write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- e2l_h2osoi_ice(c,j)=', e2l_h2osoi_ice(c,j)
#endif
       end do


       ! Find maximum water balance error over the column
      !  abs_mass_error_col = max(abs_mass_error_col,                     &
      !      abs(mass_beg_col(c) - mass_end_col(c) + &
      !      total_mass_flux_col(c)*dt))
       mass_bal_error_col(c) = mass_beg_col(c) - mass_end_col(c) + total_mass_flux_col(c)*dt
       total_mass_bal_error = total_mass_bal_error + mass_bal_error_col(c)
       e2l_qrecharge     (c) = 0._r8

       e2l_wtd(c) = l2e_zi(c,nlevmapped)
#ifdef DEBUG_ELMPFEH
      write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] water balance check col'
      write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- c=', c
      write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- mass_bal_error_col(c) =', mass_bal_error_col(c)
      write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- mass_beg_col(c) =', mass_beg_col(c)
      write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- mass_end_col(c) =', mass_end_col(c)
      write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- total_mass_flux_col(c)*dt =', total_mass_flux_col(c)*dt
      write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] |- total_mass_bal_error =', total_mass_bal_error
#endif
    end do

    ! Save soil liquid pressure from VSFM for all (active+nonactive) cells.
    ! soilp_col is used for restarting VSFM.
    do c = begc, endc
       do j = 1, nlevgrnd
          c_idx = (c - begc)*nlevgrnd + j
          e2l_soilp(c,j) = 0._r8!soilp_col_1d(c_idx)
       end do
    end do

    ! e2l_drain_perched (:) = 0._r8
    ! e2l_drain         (:) = 0._r8
    ! e2l_qrgwl         (:) = 0._r8
    ! e2l_rsub_sat      (:) = 0._r8
#ifdef DEBUG_ELMPFEH
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] check elm_pf_idata update3'
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] area_elm_loc = ', area_elm_loc
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] sat_elm_loc = ', sat_elm_loc
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] mass_elm_loc = ', mass_elm_loc
    write(*,*) '[YX DEBUG][ExternalModelPFLOTRANMod::EM_PFLOTRAN_Solve_Soil_Hydro] watsat_elm_loc = ', watsat_elm_loc
    !stop
#endif

    call VecRestoreArrayF90(elm_pf_idata%area_top_face_elms, area_elm_loc, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%sat_elms   , sat_elm_loc   , ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%mass_elms  , mass_elm_loc  , ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(elm_pf_idata%watsat_elm, watsat_elm_loc, ierr); CHKERRQ(ierr)

#ifdef PRINT_INTERNALFLOW
    call pflotranModelGetInternalflow( this%pflotran_m )
#endif

    deallocate(frac_ice                    )
    deallocate(total_mass_flux_col         )
    deallocate(total_mass_flux_et_col      )
    deallocate(total_mass_flux_infl_col    )
    deallocate(total_mass_flux_dew_col     )
    deallocate(total_mass_flux_drain_col   )
    deallocate(total_mass_flux_snowlyr_col )
    deallocate(total_mass_flux_sub_col     )
    deallocate(total_mass_flux_lateral_col )
    deallocate(total_mass_flux_seepage_col )
    deallocate(qflx_seepage                )
    deallocate(mass_prev_col          )
    deallocate(dmass_col              )
    deallocate(mass_beg_col                )
    deallocate(mass_end_col                )
    deallocate(mass_bal_error_col          )

    deallocate(mflx_et_col_1d              )
    deallocate(mflx_drain_col_1d           )
    deallocate(mflx_infl_col_1d            )
    deallocate(mflx_dew_col_1d             )
    deallocate(mflx_sub_snow_col_1d        )
    deallocate(mflx_snowlyr_col_1d         )
    deallocate(t_soil_col_1d               )

    deallocate(mass_col_1d            )
    deallocate(fliq_col_1d            )
    deallocate(smpl_col_1d            )
    deallocate(soilp_col_1d           )
    deallocate(sat_col_1d             )

  end subroutine EM_PFLOTRAN_Solve_Soil_Hydro

  subroutine EM_PFLOTRAN_Finalize(this)
  !
  ! Finalizes the ELM-PFLOTRAN coupling
  !
  ! Author: Yi Xiao
  ! Date: 8/23/2024
  !
    use pflotran_model_module         , only : pflotranModelDestroy, pflotranModelStepperRunFinalize
    use elm_pflotran_interface_data   , only : ELMPFLOTRANIDataDestroy

    implicit none
    ! !ARGUMENTS:
    class(em_pflotran_type)              :: this
    PetscErrorCode                       :: ierr

    !! Finalize PFLOTRAN Stepper
    !call pflotranModelStepperRunFinalize(this%pflotran_m)

    call ELMPFLOTRANIDataDestroy()

    call pflotranModelDestroy(this%pflotran_m)

    !call MPI_Finalize(ierr);CHKERRQ(ierr)
  end subroutine EM_PFLOTRAN_Finalize

#endif
end module ExternalModelPFLOTRANMod
