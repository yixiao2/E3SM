module elm_finalizeMod

  !-----------------------------------------------------------------------
  ! Performs land model cleanup
  !
  !
  implicit none
  save
  public   ! By default everything is public
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------
  !
  public :: final
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine final( )
    !
    ! !DESCRIPTION:
    ! Finalize land surface model
    !
#ifdef USE_PETSC_LIB
#include <petsc/finclude/petsc.h>
#endif
    ! !USES:
    use elm_varctl             , only : use_vsfm, use_pflotran_hmode_via_emi
#ifdef USE_PETSC_LIB
    use petscsys
    use ExternalModelInterfaceMod, only : EMI_Finalize_For_PFLOTRAN
#endif
    !
    ! !ARGUMENTS
    implicit none
    !

#ifdef USE_PETSC_LIB
    PetscErrorCode        :: ierr

    if (use_vsfm) then
      call PetscFinalize(ierr)
    elseif (use_pflotran_hmode_via_emi) then
#ifdef DEBUG_ELMPFEH
      write(*,*) '[YX DEBUG][elm_finalizeMod::final]'
      stop
#endif
      call EMI_Finalize_For_PFLOTRAN()
      call PetscFinalize(ierr)
    endif
#endif

  end subroutine final

end module elm_finalizeMod
