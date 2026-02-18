module qmmm_files_m
  !! This module contains just information related to input and output files
  !! used in the QMMM driver.

  implicit none
  type, public :: file_names_t
    !! Data structure to store filenames.
    character(len=224) :: slabel
      !! Backup of the system label to avoid multiple passes.
    character(len=224) :: qmlabel
      !! System label for the QM siesta program.
    character(len=255) :: ene
      !! Energy and temperature output file.
    character(len=255) :: cpdb
      !! PDB file when running restrained optimizations.
    character(len=255) :: cene
      !! Energy output file when running restrained optimizations.
    character(len=255) :: pdbi
      !! Initial structure PDB file.
    character(len=255) :: pdbl
      !! Last structure PDB file.
    character(len=255) :: mcrd
      !! MDCRD file with the trajectory.
    character(len=255) :: wrt
      !! File to write constraints.
    character(len=255) :: input
      !! Input file (debug or standard input)
  contains
    procedure :: set_filenames
  end type file_names_t

  type(file_names_t), public :: qmmm_files
    !! Contains all filenames.

  private
contains

  subroutine set_filenames( self, slabel )
    !! Sets filenames for all outputs in order to avoid re-setting
    !! them every time.
    implicit none
    character(len=*)   , intent(in) :: slabel
      !! System label (root name).
    class(file_names_t), intent(inout) :: self
      !! Data structure for filenames.

    self%slabel  = trim(slabel)
    self%qmlabel = trim(slabel) // '.siesta'
    self%ene     = trim(slabel) // '.ene'
    self%cpdb    = trim(slabel) // '.ctr.pdb'
    self%cene    = trim(slabel) // '.ctr.ene'
    self%pdbi    = trim(slabel) // '.init.pdb'
    self%mcrd    = trim(slabel) // '.mdcrd'
    self%pdbl    = trim(slabel) // '.last.pdb'
    self%wrt     = trim(slabel) // '.wrt'
  end subroutine set_filenames

end module qmmm_files_m