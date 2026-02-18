! This is a small library to interface SIESTA with other MM drivers using
! simple subroutine interfaces. It handles launching SIESTA, stablishing
! communications via sockets (default) or pipes.
!
! It is important to note that SIESTA does NOT handle LJ interactions
! between the QM and MM parts, that should be taken care of by the MM
! driver itself.
!
! Relevant data for SIESTA launch is read from a siestaqmmm.fdf file
! within the working directory.
subroutine get_siesta_forces_and_stress( n_qm, n_mm, r_qm, r_mm, pc, cell, &
                                         last, energy, f_qm, f_mm, stress  )
  !! This subroutine acts as the main interfaces with other MM software.
  !! WARNING: Forces and stressses are strict outputs, meaning that they
  !!          will be overwritten.
  use precision   , only : dp
  use siestaqmmm_m, only : send_coordinates, recv_forces, &
                           close_communicators
  implicit none
  integer , intent(in)  :: n_qm
    !! Total number of atoms in the QM region.
  integer , intent(in)  :: n_mm
    !! Total number of atoms in the MM region.
  real(dp), intent(in)  :: r_qm(3,1:n_qm)
    !! Positions of atoms in the QM region.
  real(dp), intent(in)  :: r_mm(3,1:n_mm)
    !! Positions of atoms in the MM region.
  real(dp), intent(in)  :: pc(1:n_mm)
    !! Classical partial charges of atoms in the MM region.
  real(dp), intent(in)  :: cell(3,3)
    !! Size of the unit cell.
  logical , intent(in)  :: last
    !! Indicates whether this is the last MD step.
  real(dp), intent(out) :: energy
    !! QM/MM energy.
  real(dp), intent(out) :: f_qm(3,1:n_qm)
    !! Forces over QM atoms.
  real(dp), intent(out) :: f_mm(3,1:n_mm)
    !! Forces over MM atoms.
  real(dp), intent(out) :: stress(3,3)
    !! Cell stress.

  ! Sends coordinates to SIESTA.
  ! This routine also launches SIESTA if necessary.
  call send_coordinates( n_qm, n_mm, r_qm, r_mm, pc, cell, last )

  ! Receives forces and stresses from SIESTA.
  call recv_forces( n_qm, n_mm, f_qm, f_mm, stress, energy )

  ! If this is the last step, close the pipes or files.
  if ( last ) call close_communicators( )

end subroutine get_siesta_forces_and_stress

subroutine set_qmmm_input_file( input_file )
  !! Sets the input filename for this library; it defaults to
  !! siestaqmmm.fdf otherwise.
  use siestaqmmm_m, only : set_input_file

  implicit none
  character(len=*), intent(in) :: input_file
  !! Optional input filename for the options used in this library.

  call set_input_file( input_file )
end subroutine set_qmmm_input_file