!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module io_ext_m

   implicit none

   public

   public :: io_fc_indices
   public :: io_ext_fc
   public :: io_ext_md

contains

   subroutine io_fc_indices(istep, ia1, fc_ia, fc_istep)
      !< The current step of the simulation, starting at 0 (initial configuration)
      integer, intent(in) :: istep
      !< Which atom that is started displacing.
      integer, intent(in) :: ia1
      !< The current atom being displaced
      integer, intent(out) :: fc_ia
      !< The corrected step, more or less (% 6)
      integer, intent(out) :: fc_istep

      if ( istep == 0 ) then
         fc_istep = 0
         fc_ia = ia1
      else
         fc_istep = mod(istep-1,6) + 1
         fc_ia = (istep - mod(istep-1,6))/6 + ia1
      end if

   end subroutine

   !< Create an extension character that unifies the way we built FC file names
   !<
   !< This will create filenames:
   !<
   !< - 00000.EXT
   !< - 00010-1.EXT
   !< - ...
   !< - 00010-6.EXT
   !< - 00011-1.EXT
   !< - ...
   !< - 00011-6.EXT
   !<
   !< Where `10` is the initial atom being displaced.
   subroutine io_ext_fc(ext, istep, ia1, full_ext)
      !< Default extension name, will be suffixed for the correct extension.
      character(len=*), intent(in) :: ext
      !< The current step of the simulation, starting at 0 (initial configuration)
      integer, intent(in) :: istep
      !< Which atom that is started displacing.
      integer, intent(in) :: ia1

      !< Resulting full extension depending on `istep` and `ia1`.
      character(len=*), intent(out) :: full_ext

      integer :: fc_ia
      integer :: fc_istep

      call io_fc_indices(istep, ia1, fc_ia, fc_istep)

      if ( fc_istep == 0 ) then
         ! This is our initial configuration before any displacements.
         ! The output extension will then be .00000.<EXT>
         write(full_ext, '(i5.5,2a)') 0, '.', trim(ext)
      else
         ! Subsequent displacements
         write(full_ext, '(i5.5,a,i1,2a)') fc_ia, '-', fc_istep, '.', trim(ext)
      end if

   end subroutine

   !< Create an extension character that unifies the way we built MD file names
   !<
   !< This will create filenames:
   !<
   !< - 0.EXT
   !< - 1.EXT
   !< - ...
   !< - 11.EXT
   !<
   !< Where `11` is the final MD step.
   subroutine io_ext_md(ext, istep, full_ext)
      !< Default extension name, will be suffixed for the correct extension.
      character(len=*), intent(in) :: ext
      !< The current step of the simulation, starting at 0 (initial configuration)
      integer, intent(in) :: istep

      !< Resulting full extension depending on `istep`.
      character(len=*), intent(out) :: full_ext

      write(full_ext, '(i0,2a)') istep, '.', trim(ext)

   end subroutine

end module
