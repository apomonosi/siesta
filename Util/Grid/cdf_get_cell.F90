! 
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!--------------------
program cdf_get_cell
#ifndef CDF
  print *, "No netCDF support at compile time"
#else

use m_grid

type(grid_t)   :: gp

call get_cdf_grid("Rho.grid.nc",gp)

print *, "Unit cell: ", gp%cell
print *, "n(:):", gp%n

#endif /* CDF */
end program cdf_get_cell

