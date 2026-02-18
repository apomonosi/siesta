This directory contains sources for comprehensive handling of the units.

These sources are client-side modules containing the materials-physics
"units table" which was formerly located in the fdf package.

Now all units handling is provided by the hosting program (Siesta)
and will enable easy updates of unit tables etc.

To use units through fdf we pass a function `inquire_unit` as
a procedure pointer to allow fdf to call client code units handling.

The function has this interface
  
```fortran
  subroutine inquire_unit(unit_str, stat, phys_dim, unit_name, unit_value)

    character(len=*), intent(in)   :: unit_str   ! unit specification
    character(len=*), intent(out)  :: phys_dim   ! physical dimension (e.g. 'mass')
    character(len=*), intent(out)  :: unit_name  ! unit name (e.g. 'g')
    real(dp), intent(out)          :: unit_value ! actual value (e.g. 1.e-3)
    integer, intent(out)           :: stat       ! status code

```

This also means that the client-side implementation of the function
should handle all the errors and warning without fdf knowing about fall-backs
etc.

To produce a new table please edit the units.py and add the new data format
there.
Generically this small python code should be used to generate the unit-tables
as it can create consistent data-types and output for future works.
