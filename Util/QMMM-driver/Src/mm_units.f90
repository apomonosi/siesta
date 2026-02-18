module mm_units
    !! Contains references to units commonly used in classical MD,
    !! and converts them to SIESTA units.
    use precision, only : dp
    use units    , only : eV, Ang, pi

    ! Converts kcal/mol to Hartree, eV or Siesta units.
    implicit none
    real(dp), parameter :: kcal_mol_eV = 23.060541945329334_dp
    real(dp), parameter :: kcal_mol    = kcal_mol_eV / eV

    ! 1.10^-10 m/Ang * 8.8541878E-12_dp C^2/Jm / 1.602177E-19_dp C/e
    ! = 0.00552634808 e/(eV*Ang)
    ! Final result is in e/(Ry*Bohr).
    real(dp), parameter :: eps0 = (8.8541878_dp / 1602.177_dp) / (Ang * eV)


    real(dp), parameter :: rads = 180.0_dp / pi
contains

end module mm_units