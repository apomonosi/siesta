module qmmm_lj
  !! Calculates QM-MM Lennard-Jones interactions.
  !! Also QM-QM via Grimme's model.
  implicit none

  public :: ljef

  ! GRIMME SHOULD NOT BE HERE AT ALL.

contains
  subroutine ljef( na_qm, natot, r, qm_atoms, mm_atoms, f, stress, Es, &
                   rcut_qmmm, cell, lattice_type )
    !! Calculates the LJ interaction energy and forces between the QM and
    !! the MM regions, and QM-QM via a simple Grimme's model.
    use precision        , only : dp
    use mm_topology      , only : qm_atom_t, mm_atom_t
    use qmmm_mm_neighbour, only : qmmm_veclist, qmmm_veclistxat, qmmm_nr
    use units            , only : Ang

    implicit none
    integer         , intent(in)    :: na_qm
      !! Number of QM atoms.
    integer         , intent(in)    :: natot
      !! Total number of atoms, QM+MM.
    real(dp)        , intent(in)    :: r(3,natot)
      !! Atomic positions.
    type(qm_atom_t) , intent(in)    :: qm_atoms(na_qm)
      !! Data for QM atoms.
    type(mm_atom_t) , intent(in)    :: mm_atoms(natot-na_qm)
      !! Data for MM atoms.
    real(dp)        , intent(in)    :: cell(3,3)
      !! Cell vectors.
    real(dp)        , intent(in)    :: rcut_qmmm
      !! Cut-off radius for QM-MM interactions.
    character(len=1), intent(in)    :: lattice_type
      !! Type of lattice.
    real(dp)        , intent(inout) :: f(3,natot)
      !! Atomic forces.
    real(dp)        , intent(inout) :: Es
      !! Output for Lennard-Jones energy.
    real(dp)        , intent(inout) :: stress(3,3)
      !! Cell stress.

    integer  :: n_pointer, iat, jat, kat, iCrd, iVec
    real(dp) :: lj_a, lj_b, lj_e, lj_s, Elj, flj(3), drij(3), dd, fej, &
                rcut_qmmm2, stress_fact

    ! External from siesta.
    real(dp), external  :: volcel

    Elj = 0.0_dp
    flj = 0.0_dp

    ! Here some work can be done to decouple the rcut_qmmm below from the
    ! cutoff used to define the real and reciprocal ewald sums in pcpot
    ! and mmforce.
    rcut_qmmm2  = ( rcut_qmmm * Ang ) ** 2
    stress_fact = 1.0_dp / volcel(cell)

    ! loop over QM atoms without considering the link atoms.
    ! n_pointer: points the first neighbour atom of i in the neighbour list.
    n_pointer = 1
    do iat = 1, na_qm
      do kat = n_pointer, qmmm_veclistxat(iat)
        jat = qmmm_veclist(kat)
        if ( jat <= na_qm ) cycle

        if ( lattice_type == 'D' ) then
          drij(1) = r(1,iat) - r(1,jat) + qmmm_nr(1,kat) * cell(1,1)
          drij(2) = r(2,iat) - r(2,jat) + qmmm_nr(2,kat) * cell(2,2)
          drij(3) = r(3,iat) - r(3,jat) + qmmm_nr(3,kat) * cell(3,3)
        else
          drij(1) = r(1,iat) - r(1,jat)
          drij(2) = r(2,iat) - r(2,jat)
          drij(3) = r(3,iat) - r(3,jat)

          do iCrd = 1, 3
          do iVec = 1, 3
            drij(iCrd) = drij(iCrd) + qmmm_nr(iVec,kat) * cell(iCrd,iVec)
          enddo
          enddo
        endif

        dd = sqrt( drij(1) * drij(1) + drij(2) * drij(2) + drij(3) * drij(3) )
        if ( (dd * dd) < rcut_qmmm2 ) then
          ! Energy and forces from LJ term
          ! Definition of the mixing rule for QM and MM epsilon and sigma.

          lj_e = sqrt( mm_atoms(jat-na_qm)%lj_Em * qm_atoms(iat)%lj_Em )
          if ( lj_e < 1.0e-14_dp ) cycle

          lj_s = 0.5_dp * ( mm_atoms(jat-na_qm)%lj_Rm + qm_atoms(iat)%lj_Rm )
          lj_b = 4.0_dp * lj_e * lj_s ** 6
          lj_a = 4.0_dp * lj_e * lj_s ** 12

          !! For QM-MM we use regular LJ.
          Elj = Elj + lj_a / dd**12 - lj_b / dd**6
          fej = -12.0_dp * lj_a / dd**14 + 6.0_dp * lj_b / dd**8

          flj(1) = -2.0_dp * fej * drij(1)
          flj(2) = -2.0_dp * fej * drij(2)
          flj(3) = -2.0_dp * fej * drij(3)

          do iCrd = 1, 3
            f(iCrd,iat) = f(iCrd,iat) + flj(iCrd)
            f(iCrd,jat) = f(iCrd,jat) - flj(iCrd)

            do iVec = 1, 3
              stress(iVec,iCrd) = stress(iVec,iCrd) &
                             + stress_fact * drij(iVec) * flj(iCrd)
            enddo
          enddo
        endif
      enddo !! Neighbour QM+MM atoms

      n_pointer = qmmm_veclistxat(iat) + 1
    enddo !! QM atoms
    Es = 2.0_dp * Elj
  end subroutine ljef
end module qmmm_lj
