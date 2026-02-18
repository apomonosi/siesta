module graphite_m
  !! This module deals with auxiliary routines for MM graphene setups.
  use mm_topology, only : atom_connect_t, mm_atom_t

  implicit none
  private
  public :: find_graphite_layers

contains
  subroutine find_graphite_layers( nac, mm_atoms, mm_connectivity )

    !! For each MM atom, this subroutine assigns the graphite layer to which
    !! each MM atom belongs.
    implicit none
    integer         , intent(in)    :: nac
      !! Number of MM atoms.
    type(mm_atom_t) , intent(inout) :: mm_atoms(nac)
      !! MM atom information.
    type(atom_connect_t), intent(inout) :: mm_connectivity(nac)
    !! Data structure containing the connectivity data for each MM atom.

    integer :: iat, num_of_graphite_layers

    num_of_graphite_layers = 0
    do iat = 1, nac
      if ( (mm_atoms(iat)%graph_layer /= 0) .or. &
           (mm_atoms(iat)%attype /= 'C#') ) cycle

      num_of_graphite_layers    = num_of_graphite_layers +1
      mm_atoms(iat)%graph_layer = num_of_graphite_layers
      call check_graphite_connectivity( iat, nac, mm_atoms, mm_connectivity )
    enddo
  end subroutine find_graphite_layers

  recursive subroutine check_graphite_connectivity( iat, nac, mm_atoms, mm_connectivity )
    !! Assigns the graphene layer number of a given atom by recursively exploring
    !! its connectivity.
    implicit none
    integer         , intent(in)    :: iat
      !! Currrent MM atom index.
    integer         , intent(in)    :: nac
      !! Number of MM atoms.
    type(mm_atom_t) , intent(inout) :: mm_atoms(nac)
      !! MM atom information.
    type(atom_connect_t), intent(inout) :: mm_connectivity(nac)
      !! Data structure containing the connectivity data for each MM atom.

    integer :: jat

    do jat = 1, mm_connectivity(iat)%nbonds
      if ( (mm_atoms( mm_connectivity(iat)%bond_at(jat) )%graph_layer /= 0) &
            .or. (mm_atoms(iat)%attype /= 'C#') ) cycle
      mm_atoms( mm_connectivity(iat)%bond_at(jat) )%graph_layer = &
      mm_atoms(iat)%graph_layer
      call check_graphite_connectivity( mm_connectivity(iat)%bond_at(jat), &
                                        nac, mm_atoms, mm_connectivity )
    enddo
  end subroutine check_graphite_connectivity
end module graphite_m