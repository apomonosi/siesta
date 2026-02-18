! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_setup_hamiltonian
      private
      public :: setup_hamiltonian
      CONTAINS

      subroutine setup_hamiltonian( iscf )

      use option_charges_m, only: OPT_CHARGE_SCF
      USE siesta_options
      use sparse_matrices, only: H_2D, S_1D, xij_2D
      use sparse_matrices, only: H_kin_1D, H_vkb_1D
      use sparse_matrices, only: H_dftu_2D, H_dftu_so_2D
      use sparse_matrices, only: H_so_on_2D, H_so_off_2D
      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use sparse_matrices, only: H, S, Hold
      use sparse_matrices, only: Dscf, Escf, xijo
      use class_dSpData1D,  only: val
      use class_dSpData2D,  only: val
      use class_zSpData2D,  only: val

      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     .                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     .                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use dftu_specs,   only: switch_dftu     ! This variable determines whether
                                              !   the subroutine to compute the
                                              !   Hubbard terms should be called
                                              !   or not
      use m_dftu,       only: hubbard_term    ! Subroutine that compute the
                                              !   Hubbard terms
      use m_dhscf,      only: dhscf
      use m_stress
      use m_energies
      use parallel, only: Node, Nodes
      use m_steps, only: istp
      use m_ntm

      use m_spin,          only: spin 

      use m_dipol
      use alloc, only: re_alloc, de_alloc
      use io_hsx_m, only: write_hsx
      use sys, only: die, bye
      use partial_charges_m, only: want_partial_charges
      use files, only : filesOut_t    ! derived type for output file names
      use m_rhog, only: rhog_in, rhog
      use m_dftu_so, only: dftu_so_hamil
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
      use mpi_siesta
#endif
      use kpoint_scf_m, only: kpoint_scf

#ifdef SIESTA__WANNIER90
      use w90_wrapper_types,   only: compute_chempotwann_after_scf
      use w90_wrapper_types,   only: first_chempotwann
      use m_chempotwann,     only: compute_H_chempotwann
      use m_chempotwann,     only: add_Hamiltonian_chempotwann
      use sparse_matrices,   only: DM_2D
      use sparse_matrices,   only: EDM_2D
      use sparse_matrices,   only: H_chempotwann_2D
      use sparse_matrices,   only: DM_history
      use class_dSpData2D,   only: val
      use m_new_dm,          only: new_dm
#endif

      implicit none
      integer, intent(in) :: iscf
      real(dp)            :: stressl(3,3)
      real(dp), pointer   :: fal(:,:)   ! Local-node part of atomic F
#ifdef MPI
      real(dp)            :: buffer1
      integer             :: MPIError
#endif
      integer             :: io, is, ispin
      integer             :: ifa     ! Calc. forces?      0=>no, 1=>yes
      integer             :: istr    ! Calc. stress?      0=>no, 1=>yes
      integer             :: ihmat   ! Calc. hamiltonian? 0=>no, 1=>yes
      real(dp)            :: g2max
      type(filesOut_t)    :: filesOut  ! blank output file names
      logical             :: use_rhog_in

      real(dp), pointer   :: H_vkb(:), H_kin(:), H_dftu(:,:)
      real(dp), pointer   :: H_so_on(:,:)
      complex(dp), pointer:: H_dftu_so(:,:)
      complex(dp), pointer:: H_so_off(:,:)

      complex(dp):: Dc
      integer :: ind, i, j, jo

#ifdef SIESTA__WANNIER90
      integer             :: index_manifold
      real(dp), pointer   :: H_chempotwann(:,:)
#endif

!------------------------------------------------------------------------- BEGIN

      call timer('setup_H',1)

      ! Nullify pointers
      nullify(fal)

!$OMP parallel default(shared), private(ispin,io)

!     Save present H matrix
      do ispin = 1, spin%H
!$OMP do
         do io = 1,maxnh
            Hold(io,ispin) = H(io,ispin)
         end do
!$OMP end do nowait
      end do

!$OMP single
      H_kin => val(H_kin_1D)
      H_vkb => val(H_vkb_1D)
      
      if ( spin%SO_onsite ) then
        ! Sadly some compilers (g95), does
        ! not allow bounds for pointer assignments :(
        H_so_on => val(H_so_on_2D)
        
      else if ( spin%SO_offsite ) then
        H_so_off => val(H_so_off_2D)
        
      end if
!$OMP end single ! keep wait

      ! Initialize diagonal Hamiltonian
      do ispin = 1, spin%spinor
!$OMP do
        do io = 1,maxnh
          H(io,ispin) = H_kin(io) + H_vkb(io)
        end do
!$OMP end do nowait
      end do

      if ( spin%SO_onsite ) then
        do ispin = 3 , spin%H
!$OMP do
          do io = 1,maxnh
            H(io,ispin) = H_so_on(io,ispin-2)
          end do
!$OMP end do nowait
        end do
         
      else

        do ispin = 3 , spin%H
!$OMP do
          do io = 1,maxnh
            H(io,ispin) = 0._dp
          end do
!$OMP end do nowait
        end do

      end if

! ..................

! Non-SCF part of total energy .......................................
! Note that these will be "impure" for a mixed Dscf

! If mixing the charge, Dscf is the previous step's DM_out. Since
! the "scf" components of the energy are computed with the (mixed)
! charge, this introduces an inconsistency. In this case the energies
! coming out of this routine need to be corrected.
! 
!$OMP single
      Ekin = 0.0_dp
      Enl = 0.0_dp
      Eso = 0.0_dp
!$OMP end single ! keep wait
      
!$OMP do reduction(+:Ekin,Enl)
      do io = 1,maxnh
        do ispin = 1, spin%spinor
          Ekin = Ekin + H_kin(io) * Dscf(io,ispin)
          Enl  = Enl  + H_vkb(io) * Dscf(io,ispin)
        end do
      end do
!$OMP end do nowait

!
!  Dc IS NOT the dense matrix, it is just a complex number 
! (per each io index) used as an artifact to multiply the 
! elements of the H_SO times the corresponding elements of 
! DM in a such way that the result gives Re{Tr[H_SO*DM]}.
!

      if ( spin%SO_offsite ) then
!!       For debugging
!        do io = 1, no_u
!          do j = 1, numh(io)
!            ind = listhptr(io) + j
!            jo = listh(ind)
!            if( io .eq. jo ) then
!              H_so_off(ind,*)= cmplx(0.0_dp,0.0_dp,kind=dp)
!            endif
!          enddo 
!        enddo 
!!       End debugging

        do io = 1, maxnh
!-------- Eso(u,u)
          Dc = cmplx(Dscf(io,1),Dscf(io,5), dp)
          Eso = Eso + real( H_so_off(io,1)*Dc, dp)
!-------- Eso(d,d)
          Dc = cmplx(Dscf(io,2),Dscf(io,6), dp)
          Eso = Eso + real( H_so_off(io,2)*Dc, dp)
!-------- Eso(u,d)
          Dc = cmplx(Dscf(io,3),Dscf(io,4), dp)
          Eso = Eso + real( H_so_off(io,4)*Dc, dp)
!-------- Eso(d,u)
          Dc = cmplx(Dscf(io,7),-Dscf(io,8), dp)
          Eso = Eso + real( H_so_off(io,3)*Dc, dp)
        end do

      else if ( spin%SO_onsite ) then
        
!$OMP do reduction(+:Eso)
         do io = 1, maxnh
           Eso = Eso + H_so_on(io,1)*Dscf(io,7) +
     &         H_so_on(io,2)*Dscf(io,8)+ H_so_on(io,5)*Dscf(io,3) +
     &         H_so_on(io,6)*Dscf(io,4)- H_so_on(io,3)*Dscf(io,5) -
     &         H_so_on(io,4)*Dscf(io,6)
         end do
!$OMP end do nowait

      end if

!$OMP end parallel

#ifdef MPI
      ! Global reduction of Ekin, Enl
      call globalize_sum(Ekin,buffer1)
      Ekin = buffer1
      call globalize_sum(Enl,buffer1)
      Enl = buffer1
      if ( spin%SO ) then
         ! Global reduction of Eso 
         call globalize_sum(Eso,buffer1)
         Eso = buffer1
      end if
#endif

!     Non-SCF part of total energy
      call update_E0()

! Hubbard term for DFT+U: energy, forces, stress and matrix elements ....
      if( switch_dftu ) then
        call re_alloc(fal, 1, 3, 1, na_u, 'fal', 'setup_hamiltonian')      
        if ( spin%NCol ) then
          call die('DFT+U cannot be used with non-collinear spin.')
        else if ( spin%SO ) then
          H_dftu_so => val(H_dftu_so_2D)
          call dftu_so_hamil( H_dftu_so, fal, stressl )
!        When SO coupling is considered, the density matrix and the Hamiltonian
!        must be globally Hermitian (see last sentence of Section 7 of the 
!        technical SIESTA paper in JPCM 14, 2745 (2002).
!        However, in diag3k, the sign of the imaginary part of the (up,down) 
!        matrix elements is changed:
!        Dnew(ind,1) = Dnew(ind,1) + real(D11,dp)
!        Dnew(ind,2) = Dnew(ind,2) + real(D22,dp)
!        Dnew(ind,3) = Dnew(ind,3) + real(D12,dp)
!        Dnew(ind,4) = Dnew(ind,4) - aimag(D12)
!        Dnew(ind,5) = Dnew(ind,5) + aimag(D11)
!        Dnew(ind,6) = Dnew(ind,6) + aimag(D22)
!        Dnew(ind,7) = Dnew(ind,7) + real(D21,dp)
!        Dnew(ind,8) = Dnew(ind,8) + aimag(D21)
!        In the subroutines to compute the corresponding DFT+U matrix elements:
!        - We change locally the sign of this imaginary part, so the 
!          density matrix recovers all its properties
!        - The Potential matrix elements are computed with the "good" DM.
!        - Here, at the time of adding up the new Hamiltonian matrix elements
!          to the potential, we change the sign of the imaginary part of the
!          (up,down) component

!        
!------- H(u,u)
         H(:,1) = H(:,1) + real(H_dftu_so(:,1), dp)
         H(:,5) = H(:,5) + aimag(H_dftu_so(:,1))
!------- H(d,d)
         H(:,2) = H(:,2) + real(H_dftu_so(:,2), dp)
         H(:,6) = H(:,6) + aimag(H_dftu_so(:,2))
!------- H(u,d)
         H(:,3) = H(:,3) + real(H_dftu_so(:,3), dp)
         H(:,4) = H(:,4) - aimag(H_dftu_so(:,3))
!------- H(d,u)
         H(:,7) = H(:,7) + real(H_dftu_so(:,4), dp)
         H(:,8) = H(:,8) + aimag(H_dftu_so(:,4))

#ifdef MPI
          ! Global reduction of E_dftu_so 
          call globalize_sum(E_dftu_so,buffer1)
          E_dftu_so = buffer1
          ! Global reduction of E_correc_dc
          call globalize_sum(E_correc_dc,buffer1)
          E_correc_dc = buffer1
#endif
        else  ! collinear spin
          H_dftu => val(H_dftu_2D)
          call hubbard_term(scell, na_u, na_s, isa, xa, indxua,
     .                      maxnh, maxnh, lasto, iphorb, no_u, no_l, 
     .                      numh, listhptr, listh, numh, listhptr,listh,      
     .                      spin%spinor, Dscf, Edftu, DEdftu, H_dftu, 
     .                      fal, stressl, H, iscf,
     .                      matrix_elements_only=.true.)

#ifdef MPI
          ! Global reduction of energy terms
          call globalize_sum(Edftu,buffer1)
          Edftu = buffer1
          ! DEdftu should not be globalized
          ! as it is based on globalized occupations
#endif
          Edftu = Edftu + DEdftu

        endif
        call de_alloc( fal, 'fal', 'setup_hamiltonian' ) 
      endif
! ..................

! Add SCF contribution to energy and matrix elements ..................
      g2max = g2cut

      call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'setup_hamiltonian' )

      ifa  = 0
      istr = 0
      ihmat = 1

      want_partial_charges = 0
      if ( option_charges%hirshfeld%run(OPT_CHARGE_SCF) ) then
         want_partial_charges = OPT_CHARGE_SCF
      end if
      if ( option_charges%voronoi%run(OPT_CHARGE_SCF) ) then
         want_partial_charges = OPT_CHARGE_SCF
      end if
      use_rhog_in =  (mix_charge .and. iscf > 1)
  
      call dhscf( spin%Grid, no_s, iaorb, iphorb, no_l,
     .            no_u, na_u, na_s, isa, xa, indxua, 
     .            ntm, ifa, istr, ihmat, filesOut,
     .            maxnh, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, stress, fal, stressl,
     .            use_rhog_in)

      want_partial_charges = 0

      if ( spin%SO_offsite ) then

! H(:, [5, 6]) are not updated in dhscf, see vmat for details.
        
!------- H(u,u)
         H(:,1) = H(:,1) + real(H_so_off(:,1), dp)
         H(:,5) = H(:,5) + aimag(H_so_off(:,1))
!------- H(d,d)
         H(:,2) = H(:,2) + real(H_so_off(:,2), dp)
         H(:,6) = H(:,6) + aimag(H_so_off(:,2))
!------- H(u,d)
         H(:,3) = H(:,3) + real(H_so_off(:,3), dp)
         H(:,4) = H(:,4) + aimag(H_so_off(:,3))
!------- H(d,u)
         H(:,7) = H(:,7) + real(H_so_off(:,4), dp)
         H(:,8) = H(:,8) - aimag(H_so_off(:,4))

!!        For debugging
!         do io = 1, no_u
!           do j = 1, numh(io)
!             ind = listhptr(io) + j
!             jo = listh(ind)
!             write(6,'(a,6i7,8f12.5)') 
!     .         ' Node, Nodes, io, jo, j, ind, Dscf = ' , 
!     .           Node, Nodes, io, jo, j, ind, Dscf(ind,:)
!           enddo 
!         enddo 
!         do io = 1, no_u
!           do j = 1, numh(io)
!             ind = listhptr(io) + j
!             jo = listh(ind)
!             write(6,'(a,6i7,8f12.5)') 
!     .         ' Node, Nodes, io, jo, j, ind, H-DFTU = ' , 
!     .           Node, Nodes, io, jo, j, ind, H_dftu_so(ind,:)
!           enddo 
!         enddo 
!         do io = 1, no_u
!           do j = 1, numh(io)
!             ind = listhptr(io) + j
!             jo = listh(ind)
!             write(6,'(a,6i7,8f12.5)') 
!     .         ' Node, Nodes, io, jo, j, ind, H-SO = ' , 
!     .           Node, Nodes, io, jo, j, ind, H_so_off(ind,:)
!           enddo 
!         enddo 
!         do io = 1, no_u
!           do j = 1, numh(io)
!             ind = listhptr(io) + j
!             jo = listh(ind)
!             write(6,'(a,6i7,8f12.5)') 
!     .         ' Node, Nodes, io, jo, j, ind, H = ' , 
!     .           Node, Nodes, io, jo, j, ind, H(ind,:)
!           enddo 
!         enddo 
!!        End debugging

      endif

! Energy term coming from the chemical potential added to the Wannier functions
#ifdef SIESTA__WANNIER90
      if( w90_wrapper_processing ) then
        if( compute_chempotwann_after_scf ) then
          H_chempotwann => val(H_chempotwann_2D)
          if( first_chempotwann ) then
!           In this subroutine an extra contribution to the 
!           DFT hamiltonian is computed.
!           This extra contribution rises or lowers the diagonal energy element
!           of the Hamiltonian expressed in a Wannier function basis.
!           The magnitude of this shifting is introduced in the
!           Wannier.ChemicalPotential block
            call compute_H_chempotwann( H_chempotwann )
            first_chempotwann = .false.
          endif  ! Endif first_chempotwann
          call add_Hamiltonian_chempotwann( H_chempotwann )
!! For debugging
!#ifdef MPI
!          call MPI_barrier(MPI_Comm_world,MPIError)
!#endif
!          call die()  
!! End debugging
        endif    ! Endif compute_chempotwann_after_scf
      endif      ! Endif w90_wrapper_processing
#endif

! ..................


      ! This statement will apply to iscf = 1, for example, when
      ! we do not use rhog_in. Rhog here is always the charge used to
      ! build H, that is, rhog_in.
      if (mix_charge) rhog_in = rhog

      call de_alloc( fal, 'fal', 'setup_hamiltonian' ) 

      call timer('setup_H',2)
#ifdef SIESTA__PEXSI
      if (node==0) call memory_snapshot("after setup_H")
#endif

      if ( h_setup_only ) then

        ! Save Hamiltonian and overlap matrices ............................
        ! Only in HSX format now.  Use Util/HSX/hsx2hs to generate an HS file
        if (savehs .or. write_coop) then
!         We don't know the fermi-level of this one. So we will pass
!         a huge number to signal this.
          call write_hsx(H_2D, S_1D, huge(0._dp), qtot, temp,
     &        kpoint_scf, extension="0.HSX")
        end if

        call timer( 'all', 2 )  ! New call to close the tree
        call timer( 'all', 3 )
        call bye("H-Setup-Only requested")
        STOP
      endif

!------------------------------------------------------------------------- END
      END subroutine setup_hamiltonian
      END module m_setup_hamiltonian
