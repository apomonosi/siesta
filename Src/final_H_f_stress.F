! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_final_H_f_stress
      private
      public :: final_H_f_stress
      CONTAINS

      subroutine final_H_f_stress( istep, iscf , SCFconverged )
      use units, only: eV, Ang

      use files,        only : slabel
#ifdef NCDF_4
      use siesta_options, only: write_cdf
#endif
      use siesta_options, only: fixspin
      use siesta_options, only: savehs, write_coop, fc_save_dhs
      use siesta_options, only: g2cut, temp, nmove
      use siesta_options, only: recompute_H_after_scf
      use siesta_options, only: want_dftd3_dispersion
      use sparse_matrices, only: numh, listh, listhptr
      use sparse_matrices, only: H, Dscf, Escf, maxnh
      use sparse_matrices, only: H_dftu_2D, H_so_off_2D, H_dftu_so_2D
      use fc_matrices_m, only: dHdR_2D, fc_H_2D, dSdR_1D, fc_S_1D
      use fc_matrices_m, only: fc_nsc
      use class_dSpData1D, only: val, delete
      use class_dSpData2D, only: val, delete
      use class_zSpData2D, only: val

      use siesta_geom
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm,
     &                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     &                    rmaxo, no_l, iza
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_nlefsm,     only: nlefsm, nlefsm_SO_off
      use m_overfsm,    only: overfsm
      use m_kinefsm,    only: kinefsm
      use m_naefs,      only: naefs
      use m_dnaefs,     only: dnaefs
      use m_grdsam,     only: grdsam
      use dftu_specs,   only: switch_dftu     ! This variable determines whether
                                              !   the subroutine to compute the
                                              !   Hubbard terms should be called
                                              !   or not
      use m_dftu,       only: hubbard_term    ! Subroutine that compute the
                                              !   Hubbard terms
      use m_dftu_so, only: dftu_so_hamil
      use m_stress
      use m_energies
      use m_ntm
      use m_spin,          only: spin
      use spinorbit,       only: spinorb
      use m_dipol
      use m_forces,          only: fa
      use alloc, only: re_alloc, de_alloc
      use io_hsx_m, only: write_hsx
      use sys, only: die
      use fdf
      use kpoint_t_m, only: kpoint_t
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif
      use parallel, only: IOnode
#ifdef CDF
#ifdef NCDF_4
      use dictionary
      use m_ncdf_siesta
      use dHSdR_m
#endif
#endif
      use siesta_options, only : idyn, ia1
      use siesta_options, only : write_hs_history
      use sparse_matrices, only: H_2D, S_1D
      use files,        only : label_length
      use kpoint_scf_m, only: kpoint_scf
      use ts_kpoint_scf_m, only: ts_kpoint_scf
      use m_ts_global_vars, only : TSrun
      use m_ts_options, only : TS_HS_save, write_tshs_history
      use io_ext_m, only : io_ext_md, io_ext_fc, io_fc_indices
      use ts_io_hs_m, only : ts_write_TSHS
#ifdef SIESTA__DFTD3
      use grimme_dispersion_m, only : dftd3_energy_forces
#endif

#ifdef SIESTA__WANNIER90
      use siesta_options,      only: w90_wrapper_processing
      use w90_wrapper_types, only: compute_chempotwann_after_scf
      use m_chempotwann,       only: add_Hamiltonian_chempotwann
      use sparse_matrices,     only: H_chempotwann_2D
      use class_dSpData2D,     only: val
#endif

#ifdef FINAL_CHECK_HS
      use m_compute_max_diff, only: compute_max_diff
#endif

      implicit none

      ! MD-step, SCF-step
      integer,   intent(in) :: istep, iscf
      logical,   intent(in) :: SCFconverged
      real(dp)              :: stressl(3,3)
      real(dp), pointer     :: fal(:,:)   ! Local-node part of atomic F
      integer               :: ifa     ! Calculate forces? 0=>no, 1=>yes
      integer               :: istr    ! Calculate stress? 0=>no, 1=>yes
      real(dp)              :: g2max
      ! To avoid overwriting the current Hamiltonian and S
      real(dp), pointer :: H_new(:,:), H_tmp(:,:)
      real(dp), pointer :: S_tmp(:)
      complex(dp), pointer  :: H_so_off(:,:)
      real(dp), pointer :: H_dftu(:,:)
      complex(dp), pointer  :: H_dftu_so(:,:)
#ifdef FINAL_CHECK_HS
      real(dp) :: diff_H
#endif
#ifdef SIESTA__WANNIER90
      real(dp), pointer   :: H_chempotwann(:,:)
#endif
#ifdef MPI
      real(dp)              :: stresstmp(3,3)
      real(dp), pointer     :: fatmp(:,:)
      real(dp)              :: buffer1
#endif
      ! Just make the filename long enough!
      character(len=label_length+56) :: fname
      integer :: fc_istep, fc_ia
      logical :: not_using_auxcell, dummy_logical
      logical :: lwrite_hsx, lwrite_tshs, lwrite_hs
#ifdef CDF
#ifdef NCDF_4
      type(dictionary_t) :: dic_save
#endif
#endif


!------------------------------------------------------------------------- BEGIN
      call timer( 'FINAL_HF', 1 )

      not_using_auxcell = (no_s == no_u)

!     Initialize Hamiltonian ........................................
      nullify( H_new )
      call re_alloc(H_new, 1,maxnh, 1, spin%H,
     $             'H_new','final_H_f_stress')

      H_new(:,:) = 0.0_dp

!     Initialize forces and stress ...................
      nullify(fal)
      call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'final_H_f_stress' )

      fa(1:3,1:na_u)   = 0.0_dp
      fal(1:3,1:na_u)  = 0.0_dp
      stress(1:3,1:3)  = 0.0_dp
      stressl(1:3,1:3) = 0.0_dp


! Neutral-atom: energy, forces and stress ............................
      ! These routines parallelize fa and stress over atoms
      call naefs( na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, Ena, fal, stressl, forces_and_stress=.true. )
      call dnaefs( na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, DEna, fal, stressl, forces_and_stress=.true. )
      Ena = Ena + DEna

!     Kinetic: energy, forces, stress and matrix elements .................
      ! Forces and stresses parallelized over orbitals
      call kinefsm( na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     &              maxnh, maxnh, lasto, iphorb, isa,
     &              numh, listhptr, listh, numh, listhptr, listh,
     &              spin%spinor, Dscf, Ekin, fal, stressl, H_new,
     &              matrix_elements_only=.false. )
#ifdef MPI
! Global reduction of energy terms
      call globalize_sum(Ekin,buffer1)
      Ekin = buffer1
#endif
! ..................

! Non-local-pseudop: energy, forces, stress and matrix elements .......

      Eso = 0.0_dp

      if ( .not.spin%SO_offsite ) then
      ! Forces and stresses parallelized over orbitals
       call nlefsm( scell, na_u, na_s, isa, xa, indxua,
     &              maxnh, maxnh, lasto, lastkb, iphorb, iphKB,
     &              numh, listhptr, listh, numh, listhptr, listh,
     &              spin%spinor, Dscf, Enl, fal, stressl, H_new,
     &              matrix_elements_only=.false. )

#ifdef MPI
! Global reduction of energy terms
       call globalize_sum(Enl,buffer1)
       Enl = buffer1
#endif
      else
        H_so_off => val(H_so_off_2D)
      ! Forces and stresses parallelized over orbitals
       call nlefsm_SO_off(scell, na_u, na_s, isa, xa, indxua,
     &                maxnh, maxnh, lasto, lastkb, iphorb, iphKB,
     &                numh, listhptr, listh, numh, listhptr, listh,
     &                spin%Grid,
     &                Enl, Eso, fal,
     &                stressl, H_new, H_so_off,
     &                matrix_elements_only=.false.,
     $                ill_defined_sr_so_split=dummy_logical)
#ifdef MPI
! Global reduction of energy terms
       call globalize_sum(Enl,buffer1)
       Enl = buffer1
       call globalize_sum(Eso,buffer1)
       Eso = buffer1
#endif
      endif

      if ( spin%SO_onsite ) then
         call spinorb(no_u,no_l,iaorb,iphorb,isa,indxuo,
     &        maxnh,numh,listhptr,listh,Dscf,H_new(:,3:),Eso)
      endif

!     Non-SCF part of total energy
      call update_E0()

! Hubbard term for LDA+U: energy, forces, stress and matrix elements ....
      if( switch_dftu ) then
        if ( spin%SO ) then
          H_dftu_so => val(H_dftu_so_2D)
          ! Forces and stresses parallelized over orbitals
          call dftu_so_hamil( H_dftu_so, fal, stressl )
#ifdef MPI
          ! Global reduction of energy terms
          call globalize_sum(E_dftu_so,buffer1)
          E_dftu_so = buffer1
          ! Global reduction of E_correc_dc
          call globalize_sum(E_correc_dc,buffer1)
          E_correc_dc = buffer1
#endif
        else
          H_dftu => val(H_dftu_2D)
          ! Forces and stresses parallelized over orbitals
          call hubbard_term(scell, na_u, na_s, isa, xa, indxua,
     .         maxnh, maxnh, lasto, iphorb, no_u, no_l,
     .         numh, listhptr, listh, numh, listhptr, listh,
     .         spin%spinor, Dscf, Edftu, DEdftu, H_dftu,
     .         fal, stressl, H_new, iscf,
     .         matrix_elements_only=.false.)
#ifdef MPI
          ! Global reduction of energy terms
          call globalize_sum(Edftu,buffer1)
          Edftu = buffer1
          ! DEdftu should not be globalized as it
          ! is based on globalized occupations
#endif
          Edftu = Edftu + DEdftu
        endif
      endif
! ..................


!     Non-local-pseudop: energy, forces, stress and matrix elements
!     Add SCF contribution to energy and matrix elements
      g2max = g2cut
!     Last call to dhscf and grid-cell sampling if requested
      ifa  = 1
      istr = 1

      ! This will call dhscf with the final DM coming out of
      ! the scf cycle

      call grdsam( spin%Grid, no_s, iaorb, iphorb,
     &             no_l, no_u, na_u, na_s, isa, xa, indxua,
     &             ucell, ntm, ifa, istr, maxnh,
     &             maxnh, numh, listhptr, listh, Dscf, Datm, H_new,
     &             Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     &             Exc, Dxc, dipol, stress, fal, stressl )

      ! Orthonormalization forces
      call overfsm( na_u, na_s, no_l, no_s, scell, xa, indxua, rmaxo,
     &              lasto, iphorb, isa, maxnh, numh, listhptr, listh,
     &              spin, Escf, fal, stressl )

! Energy term coming from the chemical potential added to the Wannier functions
#ifdef SIESTA__WANNIER90
      if( w90_wrapper_processing ) then
        if( compute_chempotwann_after_scf ) then
          H_chempotwann => val(H_chempotwann_2D)
          call add_Hamiltonian_chempotwann( H_chempotwann )
        endif
      endif
#endif

!     Metadynamics forces
      if (lMetaForce) then
        call meta(xa,na_u,ucell,Emeta,fa,stress,.true.,.true.)
      endif

!     Add on force field contribution if required
      call twobody( na_u, xa, isa, ucell, Emm,
     &              ifa=1, fa=fa, istr=1, stress=stress )

      Evdw_d3 = 0.0_dp
#ifdef SIESTA__DFTD3
!     Add D3 corrections to forces and energy (if required) Note that
!     only the master node does the computations.  Its 'fal' and
!     'stressl' will be updated. The rest of the nodes will not add
!     anything to these arrays.

      if (want_dftd3_dispersion) then
         if (IOnode) then
            call dftd3_energy_forces( na_u, xa, iza, ucell,
     &                                cell_periodic, Evdw_d3, fal,
     &                                stressl )
         endif
#ifdef MPI
         ! Global reduction of D3 energy terms
         call globalize_sum( Evdw_d3, buffer1 )
         Evdw_d3 = buffer1
#endif
      endif
#endif

#ifdef MPI
!     Global reduction of forces and stresses
      nullify(fatmp)
      call re_alloc( fatmp, 1, 3, 1, na_u, 'fatmp', 'final_H_f_stress' )
      call globalize_sum(stressl(1:3,1:3),stresstmp(1:3,1:3))
      call globalize_sum( fal(1:3,1:na_u), fatmp(1:3,1:na_u) )
      stress(1:3,1:3) = stress(1:3,1:3) + stresstmp(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fatmp(1:3,1:na_u)
      call de_alloc( fatmp, 'fatmp', 'final_H_f_stress' )
#else
      stress(1:3,1:3) = stress(1:3,1:3) + stressl(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fal(1:3,1:na_u)
#endif

      ! If backward compatibility is requested, recover
      ! the old behavior with respect to H

#ifdef FINAL_CHECK_HS
      ! We also print-out the max-difference between the
      ! final H after SCF and the final computed H used for
      ! the forces.
      ! This value should correspond to dHmax
      call compute_max_diff(H, H_new, diff_H)
      if ( IONode ) then
         ! Print out the final differences
         write(*,"(a,f10.6)") ":!: Final H_scf - H_force (eV) : ",
     &        diff_H/eV
      end if
#endif

      ! TODO I am not sure this works with SO_offsite?
      if (recompute_H_after_scf) then
         if (ionode) then
            write(6,"(a)") ":!: Updating H after scf cycle" //
     $                     " as requested by compat. option"
         endif
         H = H_new
      endif

      call de_alloc( fal, 'fal', 'final_H_f_stress' )
      call de_alloc( H_new, 'H_new', 'final_H_f_stress' )

      ! Determine whether anything should be saved
      ! If not, return immediately
      if ( .not. SCFconverged ) then
        call timer( 'FINAL_HF', 2 )
        return
      endif

! Save Hamiltonian and overlap matrices ............................

! Only in HSX format now.  Use Util/HSX/hsx2hs to generate an HS file
      lwrite_hsx = savehs .or. write_coop .or. write_hs_history
      lwrite_tshs = TS_HS_save .or. write_tshs_history
      ! Whether *any* HS file should be written.
      lwrite_hs = lwrite_hsx .or. lwrite_tshs

#ifdef CDF
#ifdef NCDF_4
      if ( idyn == 6 ) then
         if ( fc_save_dhs ) then
           if (istep == 0) then
             call init_dHSdR_nc(trim(slabel)//'.dHSdR.nc', no_u,
     &                          spin%H)
           else
             call io_fc_indices(istep, ia1, fc_ia, fc_istep)
             select case (fc_istep)
             case(2,4,6)
               fc_H_2D(2) = H_2D
               fc_S_1D(2) = S_1D
               fc_nsc(:, 2) = nsc(:)
               call calc_dHSdR(fc_ia, fc_istep/2,
     &                         fc_H_2D, dHdR_2D,
     &                         fc_S_1D, dSdR_1D)
               ! Clean-up before write
               call delete(fc_H_2D(1))
               !call delete(fc_H_2D(2))
               call delete(fc_S_1D(1))
               !all delete(fc_S_1D(2))
               call write_dHSdR_nc(trim(slabel)//'.dHSdR.nc', fc_ia,
     &                             fc_istep/2, dHdR_2D, dSdR_1D,
     &                             product(nsc))
               ! No need to carry around these extra things
               ! It's extra memory for no reason.
               call delete(dHdR_2D)
               call delete(dSdR_1D)
               case(1,3,5)
                 fc_nsc(:, 1) = nsc(:)
                 call copy_H_m(H_2D,fc_H_2D(1),S_1D,fc_S_1D(1))
              end select
           end if
         end if
      end if

      if ( write_cdf ) then
        if ( idyn == 6) then
          ! Store FC forces
          call cdf_save_fc(trim(slabel)//'.nc', istep)

          if ( istep == 0 ) then
            dic_save = ('fa'.kv.1)//('stress'.kv.1)//('Ef'.kv.1)
            dic_save = dic_save//('NEGF_Vha'.kv.1)
            if ( lwrite_hs ) dic_save = dic_save//('H'.kv.1)
            call cdf_save_state(trim(slabel)//'.nc',dic_save)
            call delete(dic_save)
          end if

        else
          dic_save = ('fa'.kv.1)//('stress'.kv.1)//('Ef'.kv.1)
          dic_save = dic_save//('NEGF_Vha'.kv.1)
          if ( lwrite_hs ) dic_save = dic_save//('H'.kv.1)
          call cdf_save_state(trim(slabel)//'.nc',dic_save)
          call delete(dic_save)
        end if
      end if
#endif
#endif

      if ( lwrite_hs .and. fixspin ) then
! For fixed spin calculations we shift the Hamiltonian according
! to the first spin such that we have a "common" E_F, after writing
! we shift back again.
         H_tmp => val(H_2D)
         S_tmp => val(S_1D)
         call daxpy( size(S_tmp), Efs(1)-Efs(2), S_tmp(1), 1,
     &               H_tmp(1,2), 1 )
         Ef = Efs(1)
      end if
!
! Note that we save the Hamiltonian coming out of the scf cycle,
! not the one computed (from DM_out) in this routine.
! This call could be moved to a more appropriate place
!
      if ( lwrite_hsx ) then
        block
           type(kpoint_t), pointer :: kpt
           character(len=8) :: ext
           if ( TSrun ) then
              ext = "TS.HSX"
              kpt => ts_kpoint_scf
           else
              ext = "HSX"
              kpt => kpoint_scf
           end if

           if ( write_hs_history ) then
              if ( idyn == 6 ) then
                 call io_ext_fc(ext, istep, ia1, fname)
              else
                 call io_ext_md(ext, istep, fname)
              end if
           else
              fname = trim(ext)
           end if

           call write_hsx(H_2D, S_1D, Ef, Qtot, Temp, kpt,
     &          extension=fname)
        end block
      endif

      ! Initialize output variables for a *single* write call
      fc_istep = 0
      fc_ia = 0
      if ( lwrite_tshs ) then
         if ( write_tshs_history ) then
            if ( idyn == 6 ) then
               ! Correct the FC indices such that
               ! the TSHS format contains, atom displaced, index of
               ! displacement:
               !  1 = -x
               !  2 = +x
               !  3 = -y
               !  4 = +y
               !  5 = -z
               !  6 = +z
               ! The TSHS fileformat stores the displaced atom + the
               ! step of the displacement.
               call io_fc_indices(istep, ia1, fc_ia, fc_istep)
               call io_ext_fc("TSHS", istep, ia1, fname)
            else
               fc_istep = istep ! ia is still 0
               call io_ext_md("TSHS", istep, fname)
            end if
            fname = trim(slabel) // "." // trim(fname)
         else
            fname = trim(slabel) // ".TSHS"
         end if

         call ts_write_tshs(fname, .false.,
     &        not_using_auxcell, ts_kpoint_scf%is_gamma(),
     &        ucell, nsc, isc_off, na_u, no_s, spin%H,
     &        ts_kpoint_scf%k_cell, ts_kpoint_scf%k_displ,
     &        xa, lasto,
     &        H_2D, S_1D, indxuo,
     &        Ef, Qtot, Temp, fc_istep, fc_ia)
      end if

      if ( lwrite_hs .and. fixspin ) then
        ! Shift back
        call daxpy( size(S_tmp), Efs(2)-Efs(1), S_tmp(1), 1,
     &              H_tmp(1,2), 1 )
      end if

      call timer( 'FINAL_HF', 2 )
!----------------------------------------------------------------------- END
      END subroutine final_H_f_stress
      END module m_final_H_f_stress
