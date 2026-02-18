
#include "mpi_macros.f"

module dHSdR_m
   use precision, only: dp
   use class_OrbitalDistribution
   use class_Sparsity
   use class_dSpData1D
   use class_dSpData2D
   use m_spin, only: spin
#ifdef NCDF_4
   use netcdf_ncdf
   use ncdf_io_m
#endif

   implicit none

   private

   public :: copy_H_m
   public :: calc_dHSdR
#ifdef NCDF_4
   public :: init_dHSdR_nc
   public :: write_dHSdR_nc
#endif

contains

   subroutine copy_H_m(H_in_2D, H_out_2D, S_in_1D, S_out_1D)
         !! Copy Hamiltonian for current displacement (should be negative) to be used for calculating the
         !! derivative later on
      implicit none

      type(dSpData2D), intent(in) :: H_in_2D
         !! Hamiltonian of current SCF step
      type(dSpData2D), intent(inout) :: H_out_2D
         !! Copy of Hamiltonian of current SCF step
      type(dSpData1D), intent(in) :: S_in_1D
         !! Overlap matrix of current geometry step
      type(dSpData1D), intent(inout) :: S_out_1D
         !! Copy of overlap matrix of current geometry step

      type(Sparsity) :: sp_out
      type(Sparsity), pointer :: sp
      type(OrbitalDistribution), pointer :: block_dist
      real(dp), pointer :: H_in(:, :) => null(), H_out(:, :) => null()
      real(dp), pointer :: S_in(:) => null(), S_out(:) => null()

      ! Copy sparsity
      sp => spar(H_in_2D)
      block_dist => dist(H_in_2D)

      call newSparsity(sp_out, nrows(sp), nrows_g(sp), nnzs(sp), n_col(sp), list_ptr(sp), &
                       list_col(sp), name=trim(sp%data%name))

      ! Create new sparse matrices
      call newdSpData2D(sp_out, spin%H, block_dist, H_out_2D, name=trim(H_in_2D%data%name))
      call newdSpData1D(sp_out, block_dist, S_out_1D, name=trim(S_in_1D%data%name))

      ! Copy over the values
      H_in => val(H_in_2D)
      H_out => val(H_out_2D)
      call dcopy(nnzs(sp)*spin%H, H_in, 1, H_out, 1)

      S_in => val(S_in_1D)
      S_out => val(S_out_1D)
      call dcopy(nnzs(sp), S_in, 1, S_out, 1)

      ! clean-up, now there should only be 2 references (H|S_out_2|1D)
      call delete(sp_out)

   end subroutine copy_H_m

   subroutine calc_dHSdR(ia, ixyz, fc_H_2D, dHdR_2D, fc_S_1D, dSdR_1D)
      !! Calculate the derivative of H / S from negative and positive displacments

      implicit none

      integer, intent(in) :: ia
         !! Index of displaced atom
      integer, intent(in) :: ixyz
         !! Index of axis along which atom is displaced
      type(dSpData2D), intent(in) :: fc_H_2D(2)
         !! Hamiltonian for positive and negative displacement
      type(dSpData2D), intent(inout) :: dHdR_2D
         !! Derivatigve of Hamiltonian
      type(dSpData1D), intent(in) :: fc_S_1D(2)
         !! Overlap matrix for positive and negative displacement
      type(dSpData1D), intent(inout) :: dSdR_1D
         !! Derivatigve of overlap matrix

      call calc_dHdR(ia, ixyz, fc_H_2D, dHdR_2D)
      call calc_dSdR(ia, ixyz, fc_S_1D, dSdR_1D)

   end subroutine calc_dHSdR

   subroutine calc_dHdR(ia, ixyz, fc_H_2D, dHdR_2D)
      !! Calculate the derivative of H / S from negative and positive displacements
      use intrinsic_missing, only: index_sort_heap
      use siesta_options, only: fc_dHdR_tol, dx
      use fc_matrices_m, only: fc_nsc

      implicit none

      integer, intent(in) :: ia
         !! Index of displaced atom
      integer, intent(in) :: ixyz
         !! Index of axis along which atom is displaced
      type(dSpData2D), intent(in) :: fc_H_2D(2)
         !! Hamiltonian for positive and negative displacement
      type(dSpData2D), intent(inout) :: dHdR_2D
         !! Derivative of Hamiltonian

      type(OrbitalDistribution), pointer :: block_dist
      type(Sparsity), pointer :: sp_m, sp_p
      type(Sparsity) :: sp_dH
      real(dp), pointer :: H_m(:, :), H_p(:, :), dH(:, :), dHdR(:,:)
      integer, pointer :: ncol_m(:), ncol_p(:), ncol(:)
      integer, pointer :: lptr_m(:), lptr_p(:), lptr(:)
      integer, pointer :: lcol_m(:), lcol_p(:), lcol(:)
      integer :: no_l, mo_l, no_u, mo_u, n_nzs, nnzs_m, nnzs_p
      integer :: io, im, ip, iim, iip, ind, ind_m, ind_p, ispin, jm, jp
      integer, allocatable :: sidx_m(:), sidx_p(:), indices(:, :)
      logical :: discard
      real(dp) :: d2R, dH_tol
      character(len=150) :: sname

      d2R = 2._dp * dx
      dH_tol = fc_dHdR_tol * d2R

      block_dist => dist(fc_H_2D(1))
      sp_m => spar(fc_H_2D(1))
      sp_p => spar(fc_H_2D(2))

      if (.not. all(fc_nsc(:, 1) == fc_nsc(:, 2))) then
         call die("Super cell changed during FC displacement.  &
                  &Unable to calculate derivative.")
      end if

      call attach(sp_m, list_ptr=lptr_m, list_col=lcol_m, nrows=mo_l, nrows_g=mo_u, n_col=ncol_m, nnzs=nnzs_m)
      call attach(sp_p, list_ptr=lptr_p, list_col=lcol_p, nrows=no_l, nrows_g=no_u, n_col=ncol_p, nnzs=nnzs_p)
      if ((mo_l /= no_l) .or. (mo_u /= no_u)) then
         call die('calc_dHSdR: sparsity dimensions for positive and negative displacements '&
                 &'do not match')
      end if
      H_m => val(fc_H_2D(1))
      H_p => val(fc_H_2D(2))

      allocate (sidx_m(maxval(ncol_m)))
      allocate (sidx_p(maxval(ncol_p)))
      allocate (indices(2, nnzs_m + nnzs_p))
      allocate (ncol(no_l))
      allocate (lptr(no_l))

      ! Determine nnzs for derivative, and the mapping between the sparse matrices
      allocate (dH(spin%H, 1))
      ind = 0
      do io = 1, no_l
         im = 1
         ip = 1

         ! Sorted indices in current row
         sidx_m(:) = -1
         sidx_p(:) = -1
         call index_sort_heap(ncol_m(io), lcol_m(lptr_m(io) + 1:lptr_m(io) + ncol_m(io)), sidx_m)
         call index_sort_heap(ncol_p(io), lcol_p(lptr_p(io) + 1:lptr_p(io) + ncol_p(io)), sidx_p)

         lptr(io) = ind
         do while ((im <= ncol_m(io)) .and. (ip <= ncol_p(io)))

            iim = sidx_m(im)
            iip = sidx_p(ip)

            ! Current column indices in both sparsities
            jm = lcol_m(lptr_m(io) + iim)
            jp = lcol_p(lptr_p(io) + iip)

            ! Check size of matrix elements of derivative
            ind = ind + 1
            if (jm == jp ) then ! Matrix element is non-zero in both matrices
               indices(1, ind) = lptr_m(io) + iim
               indices(2, ind) = lptr_p(io) + iip

               dH(:, 1) = H_p(lptr_p(io) + iip, :) - H_m(lptr_m(io) + iim, :)

               im = im + 1
               ip = ip + 1
            else if (jm < jp) then ! Matrix element only exists for the negative displacement
               indices(1, ind) = lptr_m(io) + iim
               indices(2, ind) = -1

               dH(:, 1) = -H_m(lptr_m(io) + iim, :)
               im = im + 1
            else ! Matrix element only exists for the positive displacement
               indices(1, ind) = -1
               indices(2, ind) = lptr_p(io) + iip
               dH(:, 1) = H_p(lptr_p(io) + iip, :)
               ip = ip + 1
            end if

            ! Check size of matrix elements of derivative
            discard = .false.
            do ispin = 1, spin%H
               discard = discard .and. (abs(dH(ispin, 1)) < dH_tol)
            end do

            if (discard) then
               indices(:, ind) = -1
               ind = ind - 1
            end if
         end do ! im / ip

         ! At this point we are finished with the columns of at least one of the
         ! two sparsities, but we might have to finish parsing the other one...
         do while ((im <= ncol_m(io)))

            iim = sidx_m(im)

            ! Current column index
            jm = lcol_m(lptr_m(io) + iim)

            ind = ind + 1
            indices(1, ind) = lptr_m(io) + iim
            indices(2, ind) = -1

            dH(:, 1) = -H_m(lptr_m(io) + iim, :)
            im = im + 1

            discard = .false.
            do ispin = 1, spin%H
               discard = discard .and. (abs(dH(ispin, 1)) < dH_tol)
            end do

            if (discard) then
               indices(:, ind) = -1
               ind = ind - 1
            end if

         end do ! im

         do while ((ip <= ncol_p(io)))

            iip = sidx_p(ip)

            ! Current column index
            jp = lcol_p(lptr_p(io) + iip)

            ind = ind + 1
            indices(1, ind) = -1
            indices(2, ind) = lptr_p(io) + iip

            dH(:, 1) = H_p(lptr_p(io) + iip, :)
            ip = ip + 1

            discard = .false.
            do ispin = 1, spin%H
               discard = discard .and. (abs(dH(ispin, 1)) < dH_tol)
            end do

            if (discard) then
               indices(:, ind) = -1
               ind = ind - 1
            end if

         end do ! ip

         ncol(io) = ind - lptr(io)

      end do ! io

      n_nzs = ind

      deallocate (dH)
      allocate (lcol(n_nzs))
      allocate (dH(n_nzs, spin%H))
      dH(:,:) = 0._dp
      do ind = 1, n_nzs
         ind_m = indices(1, ind)
         ind_p = indices(2, ind)
         if (ind_p > 0) then
            lcol(ind) = lcol_p(ind_p)
            do ispin = 1, spin%H
               dH(ind, ispin) = H_p(ind_p, ispin)
            end do
         end if
         if (ind_m > 0) then
            lcol(ind) = lcol_m(ind_m)
            do ispin = 1, spin%H
               dH(ind, ispin) = dH(ind, ispin)-H_m(ind_m, ispin)
            end do
         end if
      end do

      deallocate (indices, sidx_m, sidx_p)

      write (sname, '(a,I0,a,I0)') &
         'Hamiltonian derivative for displacment of atom ', ia, 'along ', ixyz
      call newSparsity(sp_dH, no_l, no_u, n_nzs, ncol, lptr, lcol, name=trim(sname))

      deallocate (ncol, lptr, lcol)

      write (sname, '(a,I0,a,I0)') 'Hamiltonian derivative for displacment of atom ', ia, 'along ', ixyz
      call newdSpData2D(sp_dH, spin%H, block_dist, dHdR_2D, name=trim(sname))
      dHdR => val(dHdR_2D)
      dHdR(:,:) = dH(:,:) / d2R

      deallocate(dH)

   end subroutine calc_dHdR

   subroutine calc_dSdR(ia, ixyz, fc_S_1D, dSdR_1D)
      !! Calculate the derivative of H / S from negative and positive displacments
      use intrinsic_missing, only: index_sort_heap
      use siesta_options, only: fc_dSdR_tol, dx
      use fc_matrices_m, only: fc_nsc
      use atomlist, only: iaorb, indxuo

      implicit none

      integer, intent(in) :: ia
         !! Index of displaced atom
      integer, intent(in) :: ixyz
         !! Index of axis along which atom is displaced
      type(dSpData1D), intent(in) :: fc_S_1D(2)
         !! Overlap matrix for positive and negative displacement
      type(dSpData1D), intent(inout) :: dSdR_1D
         !! Derivative of overlap matrix

      type(OrbitalDistribution), pointer :: block_dist
      type(Sparsity), pointer :: sp_m, sp_p
      type(Sparsity) :: sp_dS
      real(dp), pointer :: S_m(:), S_p(:), dS(:), dSdR(:)
      integer, pointer :: ncol_m(:), ncol_p(:), ncol(:)
      integer, pointer :: lptr_m(:), lptr_p(:), lptr(:)
      integer, pointer :: lcol_m(:), lcol_p(:), lcol(:)
      integer :: no_l, mo_l, no_u, mo_u, n_nzs, nnzs_m, nnzs_p
      integer :: io, i_g, im, ip, iim, iip, ind, ind_m, ind_p, jm, jp, juo
      integer, allocatable :: sidx_m(:), sidx_p(:), indices(:, :)
      real(dp) :: d2R, dS_tol
      character(len=150) :: sname

      d2R = 2._dp * dx
      dS_tol = fc_dSdR_tol * d2R

      block_dist => dist(fc_S_1D(1))
      sp_m => spar(fc_S_1D(1))
      sp_p => spar(fc_S_1D(2))

      if (.not. all(fc_nsc(:, 1) == fc_nsc(:, 2))) then
         call die("Super cell changed during FC displacement.  &
                  &Unable to calculate derivative.")
      end if

      call attach(sp_m, list_ptr=lptr_m, list_col=lcol_m, nrows=mo_l, nrows_g=mo_u, n_col=ncol_m, nnzs=nnzs_m)
      call attach(sp_p, list_ptr=lptr_p, list_col=lcol_p, nrows=no_l, nrows_g=no_u, n_col=ncol_p, nnzs=nnzs_p)
      if ((mo_l /= no_l) .or. (mo_u /= no_u)) then
         call die('calc_dHSdR: sparsity dimensions for positive and negative displacements '&
                 &'do not match')
      end if
      S_m => val(fc_S_1D(1))
      S_p => val(fc_S_1D(2))

      allocate (sidx_m(maxval(ncol_m, 1)))
      allocate (sidx_p(maxval(ncol_p, 1)))
      allocate (indices(2, nnzs_m + nnzs_p))
      allocate (ncol(no_l))
      allocate (lptr(no_l))

      ! Determine nnzs for derivative, and the mapping between the sparse matrices
      allocate (dS(1))
      ind = 0
      do io = 1, no_l
         i_g = index_local_to_global(block_dist, io)
         im = 1
         ip = 1

         ! Sorted indices in current row
         sidx_m(:) = -1
         sidx_p(:) = -1
         call index_sort_heap(ncol_m(io), lcol_m(lptr_m(io) + 1:lptr_m(io) + ncol_m(io)), sidx_m)
         call index_sort_heap(ncol_p(io), lcol_p(lptr_p(io) + 1:lptr_p(io) + ncol_p(io)), sidx_p)

         lptr(io) = ind
         do while ((im <= ncol_m(io)) .and. (ip <= ncol_p(io)))

            iim = sidx_m(im)
            iip = sidx_p(ip)

            ! Current column indices in both sparsities
            jm = lcol_m(lptr_m(io) + iim)
            jp = lcol_p(lptr_p(io) + iip)
            juo = indxuo(min(jp,jm))


            ind = ind + 1
            if (jm == jp ) then ! Matrix element is non-zero in both matrices
               indices(1, ind) = lptr_m(io) + iim
               indices(2, ind) = lptr_p(io) + iip
               dS(1) = S_p(lptr_p(io) + iip) - S_m(lptr_m(io) + iim)
               im = im + 1
               ip = ip + 1
            else if (jm < jp) then ! Matrix element only exists for the negative displacement
               indices(1, ind) = lptr_m(io) + iim
               indices(2, ind) = -1
               dS(1) = -S_m(lptr_m(io) + iim)
               im = im + 1
            else ! Matrix element only exists for the positive displacement
               indices(1, ind) = -1
               indices(2, ind) = lptr_p(io) + iip
               dS(1) = S_p(lptr_p(io) + iip)
               ip = ip + 1
            end if

            ! Only consider elements where left or right orbital was moved and
            ! Check size of matrix elements of derivative
            if ((ia/=iaorb(i_g)).and.(ia/=iaorb(juo)) .or. &
                (abs(dS(1)) < dS_tol)) then
               indices(:, ind) = -1
               ind = ind - 1
            end if
         end do ! im / ip

         ! At this point we are finished with the columns of at least one of the
         ! two sparsities, but we might have to finish parsing the other one...
         do while ((im <= ncol_m(io)))
            iim = sidx_m(im)
            ! Current column index
            jm = lcol_m(lptr_m(io) + iim)
            juo = indxuo(jm)

            ind = ind + 1
            indices(1, ind) = lptr_m(io) + iim
            indices(2, ind) = -1
            dS(1) = -S_m(lptr_m(io) + iim)
            im = im + 1
            if (((ia/=iaorb(i_g)).and.(ia/=iaorb(juo))) .or. &
                (abs(dS(1)) < dS_tol)) then
               indices(:, ind) = -1
               ind = ind - 1
            end if

         end do ! im

         do while ((ip <= ncol_p(io)))

            iip = sidx_p(ip)

            ! Current column index
            jp = lcol_p(lptr_p(io) + iip)
            juo = indxuo(jp)

            ind = ind + 1
            indices(1, ind) = -1
            indices(2, ind) = lptr_p(io) + iip
            dS(1) = S_p(lptr_p(io) + iip)
            ip = ip + 1
            if ((ia/=iaorb(i_g)).and.(ia/=iaorb(juo)) .or. &
                (abs(dS(1)) < dS_tol)) then
               indices(:, ind) = -1
               ind = ind - 1
            end if

         end do ! ip

         ncol(io) = ind - lptr(io)

      end do ! io

      n_nzs = ind

      deallocate (dS)
      allocate (lcol(n_nzs))
      allocate (dS(n_nzs))
      dS(:) = 0._dp
      do ind = 1, n_nzs
         ind_m = indices(1, ind)
         ind_p = indices(2, ind)
         if (ind_p > 0) then
            lcol(ind) = lcol_p(ind_p)
            dS(ind) = S_p(ind_p)
         end if
         if (ind_m > 0) then
            lcol(ind) = lcol_m(ind_m)
            dS(ind) = dS(ind)-S_m(ind_m)
         else
            ! Do nothing
            call die('calc_dHSdR: programming error')
         end if
      end do

      deallocate (indices, sidx_m, sidx_p)

      write (sname, '(a,I0,a,I0)') &
         'Overlap matrix derivative for displacment of atom ', ia, 'along ', ixyz
      call newSparsity(sp_dS, no_l, no_u, n_nzs, ncol, lptr, lcol, name=trim(sname))

      deallocate (ncol, lptr, lcol)

      write (sname, '(a,I0,a,I0)') 'Derivative overlap matrix for displacment of atom ', ia, 'along ', ixyz
      call newdSpData1D(sp_dS, block_dist, dSdR_1D, name=trim(sname))
      dSdR => val(dSdR_1D)
      dSdR(:) = dS(:) / d2R

      deallocate(dS)

   end subroutine calc_dSdR

#ifdef NCDF_4
   subroutine init_dHSdR_nc(fname, no_u, nspin)
      !! Initialize netCDF file storing the derivatives of the Hamiltonian and overlap matrix
      use dictionary
      use siesta_options, only: fc_dHdR_tol, fc_dSdR_tol, dx, ia1, ia2
      use siesta_geom, only: na_u

      implicit none

      character(len=*), intent(in) :: fname
         !! Filename
      integer, intent(in) :: no_u
         !! Number of unit cell orbitals
      integer, intent(in) :: nspin
         !! Number of spin components in H

      ! Internal
      type(dictionary_t) :: dic
      type(hNCDF) :: ncfile, grp
      integer :: n_disp, ia

      call ncdf_create(ncfile,fname, mode=ior(NF90_WRITE,NF90_NETCDF4), overwrite=.true.)

      dic = ('DIMno_u'.kv.no_u)//('DIMspin'.kv.nspin)//('DIMna_u'.kv.na_u)
      dic = dic//('DIMxyz'.kv.3)//('DIMone'.kv.1)
      call ncdf_crt(ncfile, dic)
      call delete(dic)

      ! call ncdf_def_var(ncfile, 'isa', NF90_DOUBLE, (/'na_u'/), atts=dic)
      ! call ncdf_put_var(ncfile, 'isa', isa)

      dic = dic//('info'.kv.'Tolerance for dH/dR matrix elements')
      dic = dic//('unit'.kv.'Ry/Bohr')
      call ncdf_def_var(ncfile, 'dHdR.Tolerance', NF90_DOUBLE, (/'one'/), atts=dic)
      call ncdf_put_var(ncfile, 'dHdR.Tolerance', fc_dHdR_tol)

      dic = dic//('info'.kv.'Tolerance for dH/dR matrix elements')
      dic = dic//('unit'.kv.'1/Bohr')
      call ncdf_def_var(ncfile, 'dSdR.Tolerance', NF90_DOUBLE, (/'one'/), atts=dic)
      call ncdf_put_var(ncfile, 'dSdR.Tolerance', fc_dSdR_tol)

      dic = dic//('info'.kv.'Length of FC displacement')
      dic = dic//('unit'.kv.'Bohr')
      call ncdf_def_var(ncfile, 'FC.Displacement', NF90_DOUBLE, (/'one'/), atts=dic)
      call ncdf_put_var(ncfile, 'FC.Displacement', dx)
      call delete(dic)

      call ncdf_def_grp(ncfile, "DISPLACEMENTS", grp)
      n_disp = (ia2 - ia1 + 1)
      dic = ('DIMn_disp'.kv.n_disp)
      call ncdf_crt(grp, dic)
      call delete(dic)

      dic = dic//('info'.kv.'Indices of displaced atoms')
      call ncdf_def_var(grp, 'atom_index', NF90_INT, (/'n_disp'/), atts=dic)
      call ncdf_put_var(grp, 'atom_index', (/(ia, ia=ia1,ia2)/))

      call ncdf_close(ncfile)
   end subroutine init_dHSdR_nc
#endif

#ifdef NCDF_4
   subroutine write_dHSdR_nc(fname, ia, ixyz, dHdR_2D, dSdR_1D, n_s)
      !! Write the derivatives of H and S for the current displacement to nc file
      use siesta_options, only: cdf_comp_lvl, cdf_w_parallel
      use siesta_geom, only: nsc, isc_off
      use dictionary
#ifdef MPI
      use mpi_siesta, only: MPI_Reduce, MPI_Integer, MPI_SUM, MPI_Comm_World
#endif

      implicit none

      character(len=*), intent(in) :: fname
         !! File name
      integer :: ia
         !! Displaced atom
      integer :: ixyz
         !! Direction of displacement
      type(dSpData2D), intent(inout) :: dHdR_2D
         !! Derivative of Hamiltonian
      type(dSpData1D), intent(inout) :: dSdR_1D
         !! Derivative of overlap matrix
      integer, intent(in) :: n_s
         !! Number of super cells

      type(Sparsity), pointer :: sp
      type(OrbitalDistribution), pointer :: block_dist
      integer :: no, n_nzs, chks(3), tmp
      integer, allocatable :: gncol(:)
      type(dictionary_t) :: dic
      character(len=150) :: grpname
      type(hNCDF) :: ncfile, grp0, grp, grp_dH, grp_dS
#ifdef MPI
      integer :: MPIerror
#endif

#ifdef MPI
      if (cdf_w_parallel) then
         call ncdf_open(ncfile, fname, mode=ior(NF90_WRITE, NF90_MPIIO), &
                        comm=MPI_COMM_ID(MPI_COMM_WORLD))
      else
#endif
         call ncdf_open(ncfile, fname, mode=ior(NF90_WRITE, NF90_NETCDF4))
#ifdef MPI
      end if
#endif

      call ncdf_open_grp(ncfile, "DISPLACEMENTS", grp)

      write (grpname, '(I0)') ia
      if (ixyz == 1) then
         call ncdf_def_grp(grp, grpname, grp0)
      else
         call ncdf_open_grp(grp, grpname, grp0)
      end if

      write (grpname, '(I0)') ixyz
      call ncdf_def_grp(grp0, grpname, grp)

      call ncdf_def_dim(grp, 'n_s', n_s)

      dic = dic//('info'.kv.'Number of supercells in each unit-cell direction')
      call ncdf_def_var(grp,'nsc',NF90_INT,(/'xyz'/), atts=dic)

      dic = dic//('info'.kv.'Index of supercell coordinates')
      call ncdf_def_var(grp, 'isc_off', NF90_INT, (/'xyz', 'n_s'/), atts=dic)

      call ncdf_put_var(grp, 'nsc', nsc)
      call ncdf_put_var(grp, 'isc_off', isc_off)

      ! ------------------------------------------------------------------------
      ! Derivative of Hamiltonian
      ! ------------------------------------------------------------------------

      call ncdf_def_grp(grp, 'dH', grp_dH)

      sp => spar(dHdR_2D)
      block_dist => dist(dHdR_2D)

      no = nrows_g(sp)
#ifdef MPI
      tmp = nnzs(sp)
      call MPI_Reduce(tmp,n_nzs,1,MPI_Integer, MPI_Sum, 0, MPI_Comm_World, MPIerror)
#else
       n_nzs = nnzs(sp)
#endif
      allocate (gncol(no))
      ! Signal that it needs to be filled (first element is negative)
      gncol(1) = -1

      call ncdf_def_dim(grp_dH, 'nnzs', n_nzs)

      dic = dic//('info'.kv.'Number of non-zero elements per row')
      call ncdf_def_var(grp_dH, 'n_col', NF90_INT, (/'no_u'/), atts=dic)

      chks = (/n_nzs, 1, 1/)

      dic = dic//('info'.kv.'Supercell column indices in the sparse format')
      call ncdf_def_var(grp_dH, 'list_col', NF90_INT, (/'nnzs'/), compress_lvl=cdf_comp_lvl, &
                        atts=dic, chunks=chks)

      dic = dic//('info'.kv.'Derivative of Hamiltonian')
      dic = dic//('unit'.kv.'Ry')
      call ncdf_def_var(grp_dH, 'dH', NF90_DOUBLE, (/'nnzs', 'spin'/), compress_lvl=cdf_comp_lvl, &
                        atts=dic, chunks=chks)

      call cdf_w_Sp(grp_dH, sp, block_dist, gncol=gncol)
      call cdf_w_d2D(grp_dH, 'dH', dHdR_2D, gncol=gncol)
      call delete(dic)

      ! ------------------------------------------------------------------------
      ! Derivative of Overlap Matrix
      ! ------------------------------------------------------------------------

      call ncdf_def_grp(grp, 'dS', grp_dS)

      sp => spar(dSdR_1D)
      block_dist => dist(dSdR_1D)

      no = nrows_g(sp)
#ifdef MPI
      tmp = nnzs(sp)
      call MPI_Reduce(tmp,n_nzs,1,MPI_Integer, MPI_Sum, 0, MPI_Comm_World, MPIerror)
#else
       n_nzs = nnzs(sp)
#endif
      deallocate (gncol)
      allocate (gncol(no))
      ! Signal that it needs to be filled (first element is negative)
      gncol(1) = -1

      call ncdf_def_dim(grp_dS, 'nnzs', n_nzs)

      dic = dic//('info'.kv.'Number of non-zero elements per row')
      call ncdf_def_var(grp_dS, 'n_col', NF90_INT, (/'no_u'/), atts=dic)

      chks = (/n_nzs, 1, 1/)

      dic = dic//('info'.kv.'Supercell column indices in the sparse format')
      call ncdf_def_var(grp_dS, 'list_col', NF90_INT, (/'nnzs'/), compress_lvl=cdf_comp_lvl, &
                        atts=dic, chunks=chks)


      dic = ('info'.kv.'Derivative of overlap matrix')
      call ncdf_def_var(grp_dS, 'dS', NF90_DOUBLE, (/'nnzs'/), compress_lvl=cdf_comp_lvl, &
                        atts=dic, chunks=chks)

      call cdf_w_Sp(grp_dS, sp, block_dist, gncol=gncol)
      call cdf_w_d1D(grp_dS, 'dS', dSdR_1D, gncol=gncol)
      call delete(dic)

      call ncdf_close(ncfile)

   end subroutine write_dHSdR_nc
#endif

end module dHSdR_m
