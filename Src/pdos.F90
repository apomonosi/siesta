! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
subroutine pdos( NO, nspin, NO_L, MAXNH, &
    no_u, NUMH, LISTHPTR, LISTH, H, S, &
    E1, E2, SIGMA, NHIST, &
    XIJ, INDXUO, GAMMA, NK, KPOINT, WK)
  ! **********************************************************************
  ! Subroutine to calculate the projected density of states on the
  ! atomic orbitals for a given eigenvalue spectra
  ! Written by J. Junquera and E. Artacho, November 1999.
  ! ***********  INPUT  **************************************************
  ! INTEGER NO                  : Number of basis orbitals in the supercell
  ! integer nspin=h_spin_dim    : Number of spin components of H and D
  ! INTEGER NO_L                : Maximum number of atomic orbitals in the unit 
  !                               cell. First dimension of eo, qo, last of xij
  !                               Must be at least max(indxuo)
  !                               IN THIS PROCESSOR
  ! INTEGER MAXNH               : Maximum number of orbitals interacting 
  !                               with any orbital
  ! INTEGER no_u                : Number of total orbitals in unit cell
  ! INTEGER NUMH(NUO)           : Number of nonzero elements of each row
  !                               of hamiltonian matrix
  ! INTEGER LISTH(MAXNH)        : Nonzero hamiltonian-matrix element
  !                               column indexes for each matrix row
  ! INTEGER LISTHPTR(NUO)       : Pointer to each row (-1) of the
  !                               density matrix
  ! REAL*8  H(MAXNH,nspin) : Hamiltonian in sparse format
  ! REAL*8  S(MAXNH)            : Overlap in sparse format
  ! REAL*8  E1, E2              : Energy range for density-matrix states
  !                               (to find local density of states)
  !                               Not used if e1 > e2
  ! REAL*8  SIGMA               : Width of the gaussian to expand the eigenvalues
  ! INTEGER NHIST               : Number of subdivisions of the histogram
  ! REAL*8  XIJ(3,MAXNH)        : Vectors between orbital centers (sparse)
  !                               (not used if only gamma point)
  ! INTEGER INDXUO(NO)          : Index of equivalent orbital in unit cell
  !                               Unit cell orbitals must be the first in
  !                               orbital lists, i.e. indxuo.le.nuo, with
  !                               nuo the nuber of orbitals in the unit cell
  ! logical Gamma               : whether only the Gamma point is sampled
  ! INTEGER NK                  : Number of k points
  ! REAL*8  KPOINT(3,NK)        : k point vectors
  ! REAL*8  WK(NK)              : k point weights (must sum one)
  ! **********************************************************************

  use precision,    only : dp
  use parallel,     only : Node, Nodes, IOnode
  use fdf
  use siesta_geom,  only : xa, isa
  use m_spin,       only : spin
  use atomlist,     only : iphorb, iaorb
  use atmfuncs,     only : zetafio, mofio, lofio, cnfigfio, labelfis
  use atmfuncs,     only : pol, izofis
  use xmlf90_wxml,  only : xmlf_t, xml_OpenFile, xml_NewElement
  use xmlf90_wxml,  only : xml_AddArray, xml_AddPCData
  use xmlf90_wxml,  only : xml_AddAttribute, xml_EndElement, xml_Close
  use xml,          only : str, xml_dump_attribute
  use densematrix,  only : allocDenseMatrix, resetDenseMatrix
  use densematrix,  only : Haux, Saux, psi
  use alloc,        only : re_alloc, de_alloc
  use units,        only : eV
  use files,        only : slabel, label_length
  use m_energies, only: Ef
#ifdef MPI
  use parallelsubs, only : GetNodeOrbs
#endif
  use m_diag_option, only: ParallelOverK, Serial, print_diag

  implicit none

  integer :: NO, NSPIN, NO_L, MAXNH, NK, NHIST, NO_U

  logical, intent(in) :: Gamma
  integer :: NUMH(*), LISTH(MAXNH), LISTHPTR(*), INDXUO(NO)

  real(dp) :: H(MAXNH,NSPIN), S(MAXNH), E1, E2, SIGMA, &
      XIJ(3,MAXNH), KPOINT(3,NK), WK(NK)

  ! Dynamic arrays -------------------------------------------------------
  real(dp), dimension(:,:)  , pointer :: DTOT
  real(dp), dimension(:,:,:), pointer :: DPR

  ! Internal variables ---------------------------------------------------
  integer :: nuo, nhs, npsi, iuo, ihist, ispin, iunit1, iunit2, i

  integer :: iat, spec, ii, iorb, nspin_pdos

  logical :: orig_ParallelOverK, orig_Serial

  real(dp), dimension(:), pointer :: tmp
  real(dp), dimension(:), pointer :: eig => null()

  character*40 pos_string

  character(len=label_length+5) :: fnamepro
  character(len=label_length+4) :: fnametot

  real(dp) :: delta, ener

  external :: io_assign, io_close, timer, die

  type(xmlf_t) :: xf            ! For new XML output

#ifdef DEBUG
  call write_debug( '    PRE pdos' )
#endif

  call timer( 'pdos', 1)

  orig_Serial = Serial
  orig_ParallelOverK = ParallelOverK

  ! Find the intervals between the subdivisions in the energy scale ------
  delta = (E2 - E1) / (NHIST-1)

  ! Reset for Gamma-only PDOS calculations
  if ( Gamma ) ParallelOverK = .false.

  if ( Nodes > 1 .and. .not. ParallelOverK ) then
    Serial = .false.
  end if

  ! Find number of orbitals per unit cell 
#ifdef MPI
  call GetNodeOrbs(no_u,Node,Nodes,nuo)
#else
  nuo = no_u
#endif

  ! Check internal dimensions --------------------------------------------
  if ( nspin <= 2 .and. gamma) then
    nhs  = no_u * nuo
    npsi = no_u * no_l * nspin
  elseif ( nspin <= 2 .and. .not.gamma) then
#ifdef MPI
    if (ParallelOverK) then
      nhs  = 2 * no_u * no_u
      npsi = 2 * no_u * no_u
    else
#endif
      nhs  = 2 * no_u * nuo
      npsi = 2 * no_u * nuo
#ifdef MPI
    endif
#endif
  elseif (nspin >= 4) then
#ifdef MPI
    if(ParallelOverK) then
      nhs  = 2 * (2*no_u) * (2*no_u)
      npsi = 2 * (2*no_u) * (2*no_u)
    else
#endif
      nhs  = 2 * (2*no_u) * (2*nuo)
      npsi = 2 * (2*no_u) * (2*nuo)
#ifdef MPI
    endif
#endif
  else
    call die('diagon: ERROR: incorrect value of nspin')
  endif


  ! Allocate local arrays -------------------------------------------
  ! -----
  call allocDenseMatrix(nhs, nhs, npsi)
  nullify( dtot, dpr )
  if ( spin%H <= 2 ) then
    nspin_pdos = spin%H
    call re_alloc( dtot, 1, nhist, 1, nspin_pdos, 'dtot', 'pdos' )
    call re_alloc( dpr, 1, no_u, 1, nhist, 1, nspin_pdos, 'dpr', 'pdos' )
    call re_alloc( eig, 1, no_u, 'eig', 'pdos' )
  else ! must be larger than 4
    nspin_pdos = 4
    call re_alloc( dtot, 1, nspin_pdos, 1, nhist, 'dtot', 'pdos' )
    call re_alloc( dpr, 1, nspin_pdos, 1, no_u, 1, nhist, 'dpr', 'pdos' )
    call re_alloc( eig, 1, no_u*2, 'eig', 'pdos' )
  end if

  ! Initialize the projected density of states ---------------------------
  dtot(:,:) = 0._dp
  dpr(:,:,:) = 0._dp
  ! On return they will only be meaning full on Node == 0

  call print_diag()
  !  Call appropiate routine ----------------------------------------------
  if (nspin <= 2 .and. gamma) then
    call pdosg( nspin, NUO, NO, MAXNH, &
        no_u, NUMH, LISTHPTR, LISTH, H, S, &
        E1, E2, NHIST, SIGMA, INDXUO, &
        HAUX, SAUX, PSI, eig, DTOT, DPR )
  elseif ( nspin <= 2 .and. .not.gamma) then
#ifdef MPI
    if (ParallelOverK) then
      call pdoskp( nspin, NUO, NO, MAXNH, &
          no_u, NUMH, LISTHPTR, LISTH, H, S, &
          E1, E2, NHIST, SIGMA, &
          XIJ, INDXUO, NK, KPOINT, WK, &
          HAUX, SAUX, PSI, eig, DTOT, DPR )
    else
#endif
      call pdosk( nspin, NUO, NO, MAXNH, &
          no_u, NUMH, LISTHPTR, LISTH, H, S, &
          E1, E2, NHIST, SIGMA, &
          XIJ, INDXUO, NK, KPOINT, WK, &
          HAUX, SAUX, PSI, eig, DTOT, DPR )
#ifdef MPI
    endif
#endif    
  elseif (nspin == 4 .and. gamma) then
    call pdos2g( NUO, NO, NO_L, MAXNH, &
        no_u, NUMH, LISTHPTR, LISTH, H, S, &
        E1, E2, NHIST, SIGMA, INDXUO, &
        HAUX, SAUX, PSI, eig, DTOT, DPR )
  elseif (nspin == 4 .and. .not. gamma) then
#ifdef MPI
    if (ParallelOverK) then
      call pdos2kp( NUO, NO, NO_L, MAXNH, &
          no_u, NUMH, LISTHPTR, LISTH, H, S, &
          E1, E2, NHIST, SIGMA, &
          XIJ, INDXUO, NK, KPOINT, WK, &
          HAUX, SAUX, PSI, eig, DTOT, DPR )
    else
#endif      
      call pdos2k( NUO, NO, NO_L, MAXNH, &
          no_u, NUMH, LISTHPTR, LISTH, H, S, &
          E1, E2, NHIST, SIGMA, &
          XIJ, INDXUO, NK, KPOINT, WK, &
          HAUX, SAUX, PSI, eig, DTOT, DPR )
#ifdef MPI
    endif
#endif    
  elseif (nspin == 8 .and. gamma) then
    call pdos3g( NUO, NO, NO_L, MAXNH, &
        no_u, NUMH, LISTHPTR, LISTH, H, S, &
        E1, E2, NHIST, SIGMA, INDXUO, &
        HAUX, SAUX, PSI, eig, DTOT, DPR )
  elseif (nspin == 8 .and. .not. gamma) then
#ifdef MPI
    if (ParallelOverK) then
      call pdos3kp( NUO, NO, NO_L, MAXNH, &
          no_u, NUMH, LISTHPTR, LISTH, H, S, &
          E1, E2, NHIST, SIGMA, &
          XIJ, INDXUO, NK, KPOINT, WK, &
          HAUX, SAUX, PSI, eig, DTOT, DPR )
    else
#endif      
      call pdos3k( NUO, NO, NO_L, MAXNH, &
          no_u, NUMH, LISTHPTR, LISTH, H, S, &
          E1, E2, NHIST, SIGMA, &
          XIJ, INDXUO, NK, KPOINT, WK, &
          HAUX, SAUX, PSI, eig, DTOT, DPR )
#ifdef MPI
    endif
#endif    
  end if

  ! Clean up memory
  call de_alloc( eig, 'eig', 'pdos' )


  if (IOnode) then
    ! Open file for write on I/O node
    fnametot = trim(slabel)//'.DOS'
    call io_assign(iunit1)
    open(unit=iunit1, file=fnametot, form='formatted', status='unknown') 

    ! Output histogram
    select case ( spin%H )
    case ( 1 )
      do ihist = 1,nhist
        ENER = E1 + (ihist-1) * delta
        write(iunit1,'(2f20.5)') ener/ev,dtot(ihist,1)*eV
      enddo
    case ( 2 )
      do ihist = 1,nhist
        ENER = E1 + (ihist-1) * delta
        write(iunit1,'(3f20.5)') ener/ev,dtot(ihist,1)*eV, dtot(ihist,2)*eV
      enddo
    case default
      ! non-collinear case
      do ihist = 1,nhist
        ENER = E1 + (ihist-1) * delta
        write(iunit1,'(5f20.5)') ener/ev,dtot(1,ihist)*eV, &
            dtot(2,ihist)*eV,dtot(3,ihist)*eV,dtot(4,ihist)*eV
      enddo
    end select

    ! Close file for write
    call io_close(iunit1)
  endif

  ! New writing
  if (IOnode) then
    call xml_OpenFile(trim(slabel)//".PDOS.xml",xf,indent=.false.)

    fnamepro = trim(slabel)//'.PDOS'
    call io_assign(iunit2)
    open(iunit2,file=fnamepro,form='formatted',status='unknown')
    call xml_NewElement(xf,"pdos")
    call xml_NewElement(xf,"nspin")
    call xml_AddPCData(xf,nspin_pdos)
    call xml_EndElement(xf,"nspin")
    call xml_NewElement(xf,"norbitals")
    call xml_AddPCData(xf,no_u)
    call xml_EndElement(xf,"norbitals")

    write(iunit2,'(a)') '<pdos>'
    write(iunit2,'(a,i1,a)') '<nspin>', nspin_pdos, '</nspin>'
    write(iunit2,'(a,i0,a)') '<norbitals>', no_u, '</norbitals>'

    ! Write fermi-level
    call xml_NewElement(xf,"fermi_energy")
    call xml_AddAttribute(xf,"units","eV")
    call xml_AddPCData(xf, Ef / eV)
    call xml_EndElement(xf,"fermi_energy")
    write(iunit2,'(a,e17.9,a)') '<fermi_energy units="eV">', Ef / eV, '</fermi_energy>'

    ! Write sampled energies
    write(iunit2,'(a)') '<energy_values units="eV">'
    call xml_NewElement(xf,"energy_values")
    call xml_AddAttribute(xf,"units","eV")

    nullify( tmp )
    call re_alloc( tmp, 1, nspin_pdos*nhist, 'tmp', 'pdos' )
    do ihist = 1,nhist
      ENER = E1 + (ihist-1) * delta
      tmp(ihist) = ener/eV
      write(iunit2,'(e17.9)') tmp(ihist)
    enddo
    write(iunit2,'(a)') '</energy_values>'

    call xml_AddArray(xf,tmp(1:nhist))
    call xml_EndElement(xf,"energy_values")

    do i = 1, no_u
      iat = iaorb(i)
      iorb = iphorb(i)
      spec = isa(iat)

      call xml_NewElement(xf,"orbital")
      write(iunit2,'(a)') '<orbital '
      call xml_AddAttribute(xf,"index",i)
      call xml_dump_attribute(iunit2,"index",str(i))
      call xml_AddAttribute(xf,"atom_index",iat)
      call xml_dump_attribute(iunit2,"atom_index",str(iat))
      call xml_AddAttribute(xf,"species",trim(labelfis(spec)))
      call xml_dump_attribute(iunit2,"species",labelfis(spec))
      call xml_AddAttribute(xf,"Z",izofis(spec))
      call xml_dump_attribute(iunit2,"Z",str(izofis(spec)))
      write(pos_string,'(3f11.6)') xa(1:3,iat)
      call xml_AddAttribute(xf,"position",trim(pos_string))
      call xml_dump_attribute(iunit2,"position",pos_string)
      call xml_AddAttribute(xf,"n",cnfigfio(spec,iorb))
      call xml_dump_attribute(iunit2,"n",str(cnfigfio(spec,iorb)))
      call xml_AddAttribute(xf,"l",lofio(spec,iorb))
      call xml_dump_attribute(iunit2,"l",str(lofio(spec,iorb)))
      call xml_AddAttribute(xf,"m",mofio(spec,iorb))
      call xml_dump_attribute(iunit2,"m",str(mofio(spec,iorb)))
      call xml_AddAttribute(xf,"z",zetafio(spec,iorb))
      call xml_dump_attribute(iunit2,"z",str(zetafio(spec,iorb)))
      call xml_AddAttribute(xf,"P",pol(spec,iorb))
      call xml_dump_attribute(iunit2,"P",str(pol(spec,iorb)))
      write(iunit2,'(a)') '>'

      !------------------------------------------------------------
      call xml_NewElement(xf,"data")
      write(iunit2,'(a)') '<data>'
      select case ( spin%H )
      case ( 1 )
        do ihist = 1, nhist
          tmp(ihist) = dpr(i,ihist,1) * eV
          write(iunit2,'(e17.9)') tmp(ihist)
        end do
        call xml_AddArray(xf,tmp(1:nhist))
      case ( 2 )
        do ihist = 1, nhist
          tmp(ihist*2-1) = dpr(i,ihist,1) * eV
          tmp(ihist*2) = dpr(i,ihist,2) * eV
          write(iunit2,'(2(tr1,e17.9))') tmp(ihist*2-1:ihist*2)
        end do
        call xml_AddArray(xf,tmp(1:2*nhist))
      case ( 4, 8 )
        do ihist = 1, nhist
          tmp(ihist*4-3) = dpr(1,i,ihist) * eV
          tmp(ihist*4-2) = dpr(2,i,ihist) * eV
          tmp(ihist*4-1) = dpr(3,i,ihist) * eV
          tmp(ihist*4  ) = dpr(4,i,ihist) * eV
          write(iunit2,'(4f10.5)') tmp(ihist*4-3:ihist*4)
        end do
        call xml_AddArray(xf,tmp(1:4*nhist))
      end select
      call xml_EndElement(xf,"data")
      write(iunit2,'(a)') '</data>'
      call xml_EndElement(xf,"orbital")
      write(iunit2,'(a)') '</orbital>'

    end do

    ! Close file
    call xml_EndElement(xf,"pdos")
    call xml_Close(xf)
    write(iunit2,'(a)') '</pdos>'
    call io_close(iunit2)

    ! Free local workspace array
    call de_alloc( tmp, 'tmp', 'pdos' )

  endif

  ! Free local arrays ----------------------------------------------------
  call de_alloc( dpr, 'dpr', 'pdos' )
  call de_alloc( dtot, 'dtot', 'pdos' )
  call resetDenseMatrix()

  ! Restore the paralleloverk options
  ParallelOverK = orig_ParallelOverK
  Serial = orig_Serial

  call timer( 'pdos', 2)

#ifdef DEBUG
  call write_debug( '    POS pdos' )
#endif

end subroutine pdos
