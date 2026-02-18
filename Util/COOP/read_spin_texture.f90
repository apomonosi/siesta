! This file has been contributed by Roberto Robles,
! and modified by Alberto Garcia to adapt it to the
! (shortcomings) of the current output format of the
! 'spin_texture' program.
!
!************************************************************
! WARNING: In this version, the number of k-points, energies,
! Efermi, and bands to process, are entered through standard
! input.
!************************************************************

! It is included to provide some perspective on ways to
! plot the spin-texture information.
! It might need some re-formatting of the default spin-texture
! data provided by the program 'spin_texture' in this directory.
! An alternative to the procedure indicated is to select a single
! band with 'spin_texture', and do some semi-automatic processing
! of the resulting file for input to gnuplot or similar tools.
!
! One needs to plot a vector field where each point k is assigned a
! vector (Sx, Sy, Sz). To achieve this satisfactorily is not
! completely trivial, and different techniques can be used.
!
!  Program to extract the spin texture from file spin_texture.dat as
!  printed by SIESTA. It requires the number of bands and their index
!  provided interactively.
!
!  If executed as "read_spin_texture_xsf" it gives the output as the
!  head of a xsf file, elsewhere it produces a .xyz file.
!
      program read_spin_texture

      double precision, allocatable :: kpoints(:,:), st(:,:,:), en(:,:)

      integer, allocatable :: bands(:)

      integer :: nkp, nen, nb, band, ext_ind, iarg, narg
      double precision:: ef

      logical :: xyz

      character(len=100) :: ofile, ab, exe_name
      character(len=10)  :: cmd_option

      xyz  = .false.
      narg = command_argument_count()

      if ( narg < 6 ) then
        write(*,*) "Wrong number of arguments."
        call print_help( )
        stop
      endif

      iarg = 0
      do while ( iarg < narg )

        iarg = iarg + 1
        call get_command_argument( iarg, cmd_option )
        if ( cmd_option(1:1) /= '-' ) cycle

        select case ( trim(cmd_option) )
        case ( "-h" )
          call print_help( )
          stop

        case ( "-nk" )
          if ( iarg >= narg ) then
            write(*,*) "Missing value for option: ", trim(cmd_option)
            stop
          endif

          iarg = iarg + 1
          call get_command_argument( iarg, cmd_option )
          read(cmd_option,*) nkp

        case ( "-nbands" )
          if ( iarg >= narg ) then
            write(*,*) "Missing value for option: ", trim(cmd_option)
            stop
          endif

          iarg = iarg + 1
          call get_command_argument( iarg, cmd_option )
          read(cmd_option,*) nen

        case ( "-efermi" )
          if ( iarg >= narg ) then
            write(*,*) "Missing value for option: ", trim(cmd_option)
            stop
          endif

          iarg = iarg + 1
          call get_command_argument( iarg, cmd_option )
          read(cmd_option,*) ef

        case ( "-xsf" )
          xyz = .false.

        case ( "-xyz" )
          xyz = .true.

        end select

      enddo

      if ( xyz ) then
        write(*,"(A)") "Output will be written in XYZ format."
      else
        write(*,"(A)") "Output will be written in XSF format."
      endif
      write(*,"(A36,I6,I6,f10.4)") &
        "k-points, bands per k-point, Efermi:", nkp, nen, ef

      allocate( kpoints(3,nkp) )
      allocate( st(3,nkp,nen) )
      allocate( en(nkp,nen) )
      allocate( bands(nen) )

      open(unit=9,file='spin_texture.dat',status='old')

      !     Skip two lines
      read(9,*)

      do ik = 1, nkp
        read(9,*) ! Skip empty line.
        read(9,'(14x,3f12.6)') (kpoints(j,ik),j=1,3)
        !write(*,'(i4,3f12.6)') ik,(kpoints(j,ik),j=1,3)
        read(9,*) ! Skip "ie     e(eV)      Sx      Sy      Sz" line.
        do ie = 1, nen
          read(9,'(7x,f12.5,3f8.4)') en(ik,ie),(st(j,ik,ie),j=1,3)
          !write(*,'(i4,f12.5,3f8.4)') ie,en(ik,ie),(st(j,ik,ie),j=1,3)
        enddo
      enddo

      write(*,'(a)') "Number of bands: "
      read(*,'(i4)') nb
      do ib=1,nb
        write(*,'(a,i4)') "Index of band # ", ib
        read(*,'(i4)') bands(ib)
        if (bands(ib) .gt. nen) then
          write(*,*) "Band index can not be bigger than ", nen
          stop
        endif

        write(ab,'(i0)') bands(ib)
        ofile = "st_band_"//trim(ab)//".xsf"
        if (xyz) ofile = "st_band_"//trim(ab)//".xyz"

        ifile = 67

        open(ifile,file=ofile)

        if (xyz) write(ifile,'(i4/)') nkp
        if (.not.xyz) write(ifile,'(a/a)') "set xsfStructure {","ATOMS"

        do ik=1,nkp
          if(xyz) write(ifile,'(a,3f12.6,3f8.4)') "X ",&
               (kpoints(j,ik)*10,j=1,2), (en(ik,bands(ib))-ef), &
               (st(j,ik,bands(ib)),j=1,3)
          if(.not.xyz) write(ifile,'(i3,3f12.6,3f8.4)') 0,&
               (kpoints(j,ik)*10,j=1,2), (en(ik,bands(ib))-ef), &
               (st(j,ik,bands(ib)),j=1,3)
        enddo
        if (.not.xyz) write(ifile,'(a)') "}"
        close(ifile)
     enddo


     contains

    subroutine print_help()
      implicit none

      write(*,*)
      write(*,*) " read_spin_texture usage: "
      write(*,*) " "
      write(*,*) " read_spin_texture -nk [] -bands [] -efermi [] (-xyz/-xsf)"
      write(*,*) " Mandatory inputs: "
      write(*,*) "    -nk     : Number of k-points."
      write(*,*) "    -nbands : Bands per k-point."
      write(*,*) "    -efermi : Energy of the fermi level."
      write(*,*) " "
      write(*,*) " Optional inputs: "
      write(*,*) "    -xsf     : Write output in xsf format (default)."
      write(*,*) "    -xyz     : Write output in xyz format."
    end subroutine print_help

end program read_spin_texture
