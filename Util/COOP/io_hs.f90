! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module io_hs

implicit none
public :: read_hs_file

CONTAINS

  function HSX_version(fname) result(version)
    character(len=*), intent(in) :: fname
    integer :: version
    integer :: iu
    integer :: na_u, no_u, no_s, nspin, n_nzs, err

    external :: io_assign, io_close

    ! Initialize
    version = 0

    ! Open file
    call io_assign( iu )
    open( iu, file=fname, form='unformatted', status='unknown' )

    read(iu,iostat=err) na_u, no_s, nspin, n_nzs
    if ( err == 0 ) then
       ! we can successfully read 4 integers
       version = 0
    else
       rewind(iu)
       read(iu,iostat=err) version
    end if

    call io_close(iu)

  end function HSX_version


  subroutine read_hs_file(fname)
    use main_vars

    character(len=*), intent(in) :: fname

    integer :: version

    version = HSX_version(fname)
    if ( version == 0 ) then
      call read_hs_file_old(fname)
    else if ( (version == 1) .or. (version == 2) ) then
      call read_hs_file_version1_2(fname)
    else
      STOP "unknown HSX file version [0, 2]"
    end if
  end subroutine read_hs_file

  subroutine read_hs_file_version1_2(fname)
    use main_vars
    use precision, only: sp

    character(len=*), intent(in) :: fname

    integer :: numx, ind
    real(dp) :: fucell(3,3), kdispl(3)
    logical :: is_dp
    integer :: fnsc(3), kcell(3,3)
    integer :: ja, jo, version
    integer, allocatable :: isc_off(:,:), lasto(:)
    real(dp), allocatable :: xa(:,:)
    real(sp), allocatable :: rbuf(:)

    write(6,"(1x,a)",advance='no') trim(fname)
    open(hs_u,file=trim(fname),status='old',form='unformatted')

    ! skip version
    read(hs_u, iostat=iostat) version
    if (iostat /= 0) STOP "version"

    if ( version > 2 ) then
      STOP "in correct call [version]"
    end if

    ! now read whether this is double precision or not
    read(hs_u, iostat=iostat) is_dp
    if ( iostat /= 0 ) STOP "is_dp"

    read(hs_u,iostat=iostat) na_u, nnao, nspin, nspecies, fnsc
    if ( iostat /= 0 ) STOP "geometry dimensions"
    if (nnao /= nao) STOP "norbs inconsistency"
    no_u = nnao
    no_s = nnao * product(fnsc)

    read(hs_u,iostat=iostat) fucell, Efermi, qtot, temp_in_file

    ! Allocate the arrays for the atomic information
    ! First local arrays that won't be shared
    allocate(xa(3,na_u), lasto(0:na_u), isc_off(3, product(fnsc)))
    lasto(0) = 0
    ! allocate shared variables
    allocate(isa(na_u), iaorb(no_u), iphorb(no_u))
    allocate(label(nspecies), zval(nspecies), no(nspecies))

    ! Now read data
    read(hs_u,iostat=iostat) isc_off, xa, isa, lasto(1:)
    if ( iostat /= 0 ) STOP "geometry coordinates"
    read(hs_u,iostat=iostat) (label(is), zval(is), no(is), is=1,nspecies)
    if ( iostat /= 0 ) STOP "species information"

    ! Allocate remaning orbital information
    naoatx = maxval(no)
    allocate(nquant(nspecies,naoatx), lquant(nspecies,naoatx))
    allocate(zeta(nspecies,naoatx))
    do is = 1, nspecies
      read(hs_u,iostat=iostat) (nquant(is,io), lquant(is,io), zeta(is,io), io=1,no(is))
      if ( iostat /= 0 ) STOP "specific specie information"
    end do

    if ( version == 2 ) then
       read(hs_u,iostat=iostat) kcell, kdispl
       if ( iostat /= 0 ) STOP "k-information"
    end if

    ! Populate iaorb and iphorb
    ind = 0
    do ia = 1, na_u
      do io = 1, lasto(ia) - lasto(ia-1)
        ind = ind + 1
        iaorb(ind) = ia
        iphorb(ind) = io
      end do
    end do

    ! Now create indxuo (siesta always produces the same order)
    allocate(indxuo(no_u * product(fnsc)))
    ind = 0
    do is = 1 , product(fnsc)
      do i = 1, no_u
        ind = ind + 1
        indxuo(ind) = i
      end do
    end do

    if (wfs_x.and.(nsp /= nspin)) STOP " nspin not the same for WFS and HSX"
    nsp = nspin
    h_spin_dim = nspin
    allocate(numh(no_u))
    read(hs_u,iostat=iostat) numh ! numh
    if (iostat /= 0) STOP "numh(io)"

    ! Create pointer
    allocate(listhptr(no_u))
    listhptr(1) = 0
    do io=2,no_u
      listhptr(io) = listhptr(io-1) + numh(io-1)
    end do

    ! Now read sparse data
    nh = listhptr(no_u) + numh(no_u)
    allocate(listh(nh))

    do io=1,no_u
      read(hs_u,iostat=iostat) listh(listhptr(io)+1:listhptr(io)+numh(io))
      if (iostat /= 0) STOP "listh"
    end do


    ! Now we need to re-create the xij array
    allocate(xij(3,nh),dij(nh))

    ! Create xij and dij
    do ia = 1 , na_u
      do io = lasto(ia-1) + 1, lasto(ia)
        do ind = listhptr(io) + 1, listhptr(io) + numh(io)
          is = (listh(ind) - 1) / no_u + 1
          ja = iaorb(ucorb(listh(ind), no_u))

          xij(:, ind) = xa(:,ja) - xa(:,ia) + &
              isc_off(1,is) * fucell(:,1) + &
              isc_off(2,is) * fucell(:,2) + &
              isc_off(3,is) * fucell(:,3)

          dij(ind) = sqrt(dot_product(xij(:,ind),xij(:,ind))) / Ang
        end do
      end do
    end do

    ! Clean-up
    deallocate(xa, isc_off, lasto)

    ! Read H and S
    allocate(hamilt(nh,nspin))
    allocate(Sover(nh))

    if ( is_dp ) then
      do is = 1, nspin
        do io = 1, no_u
          read(hs_u,iostat=iostat) hamilt(listhptr(io)+1:listhptr(io)+numh(io),is)
          if (iostat /= 0) STOP "H(dp)"
        end do
      end do

      do io = 1, no_u
        read(hs_u,iostat=iostat) Sover(listhptr(io)+1:listhptr(io)+numh(io))
        if (iostat /= 0) STOP "S(dp)"
      end do

    else
      allocate(rbuf(maxval(numh)))

      do is = 1, nspin
        do io = 1, no_u
          read(hs_u,iostat=iostat) rbuf(1:numh(io))
          if (iostat /= 0) STOP "H(sp)"
          hamilt(listhptr(io)+1:listhptr(io)+numh(io),is) = rbuf(1:numh(io))
        end do
      end do

      do io = 1, no_u
        read(hs_u,iostat=iostat) rbuf(1:numh(io))
        if (iostat /= 0) STOP "S(sp)"
        Sover(listhptr(io)+1:listhptr(io)+numh(io)) = rbuf(1:numh(io))
      end do

      deallocate(rbuf)
    end if

    close(hs_u)

  contains

    elemental function UCORB(a,p)
      integer, intent(in) :: a,p
      integer :: UCORB
      UCORB = MOD(a-1,p) + 1
    end function

  end subroutine read_hs_file_version1_2


  subroutine read_hs_file_old(fname)
    use main_vars
    use precision, only: sp

    character(len=*), intent(in) :: fname

    integer, allocatable  :: ibuff(:)
    real(sp), allocatable  :: hbuff(:)
    real(sp), allocatable  :: buff3(:,:)

    integer numx, ind
    logical lacking_indxuo

    write(6,"(1x,a)",advance='no') trim(fname)
    open(hs_u,file=trim(fname),status='old',form='unformatted')

    read(hs_u,iostat=iostat) nnao, no_s, nspin, nh
    print *, "nnao, no_s, nspin, nh:",  nnao, no_s, nspin, nh
    if (iostat /= 0) STOP "nnao, no_s..."
    if (nnao /= nao) STOP "norbs inconsistency"
    no_u = nao

    ! In modern versions of HSX files this should be always .false., that is,
    ! files always include the indxuo array, even if it is trivial.

    read(hs_u,iostat=iostat) lacking_indxuo
    if (iostat /= 0) STOP "lacking_indxuo"
    IF (DEBUG) PRINT *, "LACKING_INDXUO=", lacking_indxuo

    if (.not. lacking_indxuo) then
      ! read it
      allocate(indxuo(no_s))
      read(hs_u) (indxuo(i),i=1,no_s)
    else
      allocate(indxuo(no_u))
      ! build it
      do i=1,no_u
        indxuo(i) = i
      enddo
    endif

    if (debug) print *, "HS read: nh, nsp, nnao: ", nh, nspin, nnao
    if (nnao.ne.nao) STOP " nnao .ne. nao in HS"

    if (wfs_x.and.(nspin.ne.nsp)) STOP " nspin .ne. nsp in HS"
    nsp=nspin
    h_spin_dim = nspin
    allocate (numh(nao), listhptr(nao), listh(nh))

    allocate (hamilt(nh,nspin))
    allocate (Sover(nh))
    allocate (xij(3,nh),dij(nh))

    read(hs_u,iostat=iostat) (numh(io), io=1,no_u)         ! numhg
    if (iostat /= 0) STOP "numh(io)"
    do io=1,no_u
      if (debug) print *, "numhg ", io, numh(io)
    enddo

    numx = maxval(numh(:))
    allocate(ibuff(numx), hbuff(numx), buff3(3,numx))


    ! Create listhptr
    listhptr(1)=0
    do io=2,no_u
      listhptr(io)=listhptr(io-1)+numh(io-1)
    enddo
    if (listhptr(no_u)+numh(no_u).gt.nh) STOP "nh overflow in HS"

    do io=1,no_u
      read(hs_u,iostat=iostat) (ibuff(im), im=1,numh(io))
      if (iostat /= 0) STOP "listh"
      do im=1,numh(io)
        listh(listhptr(io)+im) = ibuff(im)
        if (debug) print *, "listh ", io, im, listh(listhptr(io)+im)
      enddo
    enddo

    do is=1,nspin
      do io=1,no_u
        read(hs_u,iostat=iostat) (hbuff(im), im=1,numh(io))
        if (iostat /= 0) STOP "Hamilt"
        do im=1,numh(io)
          hamilt(listhptr(io)+im,is) = hbuff(im)
          if (debug) print *, "Hamilt ", io, im, hbuff(im)
        enddo
      enddo
    enddo
    !
    !       Read overlap matrix
    !
    do io=1,no_u
      read(hs_u,iostat=iostat) (hbuff(im), im=1,numh(io))
      if (iostat /= 0) STOP "Overlap matrix read error"
      do im=1,numh(io)
        Sover(listhptr(io)+im) = hbuff(im)
        if (debug) print *, "S ", io, im, hbuff(im)
      enddo
    enddo

    read(hs_u,iostat=iostat) qtot, temp_in_file
    if (debug) print *, "QTOT, Temp in file: ", qtot, temp_in_file
    if (iostat /= 0) then
      if (debug) print *, "iostat:", iostat
      STOP "qtot, temp in file"
    endif

    !
    !        Always read xijk
    !
    do io=1,no_u
      read(hs_u,iostat=iostat) ((buff3(k,im), k=1,3), im=1,numh(io))
      if (iostat /= 0) STOP "xij(k)"
      do im=1,numh(io)
        ind = listhptr(io)+im
        if (debug) print *, "xijk ", buff3(:,im)
        xij(1:3,ind) = buff3(1:3,im)
        dij(ind) = sqrt(dot_product(buff3(:,im),buff3(:,im))) / Ang
      enddo
    enddo
    !
    !        Read auxiliary info
    !
    read(hs_u) nspecies
    allocate(label(nspecies), zval(nspecies), no(nspecies))
    read(hs_u) (label(is),zval(is),no(is), is=1,nspecies)
    naoatx = maxval(no(1:nspecies))
    allocate (nquant(nspecies,naoatx), lquant(nspecies,naoatx), &
        zeta(nspecies,naoatx))
    do is=1, nspecies
      do io=1, no(is)
        read(hs_u) nquant(is,io), lquant(is,io), zeta(is,io)
      enddo
    enddo
    read(hs_u) na_u
    allocate(isa(na_u))
    allocate(iaorb(no_u), iphorb(no_u))
    read(hs_u) (isa(ia), ia=1,na_u)
    read(hs_u) (iaorb(io), iphorb(io), io=1,no_u)

    close(hs_u)
    deallocate(ibuff, hbuff, buff3)

  end subroutine read_hs_file_old

end module io_hs
