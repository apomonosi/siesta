subroutine print_spin(qspin)
  !
  ! Prints spin in output and CML files
  !
  use m_spin,          only: spin

  use siesta_cml
  use parallel,        only: IOnode
  use precision,       only: dp

  implicit none

  real(dp), intent(in)  :: qspin(spin%Grid)

  real(dp) :: qaux(spin%Grid)
  real(dp) :: Stot        ! Total spin magnitude
  real(dp) :: Svec(3)     ! Total spin vector

  ! We are left with printing out to stdout
  if ( .not. IONode ) return

  if ( spin%Grid == 2) then

    Svec(1) = 0.0_dp
    Svec(2) = 0.0_dp
    Svec(3) = qspin(1) - qspin(2)
    Stot = Svec(3)
    write(6,'(5x,a,2(f10.1,tr1),f10.5,tr1,"}",tr2,f10.5)') &
        'spin moment: {S} , |S| = { ', Svec, Stot
    if (cml_p) call cmlAddProperty(xf=mainXML,            &
        value=qspin(1)-qspin(2), dictref='siesta:stot', &
        units='siestaUnits:spin')

  else if ( spin%Grid == 4 ) then

    call spnvec( spin%Grid, qspin, qaux, Stot, Svec )
    write(6,'(5x,a,2(f10.5,tr1),f10.5,tr1,"}",tr2,f10.5)') &
        'spin moment: {S} , |S| = { ', Svec, Stot
    if (cml_p) then
      call cmlAddProperty(xf=mainXML, value=Stot,  &
          dictref='siesta:stot', units='siestaUnits:spin')
      call cmlAddProperty(xf=mainXML, value=Svec,  &
          dictref='siesta:svec', units='siestaUnits:spin')
    end if !cml_p

  end if

end subroutine print_spin
