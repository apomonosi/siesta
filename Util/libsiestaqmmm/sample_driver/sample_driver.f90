! ---
! Copyright (C) 1996-2021       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This program is a sample driver that relies on the libsiestaqmmm library
! to call SIESTA as a QM workhorse for QMMM.
!
! It should be just compiled and linked with the library, and then run. Note
! that FDF is also a dependency.
!
! For example:
!    gfortran -o sampledriver sample_driver.f90 \
!             -L/PATH/TO/SIESTA/LIBS/ -lsiestaqmmm -lfdf
!
! See the siestaqmmm.fdf file in this same folder for a sample input file
! for this driver.

program sampledriver
    implicit none
    integer      :: nat = 3, npc = 3
    logical      :: is_last
    real(kind=8) :: ucell(3,3), qm_crd(3,3), mm_crd(3,3), mm_crg(3), &
                    Etot, qm_frc(3,3), mm_frc(3,3), stress(3,3)

    integer :: istep

    real(kind=8), parameter :: ANG = 1.0d0 / 0.529177d0
    real(kind=8), parameter :: EV  = 1.0d0 / 13.60580d0

    external :: get_siesta_forces_and_stress

    ucell(1:3,1) = (/20.0d0,  0.0d0,  0.0d0/)
    ucell(1:3,2) = (/ 0.0d0, 20.0d0,  0.0d0/)
    ucell(1:3,3) = (/ 0.0d0,  0.0d0, 20.0d0/)

    qm_crd(1:3,1) = (/1.000d0, 0.000d0, 0.0d0/)
    qm_crd(1:3,2) = (/1.757d0, 0.586d0, 0.0d0/)
    qm_crd(1:3,3) = (/0.243d0, 0.586d0, 0.0d0/)

    mm_crd(1:3,1) = (/1.000d0, 2.500d0, 0.0d0/)
    mm_crd(1:3,2) = (/1.757d0, 3.086d0, 0.0d0/)
    mm_crd(1:3,3) = (/0.243d0, 3.086d0, 0.0d0/)

    mm_crg(:) = (/-0.834d0, 0.417d0, 0.417d0/)

    write( *, '(A)') "Calling SIESTA as QMMM subroutine for a single step."
    is_last = .true.
    call get_siesta_forces_and_stress( nat, npc, qm_crd, mm_crd, mm_crg,    &
                                       ucell, is_last, Etot, qm_frc, mm_frc,&
                                       stress )

    open( unit = 36155, file = "driver_results.out" )
    write( 36155,'(A,F17.6)' ) "Total Energy (eV) = ", Etot / EV
    write( 36155,'(/,A,3(/,A,3F12.6))' ) "Stress tensor (eV/Ang**3):", &
      "      ", stress(:,1) * ANG * ANG * ANG / EV, &
      "      ", stress(:,2) * ANG * ANG * ANG / EV, &
      "      ", stress(:,3) * ANG * ANG * ANG / EV

    write( 36155, * )
    write( 36155,'(A,I6)' ) "QM Atoms = ", nat
    write( 36155,'(/,A,3(/,A,3F12.6))' ) "Forces on QM (eV/Ang):", &
      "      ", qm_frc(:,1) * ANG / EV, &
      "      ", qm_frc(:,2) * ANG / EV, &
      "      ", qm_frc(:,3) * ANG / EV

    write( 36155, * )
    write( 36155,'(A,I6)' ) "MM Atoms = ", npc
    write( 36155,'(/,A,3(/,A,3F12.6))' ) "Forces on MM (eV/Ang):", &
      "      ", mm_frc(:,1) * ANG / EV, &
      "      ", mm_frc(:,2) * ANG / EV, &
      "      ", mm_frc(:,3) * ANG / EV

    close( 36155 )

    write( *, '(A)') "Calling SIESTA as QMMM subroutine for a five steps."

    is_last = .false.
    do istep = 1, 5
      if ( istep == 5 ) is_last = .true.
      call get_siesta_forces_and_stress( nat, npc, qm_crd, mm_crd, mm_crg,    &
                                         ucell, is_last, Etot, qm_frc, mm_frc,&
                                         stress )

      write( *, '(A12,I2,A12,F16.6,A3)') &
        "==> Step no.", istep, " - Energy = ", Etot / EV, " eV"
      qm_crd(:,:) = qm_crd(:,:) + 0.1d0 * qm_frc(:,:)
      mm_crd(:,:) = mm_crd(:,:) + 0.1d0 * mm_frc(:,:)
      ucell(:,:)  = ucell(:,:)  + 0.1d0 * stress(:,:)
    enddo
end program sampledriver
