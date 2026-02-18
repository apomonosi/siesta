!
! Copyright (C) 1996-2022	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module units_legacy_m

  use units_common_m

  implicit none
  private

  integer, parameter :: dp = selected_real_kind(14,100)

  ! pi to 50 digits
  real(dp), parameter, public :: pi = 3.14159265358979323846264338327950288419716939937510_dp
  real(dp), parameter, public :: pi2 = pi * 2._dp
  real(dp), parameter, public :: deg = pi / 180.0_dp

  real(dp), parameter, public :: kg     = 1.28460971172331e27_dp
  real(dp), parameter, public :: Joule  = 1._dp/2.1798723611035e-18_dp
  real(dp), parameter, public :: h_planck = 6.62607015e-19_dp * Joule ! in Ry*fs
  real(dp), parameter, public :: hbar     = h_planck / pi2 ! in Ry*fs
  real(dp), parameter, public :: c_light  = 299792458.0_dp ! in m/s

  real(dp), parameter, public :: Ang    = 1._dp / 0.529177_dp
  real(dp), parameter, public :: eV     = 1._dp / 13.60580_dp
  real(dp), parameter, public :: kBar   = 1._dp / 1.47108e5_dp
  real(dp), parameter, public :: Kelvin = eV / 11604.45_dp
  real(dp), parameter, public :: Debye  = 0.393430_dp
  real(dp), parameter, public :: amu    = 2.133107_dp
  real(dp), parameter, public :: Ryd_time = 1._dp/0.04837769_dp

  real(dp), parameter, public :: GPa = kBar * 10

  public :: inquire_unit

  ! Data generated from table: legacy
  integer, parameter :: nu = 134
  character(8), save :: dimm(nu)
  character(10), save :: name(nu)
  real(dp), save :: unit(nu)
  data dimm(1)/'mass    '/, name(1)/'g         '/, unit(1)/1.d-3/
  data dimm(2)/'mass    '/, name(2)/'kg        '/, unit(2)/1.d0/
  data dimm(3)/'mass    '/, name(3)/'amu       '/, unit(3)/1.66054d-27/

  data dimm(4)/'length  '/, name(4)/'m         '/, unit(4)/1.d0/
  data dimm(5)/'length  '/, name(5)/'cm        '/, unit(5)/1.d-2/
  data dimm(6)/'length  '/, name(6)/'nm        '/, unit(6)/1.d-9/
  data dimm(7)/'length  '/, name(7)/'pm        '/, unit(7)/1.d-12/
  data dimm(8)/'length  '/, name(8)/'Ang       '/, unit(8)/1.d-10/
  data dimm(9)/'length  '/, name(9)/'Bohr      '/, unit(9)/0.529177d-10/

  data dimm(10)/'energy  '/, name(10)/'J         '/, unit(10)/1.d0/
  data dimm(11)/'energy  '/, name(11)/'kJ        '/, unit(11)/1.d3/
  data dimm(12)/'energy  '/, name(12)/'erg       '/, unit(12)/1.d-7/
  data dimm(13)/'energy  '/, name(13)/'meV       '/, unit(13)/1.60219d-22/
  data dimm(14)/'energy  '/, name(14)/'eV        '/, unit(14)/1.60219d-19/
  data dimm(15)/'energy  '/, name(15)/'mRy       '/, unit(15)/2.17991d-21/
  data dimm(16)/'energy  '/, name(16)/'Ry        '/, unit(16)/2.17991d-18/
  data dimm(17)/'energy  '/, name(17)/'mHa       '/, unit(17)/4.35982d-21/
  data dimm(18)/'energy  '/, name(18)/'mHartree  '/, unit(18)/4.35982d-21/
  data dimm(19)/'energy  '/, name(19)/'Ha        '/, unit(19)/4.35982d-18/
  data dimm(20)/'energy  '/, name(20)/'Hartree   '/, unit(20)/4.35982d-18/
  data dimm(21)/'energy  '/, name(21)/'K         '/, unit(21)/1.38066d-23/
  data dimm(22)/'energy  '/, name(22)/'Kelvin    '/, unit(22)/1.38066d-23/
  data dimm(23)/'energy  '/, name(23)/'kcal/mol  '/, unit(23)/6.94780d-21/
  data dimm(24)/'energy  '/, name(24)/'kJ/mol    '/, unit(24)/1.6606d-21/
  data dimm(25)/'energy  '/, name(25)/'Hz        '/, unit(25)/6.6262d-34/
  data dimm(26)/'energy  '/, name(26)/'THz       '/, unit(26)/6.6262d-22/
  data dimm(27)/'energy  '/, name(27)/'cm-1      '/, unit(27)/1.986d-23/
  data dimm(28)/'energy  '/, name(28)/'cm^-1     '/, unit(28)/1.986d-23/
  data dimm(29)/'energy  '/, name(29)/'cm**-1    '/, unit(29)/1.986d-23/

  data dimm(30)/'time    '/, name(30)/'s         '/, unit(30)/1.d0/
  data dimm(31)/'time    '/, name(31)/'ns        '/, unit(31)/1.d-9/
  data dimm(32)/'time    '/, name(32)/'ps        '/, unit(32)/1.d-12/
  data dimm(33)/'time    '/, name(33)/'fs        '/, unit(33)/1.d-15/
  data dimm(34)/'time    '/, name(34)/'min       '/, unit(34)/60.d0/
  data dimm(35)/'time    '/, name(35)/'mins      '/, unit(35)/60.d0/
  data dimm(36)/'time    '/, name(36)/'hour      '/, unit(36)/3600.d0/
  data dimm(37)/'time    '/, name(37)/'hours     '/, unit(37)/3600.d0/
  data dimm(38)/'time    '/, name(38)/'day       '/, unit(38)/86400.d0/
  data dimm(39)/'time    '/, name(39)/'days      '/, unit(39)/86400.d0/

  data dimm(40)/'force   '/, name(40)/'N         '/, unit(40)/1.d0/
  data dimm(41)/'force   '/, name(41)/'eV/Ang    '/, unit(41)/1.60219d-9/
  data dimm(42)/'force   '/, name(42)/'Ry/Bohr   '/, unit(42)/4.11943d-8/
  data dimm(43)/'force   '/, name(43)/'Ha/Bohr   '/, unit(43)/8.23886d-08/

  data dimm(44)/'pressure'/, name(44)/'Pa        '/, unit(44)/1.d0/
  data dimm(45)/'pressure'/, name(45)/'GPa       '/, unit(45)/1.d9/
  data dimm(46)/'pressure'/, name(46)/'atm       '/, unit(46)/1.01325d5/
  data dimm(47)/'pressure'/, name(47)/'bar       '/, unit(47)/1.d5/
  data dimm(48)/'pressure'/, name(48)/'kbar      '/, unit(48)/1.d8/
  data dimm(49)/'pressure'/, name(49)/'Mbar      '/, unit(49)/1.d11/
  data dimm(50)/'pressure'/, name(50)/'eV/Ang**3 '/, unit(50)/1.60219d11/
  data dimm(51)/'pressure'/, name(51)/'eV/Ang^3  '/, unit(51)/1.60219d11/
  data dimm(52)/'pressure'/, name(52)/'Ry/Bohr**3'/, unit(52)/1.47108d13/
  data dimm(53)/'pressure'/, name(53)/'Ry/Bohr^3 '/, unit(53)/1.47108d13/
  data dimm(54)/'pressure'/, name(54)/'Ha/Bohr^3 '/, unit(54)/2.94216d13/
  data dimm(55)/'pressure'/, name(55)/'Ha/Bohr**3'/, unit(55)/2.94216d13/

  data dimm(56)/'surftens'/, name(56)/'N/m       '/, unit(56)/1.d0/
  data dimm(57)/'surftens'/, name(57)/'mN/m      '/, unit(57)/1.d3/
  data dimm(58)/'surftens'/, name(58)/'dyn/cm    '/, unit(58)/1.d3/
  data dimm(59)/'surftens'/, name(59)/'erg/cm**2 '/, unit(59)/1.d3/

  data dimm(60)/'charge  '/, name(60)/'c         '/, unit(60)/1.d0/
  data dimm(61)/'charge  '/, name(61)/'e         '/, unit(61)/1.602177d-19/

  data dimm(62)/'dipole  '/, name(62)/'c*m       '/, unit(62)/1.d0/
  data dimm(63)/'dipole  '/, name(63)/'D         '/, unit(63)/3.33564d-30/
  data dimm(64)/'dipole  '/, name(64)/'Debye     '/, unit(64)/3.33564d-30/
  data dimm(65)/'dipole  '/, name(65)/'e*Bohr    '/, unit(65)/8.47835d-30/
  data dimm(66)/'dipole  '/, name(66)/'e*Ang     '/, unit(66)/1.602177d-29/

  data dimm(67)/'mominert'/, name(67)/'kg*m**2   '/, unit(67)/1.d0/
  data dimm(68)/'mominert'/, name(68)/'Ry*fs**2  '/, unit(68)/2.17991d-48/

  data dimm(69)/'efield  '/, name(69)/'V/m       '/, unit(69)/1.d0/
  data dimm(70)/'efield  '/, name(70)/'V/cm      '/, unit(70)/1.d2/
  data dimm(71)/'efield  '/, name(71)/'V/um      '/, unit(71)/1.d6/
  data dimm(72)/'efield  '/, name(72)/'V/nm      '/, unit(72)/1.d9/
  data dimm(73)/'efield  '/, name(73)/'V/Ang     '/, unit(73)/1.d10/
  data dimm(74)/'efield  '/, name(74)/'V/Bohr    '/, unit(74)/1.8897268d10/
  data dimm(75)/'efield  '/, name(75)/'Ry/Bohr/e '/, unit(75)/2.5711273d11/
  data dimm(76)/'efield  '/, name(76)/'Ha/Bohr/e '/, unit(76)/5.1422546d11/
  data dimm(77)/'efield  '/, name(77)/'Har/Bohr/e'/, unit(77)/5.1422546d11/

  data dimm(78)/'angle   '/, name(78)/'deg       '/, unit(78)/1.d0/
  data dimm(79)/'angle   '/, name(79)/'rad       '/, unit(79)/5.72957795d1/

  data dimm(80)/'torque  '/, name(80)/'meV/deg   '/, unit(80)/1.0d-3/
  data dimm(81)/'torque  '/, name(81)/'meV/rad   '/, unit(81)/1.745533d-5/
  data dimm(82)/'torque  '/, name(82)/'eV/deg    '/, unit(82)/1.0d0/
  data dimm(83)/'torque  '/, name(83)/'eV/rad    '/, unit(83)/1.745533d-2/
  data dimm(84)/'torque  '/, name(84)/'mRy/deg   '/, unit(84)/13.6058d-3/
  data dimm(85)/'torque  '/, name(85)/'mRy/rad   '/, unit(85)/0.237466d-3/
  data dimm(86)/'torque  '/, name(86)/'Ry/deg    '/, unit(86)/13.6058d0/
  data dimm(87)/'torque  '/, name(87)/'Ry/rad    '/, unit(87)/0.237466d0/

  data dimm(88)/'bfield  '/, name(88)/'Tesla     '/, unit(88)/1./
  data dimm(89)/'bfield  '/, name(89)/'G         '/, unit(89)/1.d-4/

  data dimm(90)/'byte    '/, name(90)/'byte      '/, unit(90)/1.00000000000000000e+00_dp/
  data dimm(91)/'byte    '/, name(91)/'kB        '/, unit(91)/1.d3/
  data dimm(92)/'byte    '/, name(92)/'kiB       '/, unit(92)/1.02400000000000000e+03_dp/
  data dimm(93)/'byte    '/, name(93)/'MB        '/, unit(93)/1.d6/
  data dimm(94)/'byte    '/, name(94)/'MiB       '/, unit(94)/1.04857600000000000e+06_dp/
  data dimm(95)/'byte    '/, name(95)/'GB        '/, unit(95)/1.d9/
  data dimm(96)/'byte    '/, name(96)/'GiB       '/, unit(96)/1.07374182400000000e+09_dp/
  data dimm(97)/'byte    '/, name(97)/'TB        '/, unit(97)/1.d12/
  data dimm(98)/'byte    '/, name(98)/'TiB       '/, unit(98)/1.09951162777600000e+12_dp/
  data dimm(99)/'byte    '/, name(99)/'PB        '/, unit(99)/1.d15/
  data dimm(100)/'byte    '/, name(100)/'PiB       '/, unit(100)/1.12589990684262400e+15_dp/

  data dimm(101)/'velocity'/, name(101)/'m/s       '/, unit(101)/1.d0/
  data dimm(102)/'velocity'/, name(102)/'cm/s      '/, unit(102)/1.d-2/
  data dimm(103)/'velocity'/, name(103)/'mm/s      '/, unit(103)/1.d-3/
  data dimm(104)/'velocity'/, name(104)/'nm/s      '/, unit(104)/1.d-9/
  data dimm(105)/'velocity'/, name(105)/'Ang/ns    '/, unit(105)/1.d-1/
  data dimm(106)/'velocity'/, name(106)/'Ang/ps    '/, unit(106)/1.d2/
  data dimm(107)/'velocity'/, name(107)/'Ang/fs    '/, unit(107)/1.d5/
  data dimm(108)/'velocity'/, name(108)/'Bohr/ns   '/, unit(108)/0.529177d-1/
  data dimm(109)/'velocity'/, name(109)/'Bohr/ps   '/, unit(109)/0.529177d2/
  data dimm(110)/'velocity'/, name(110)/'Bohr/fs   '/, unit(110)/0.529177d5/

  data dimm(111)/'freq    '/, name(111)/'1/s       '/, unit(111)/1.d0/
  data dimm(112)/'freq    '/, name(112)/'1/ns      '/, unit(112)/1.d9/
  data dimm(113)/'freq    '/, name(113)/'1/ps      '/, unit(113)/1.d12/
  data dimm(114)/'freq    '/, name(114)/'1/Ang     '/, unit(114)/1.d10/
  data dimm(115)/'freq    '/, name(115)/'1/Bohr    '/, unit(115)/1.88972687777435532e+10_dp/
  data dimm(116)/'freq    '/, name(116)/'1/fs      '/, unit(116)/1.d15/
  data dimm(117)/'freq    '/, name(117)/'s**-1     '/, unit(117)/1.d0/
  data dimm(118)/'freq    '/, name(118)/'s^-1      '/, unit(118)/1.d0/
  data dimm(119)/'freq    '/, name(119)/'s-1       '/, unit(119)/1.d0/
  data dimm(120)/'freq    '/, name(120)/'ns**-1    '/, unit(120)/1.d9/
  data dimm(121)/'freq    '/, name(121)/'ns^-1     '/, unit(121)/1.d9/
  data dimm(122)/'freq    '/, name(122)/'ns-1      '/, unit(122)/1.d9/
  data dimm(123)/'freq    '/, name(123)/'ps**-1    '/, unit(123)/1.d12/
  data dimm(124)/'freq    '/, name(124)/'ps^-1     '/, unit(124)/1.d12/
  data dimm(125)/'freq    '/, name(125)/'ps-1      '/, unit(125)/1.d12/
  data dimm(126)/'freq    '/, name(126)/'Ang**-1   '/, unit(126)/1.d10/
  data dimm(127)/'freq    '/, name(127)/'Ang^-1    '/, unit(127)/1.d10/
  data dimm(128)/'freq    '/, name(128)/'Ang-1     '/, unit(128)/1.d10/
  data dimm(129)/'freq    '/, name(129)/'Bohr**-1  '/, unit(129)/1.88972687777435532e+10_dp/
  data dimm(130)/'freq    '/, name(130)/'Bohr^-1   '/, unit(130)/1.88972687777435532e+10_dp/
  data dimm(131)/'freq    '/, name(131)/'Bohr-1    '/, unit(131)/1.88972687777435532e+10_dp/
  data dimm(132)/'freq    '/, name(132)/'fs**-1    '/, unit(132)/1.d15/
  data dimm(133)/'freq    '/, name(133)/'fs^-1     '/, unit(133)/1.d15/
  data dimm(134)/'freq    '/, name(134)/'fs-1      '/, unit(134)/1.d15/

contains

  subroutine inquire_unit(unit_str, stat, phys_dim, unit_name, unit_value)
    character(len=*), intent(in)   :: unit_str   ! unit specification
    character(len=*), intent(out)  :: phys_dim   ! physical dimension (e.g. 'mass')
    character(len=*), intent(out)  :: unit_name  ! unit name (e.g. 'g')
    real(dp), intent(out)          :: unit_value ! actual value (e.g. 1.e-3)
    integer, intent(out)           :: stat       ! status code

    call inquire_unit_table(unit_str, stat, phys_dim, unit_name, unit_value, &
        nu, dimm, name, unit)

  end subroutine inquire_unit

end module units_legacy_m
