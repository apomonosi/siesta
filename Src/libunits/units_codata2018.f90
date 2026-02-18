!
! Copyright (C) 1996-2022	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!> CODATA defined units from the 2018 table
!>
!> These units are unified to follow the CODATA-2018
!> data table.
!> In addition we have some basic definitions such as
!> pi and conversion to degree (from radians).
!> The conversion factors are following this convention:
!>
!>   Ang [Ang] == [Bohr]
!>   [Ang] = [Bohr] / Ang
!>
!> meaning that the units are conversion factors from the named unit
!> to the intrinsic Siesta unit.
module units_codata2018_m

  use units_common_m

  implicit none
  private

  integer, parameter :: dp = selected_real_kind(14,100)

  real(dp), parameter, public :: pi = 3.14159265358979323846264338327950288419716939937510_dp
  real(dp), parameter, public :: pi2 = pi * 2._dp
  real(dp), parameter, public :: deg = pi / 180.0_dp

  ! These values are not present in the *legacy* scheme.
  real(dp), parameter, public :: kg       = 1.28460971172331e27_dp
  real(dp), parameter, public :: Joule    = 1._dp/2.1798723611035e-18_dp
  real(dp), parameter, public :: h_planck = 6.62607015e-19_dp * Joule ! in Ry*fs
  real(dp), parameter, public :: hbar     = h_planck / pi2 ! in Ry*fs
  real(dp), parameter, public :: c_light  = 299792458.0_dp ! in m/s

  real(dp), parameter, public :: Ang    = 1.88972612462577017_dp
  real(dp), parameter, public :: eV     = 7.34986443513115789e-2_dp
  real(dp), parameter, public :: kBar   = 6.79786184348648780e-6_dp
  real(dp), parameter, public :: Kelvin = 6.33362312691136091e-6_dp
  real(dp), parameter, public :: Debye  = 3.93430269519899511e-1_dp
  real(dp), parameter, public :: amu    = 2.13314461165032_dp ! 1.6605390666e-27 *kg
  real(dp), parameter, public :: Ryd_time = 1._dp / hbar

  real(dp), parameter, public :: GPa = kBar * 10

  public :: inquire_unit

  ! Data generated from table: codata-2018
  integer, parameter :: nu = 138
  character(8), save :: dimm(nu)
  character(10), save :: name(nu)
  real(dp), save :: unit(nu)
  data dimm(1)/'mass    '/, name(1)/'g         '/, unit(1)/1.d-3/
  data dimm(2)/'mass    '/, name(2)/'kg        '/, unit(2)/1./
  data dimm(3)/'mass    '/, name(3)/'amu       '/, unit(3)/1.66053906660d-27/
  data dimm(4)/'mass    '/, name(4)/'da        '/, unit(4)/1.99264687991118964e-26_dp/

  data dimm(5)/'length  '/, name(5)/'m         '/, unit(5)/1./
  data dimm(6)/'length  '/, name(6)/'cm        '/, unit(6)/1.d-2/
  data dimm(7)/'length  '/, name(7)/'nm        '/, unit(7)/1.d-9/
  data dimm(8)/'length  '/, name(8)/'pm        '/, unit(8)/1.d-12/
  data dimm(9)/'length  '/, name(9)/'Ang       '/, unit(9)/1.d-10/
  data dimm(10)/'length  '/, name(10)/'Bohr      '/, unit(10)/0.529177210903d-10/

  data dimm(11)/'energy  '/, name(11)/'J         '/, unit(11)/1./
  data dimm(12)/'energy  '/, name(12)/'kJ        '/, unit(12)/1.d3/
  data dimm(13)/'energy  '/, name(13)/'erg       '/, unit(13)/1.d-7/
  data dimm(14)/'energy  '/, name(14)/'meV       '/, unit(14)/1.602176634d-22/
  data dimm(15)/'energy  '/, name(15)/'eV        '/, unit(15)/1.602176634d-19/
  data dimm(16)/'energy  '/, name(16)/'mRy       '/, unit(16)/2.1798723611035d-21/
  data dimm(17)/'energy  '/, name(17)/'Ry        '/, unit(17)/2.1798723611035d-18/
  data dimm(18)/'energy  '/, name(18)/'mHa       '/, unit(18)/4.3597447222071d-21/
  data dimm(19)/'energy  '/, name(19)/'Ha        '/, unit(19)/4.3597447222071d-18/
  data dimm(20)/'energy  '/, name(20)/'Hartree   '/, unit(20)/4.3597447222071d-18/
  data dimm(21)/'energy  '/, name(21)/'K         '/, unit(21)/1.380649d-23/
  data dimm(22)/'energy  '/, name(22)/'Kelvin    '/, unit(22)/1.380649d-23/
  data dimm(23)/'energy  '/, name(23)/'kJ/mol    '/, unit(23)/1.6605390671738467d-21/
  data dimm(24)/'energy  '/, name(24)/'kcal/mol  '/, unit(24)/6.94769545705537413e-21_dp/
  data dimm(25)/'energy  '/, name(25)/'Hz        '/, unit(25)/6.62607015d-34/
  data dimm(26)/'energy  '/, name(26)/'THz       '/, unit(26)/6.62607015d-22/
  data dimm(27)/'energy  '/, name(27)/'cm-1      '/, unit(27)/1.986445857d-23/
  data dimm(28)/'energy  '/, name(28)/'cm^-1     '/, unit(28)/1.986445857d-23/
  data dimm(29)/'energy  '/, name(29)/'cm**-1    '/, unit(29)/1.986445857d-23/

  data dimm(30)/'time    '/, name(30)/'s         '/, unit(30)/1./
  data dimm(31)/'time    '/, name(31)/'ns        '/, unit(31)/1.d-9/
  data dimm(32)/'time    '/, name(32)/'ps        '/, unit(32)/1.d-12/
  data dimm(33)/'time    '/, name(33)/'fs        '/, unit(33)/1.d-15/
  data dimm(34)/'time    '/, name(34)/'min       '/, unit(34)/60.d0/
  data dimm(35)/'time    '/, name(35)/'mins      '/, unit(35)/60.d0/
  data dimm(36)/'time    '/, name(36)/'hour      '/, unit(36)/3600.d0/
  data dimm(37)/'time    '/, name(37)/'hours     '/, unit(37)/3600.d0/
  data dimm(38)/'time    '/, name(38)/'day       '/, unit(38)/86400.d0/
  data dimm(39)/'time    '/, name(39)/'days      '/, unit(39)/86400.d0/

  data dimm(40)/'force   '/, name(40)/'N         '/, unit(40)/1./
  data dimm(41)/'force   '/, name(41)/'eV/Ang    '/, unit(41)/1.60217663399999979e-09_dp/
  data dimm(42)/'force   '/, name(42)/'Ry/Bohr   '/, unit(42)/4.11936174912694464e-08_dp/
  data dimm(43)/'force   '/, name(43)/'Ha/Bohr   '/, unit(43)/8.23872349825407855e-08_dp/

  data dimm(44)/'pressure'/, name(44)/'Pa        '/, unit(44)/1./
  data dimm(45)/'pressure'/, name(45)/'GPa       '/, unit(45)/1.d9/
  data dimm(46)/'pressure'/, name(46)/'atm       '/, unit(46)/1.01325d5/
  data dimm(47)/'pressure'/, name(47)/'bar       '/, unit(47)/1.d5/
  data dimm(48)/'pressure'/, name(48)/'kbar      '/, unit(48)/1.d8/
  data dimm(49)/'pressure'/, name(49)/'Mbar      '/, unit(49)/1.d11/
  data dimm(50)/'pressure'/, name(50)/'eV/Ang**3 '/, unit(50)/1.60217663399999969e+11_dp/
  data dimm(51)/'pressure'/, name(51)/'eV/Ang^3  '/, unit(51)/1.60217663399999969e+11_dp/
  data dimm(52)/'pressure'/, name(52)/'Ry/Bohr**3'/, unit(52)/1.47105078482607109e+13_dp/
  data dimm(53)/'pressure'/, name(53)/'Ry/Bohr^3 '/, unit(53)/1.47105078482607109e+13_dp/
  data dimm(54)/'pressure'/, name(54)/'Ha/Bohr**3'/, unit(54)/2.94210156965220977e+13_dp/
  data dimm(55)/'pressure'/, name(55)/'Ha/Bohr^3 '/, unit(55)/2.94210156965220977e+13_dp/

  data dimm(56)/'surftens'/, name(56)/'N/m       '/, unit(56)/1.d0/
  data dimm(57)/'surftens'/, name(57)/'mN/m      '/, unit(57)/1.d3/
  data dimm(58)/'surftens'/, name(58)/'dyn/cm    '/, unit(58)/1.d3/
  data dimm(59)/'surftens'/, name(59)/'erg/cm**2 '/, unit(59)/1.d3/

  data dimm(60)/'charge  '/, name(60)/'c         '/, unit(60)/1./
  data dimm(61)/'charge  '/, name(61)/'e         '/, unit(61)/1.602176634d-19/

  data dimm(62)/'dipole  '/, name(62)/'c*m       '/, unit(62)/1./
  data dimm(63)/'dipole  '/, name(63)/'e*Bohr    '/, unit(63)/8.47835362554076597e-30_dp/
  data dimm(64)/'dipole  '/, name(64)/'e*Ang     '/, unit(64)/1.60217663400000003e-29_dp/
  data dimm(65)/'dipole  '/, name(65)/'D         '/, unit(65)/3.33564095198152075e-30_dp/
  data dimm(66)/'dipole  '/, name(66)/'Debye     '/, unit(66)/3.33564095198152075e-30_dp/

  data dimm(67)/'mominert'/, name(67)/'kg*m**2   '/, unit(67)/1./
  data dimm(68)/'mominert'/, name(68)/'Ry*fs**2  '/, unit(68)/2.17987236110350002e-48_dp/

  data dimm(69)/'efield  '/, name(69)/'V/m       '/, unit(69)/1./
  data dimm(70)/'efield  '/, name(70)/'V/cm      '/, unit(70)/1.d2/
  data dimm(71)/'efield  '/, name(71)/'V/um      '/, unit(71)/1.d6/
  data dimm(72)/'efield  '/, name(72)/'V/nm      '/, unit(72)/1.d9/
  data dimm(73)/'efield  '/, name(73)/'V/Ang     '/, unit(73)/1.d10/
  data dimm(74)/'efield  '/, name(74)/'eV/Ang/e  '/, unit(74)/1.d10/
  data dimm(75)/'efield  '/, name(75)/'V/Bohr    '/, unit(75)/1.88972612462577019e+10_dp/
  data dimm(76)/'efield  '/, name(76)/'Ry/Bohr/e '/, unit(76)/2.57110337381623871e+11_dp/
  data dimm(77)/'efield  '/, name(77)/'Ha/Bohr/e '/, unit(77)/5.14220674763259521e+11_dp/
  data dimm(78)/'efield  '/, name(78)/'Har/Bohr/e'/, unit(78)/5.14220674763259521e+11_dp/

  data dimm(79)/'angle   '/, name(79)/'deg       '/, unit(79)/1./
  data dimm(80)/'angle   '/, name(80)/'rad       '/, unit(80)/5.72957795130823229e+01_dp/

  data dimm(81)/'torque  '/, name(81)/'N*m       '/, unit(81)/1./
  data dimm(82)/'torque  '/, name(82)/'meV/deg   '/, unit(82)/1.602176634d-22/
  data dimm(83)/'torque  '/, name(83)/'meV/rad   '/, unit(83)/2.79632574618201289e-24_dp/
  data dimm(84)/'torque  '/, name(84)/'eV/deg    '/, unit(84)/1.602176634d-19/
  data dimm(85)/'torque  '/, name(85)/'eV/rad    '/, unit(85)/2.79632574618201262e-21_dp/
  data dimm(86)/'torque  '/, name(86)/'mRy/deg   '/, unit(86)/2.1798723611035d-21/
  data dimm(87)/'torque  '/, name(87)/'mRy/rad   '/, unit(87)/3.80459499744788459e-23_dp/
  data dimm(88)/'torque  '/, name(88)/'Ry/deg    '/, unit(88)/2.1798723611035d-18/
  data dimm(89)/'torque  '/, name(89)/'Ry/rad    '/, unit(89)/3.80459499744788462e-20_dp/
  data dimm(90)/'torque  '/, name(90)/'Ha/deg    '/, unit(90)/4.3597447222071d-18/
  data dimm(91)/'torque  '/, name(91)/'Ha/rad    '/, unit(91)/7.60918999489594379e-20_dp/

  data dimm(92)/'bfield  '/, name(92)/'Tesla     '/, unit(92)/1./
  data dimm(93)/'bfield  '/, name(93)/'G         '/, unit(93)/1.d-4/

  data dimm(94)/'byte    '/, name(94)/'byte      '/, unit(94)/1.00000000000000000e+00_dp/
  data dimm(95)/'byte    '/, name(95)/'kB        '/, unit(95)/1.d3/
  data dimm(96)/'byte    '/, name(96)/'kiB       '/, unit(96)/1.02400000000000000e+03_dp/
  data dimm(97)/'byte    '/, name(97)/'MB        '/, unit(97)/1.d6/
  data dimm(98)/'byte    '/, name(98)/'MiB       '/, unit(98)/1.04857600000000000e+06_dp/
  data dimm(99)/'byte    '/, name(99)/'GB        '/, unit(99)/1.d9/
  data dimm(100)/'byte    '/, name(100)/'GiB       '/, unit(100)/1.07374182400000000e+09_dp/
  data dimm(101)/'byte    '/, name(101)/'TB        '/, unit(101)/1.d12/
  data dimm(102)/'byte    '/, name(102)/'TiB       '/, unit(102)/1.09951162777600000e+12_dp/
  data dimm(103)/'byte    '/, name(103)/'PB        '/, unit(103)/1.d15/
  data dimm(104)/'byte    '/, name(104)/'PiB       '/, unit(104)/1.12589990684262400e+15_dp/

  data dimm(105)/'velocity'/, name(105)/'m/s       '/, unit(105)/1.d0/
  data dimm(106)/'velocity'/, name(106)/'cm/s      '/, unit(106)/1.d-2/
  data dimm(107)/'velocity'/, name(107)/'mm/s      '/, unit(107)/1.d-3/
  data dimm(108)/'velocity'/, name(108)/'nm/s      '/, unit(108)/1.d-9/
  data dimm(109)/'velocity'/, name(109)/'Ang/ns    '/, unit(109)/1.d-1/
  data dimm(110)/'velocity'/, name(110)/'Ang/ps    '/, unit(110)/1.d2/
  data dimm(111)/'velocity'/, name(111)/'Ang/fs    '/, unit(111)/1.d5/
  data dimm(112)/'velocity'/, name(112)/'Bohr/ns   '/, unit(112)/5.29177210902999975e-02_dp/
  data dimm(113)/'velocity'/, name(113)/'Bohr/ps   '/, unit(113)/5.29177210902999988e+01_dp/
  data dimm(114)/'velocity'/, name(114)/'Bohr/fs   '/, unit(114)/5.29177210902999941e+04_dp/

  data dimm(115)/'freq    '/, name(115)/'1/s       '/, unit(115)/1.d0/
  data dimm(116)/'freq    '/, name(116)/'1/ns      '/, unit(116)/1.d9/
  data dimm(117)/'freq    '/, name(117)/'1/ps      '/, unit(117)/1.d12/
  data dimm(118)/'freq    '/, name(118)/'1/Ang     '/, unit(118)/1.d10/
  data dimm(119)/'freq    '/, name(119)/'1/Bohr    '/, unit(119)/1.88972612462577019e+10_dp/
  data dimm(120)/'freq    '/, name(120)/'1/fs      '/, unit(120)/1.d15/
  data dimm(121)/'freq    '/, name(121)/'s**-1     '/, unit(121)/1.d0/
  data dimm(122)/'freq    '/, name(122)/'s^-1      '/, unit(122)/1.d0/
  data dimm(123)/'freq    '/, name(123)/'s-1       '/, unit(123)/1.d0/
  data dimm(124)/'freq    '/, name(124)/'ns**-1    '/, unit(124)/1.d9/
  data dimm(125)/'freq    '/, name(125)/'ns^-1     '/, unit(125)/1.d9/
  data dimm(126)/'freq    '/, name(126)/'ns-1      '/, unit(126)/1.d9/
  data dimm(127)/'freq    '/, name(127)/'ps**-1    '/, unit(127)/1.d12/
  data dimm(128)/'freq    '/, name(128)/'ps^-1     '/, unit(128)/1.d12/
  data dimm(129)/'freq    '/, name(129)/'ps-1      '/, unit(129)/1.d12/
  data dimm(130)/'freq    '/, name(130)/'Ang**-1   '/, unit(130)/1.d10/
  data dimm(131)/'freq    '/, name(131)/'Ang^-1    '/, unit(131)/1.d10/
  data dimm(132)/'freq    '/, name(132)/'Ang-1     '/, unit(132)/1.d10/
  data dimm(133)/'freq    '/, name(133)/'Bohr**-1  '/, unit(133)/1.88972612462577019e+10_dp/
  data dimm(134)/'freq    '/, name(134)/'Bohr^-1   '/, unit(134)/1.88972612462577019e+10_dp/
  data dimm(135)/'freq    '/, name(135)/'Bohr-1    '/, unit(135)/1.88972612462577019e+10_dp/
  data dimm(136)/'freq    '/, name(136)/'fs**-1    '/, unit(136)/1.d15/
  data dimm(137)/'freq    '/, name(137)/'fs^-1     '/, unit(137)/1.d15/
  data dimm(138)/'freq    '/, name(138)/'fs-1      '/, unit(138)/1.d15/

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

end module units_codata2018_m
