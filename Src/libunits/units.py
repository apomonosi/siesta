#
# Small script to auto-generate the units table
# parameters used for conversion between units.
# This is merely due to the easy conversion of data
# in python, relative to manually adding elements in
# FORTRAN
from __future__ import print_function

import argparse
import math
from collections import OrderedDict

parser = argparse.ArgumentParser(description='Create table input for fdf.')
parser.add_argument('--table', metavar='table',  choices=['codata-2018', 'legacy'],
                    default='codata-2018',
                    help='specify which values should be used, only use "legacy" for testing purposes!')

parser.add_argument("--lower", default=False, action='store_true',
                    help="Lower the units in the output table")
parser.add_argument('--out-format', metavar='out_format',  choices=['fortran', 'latex-table'],
                    default='fortran',
                    help='''If "fortran" it will print out the data for direct use in fortran.h.

If "latex-table" a table enabled for LaTeX input is print out!''')

args = parser.parse_args()

table_suffix = None
if args.table == "legacy":
    table_suffix = "legacy"
elif args.table == "codata-2018":
    table_suffix = "codata2018"
else:
    raise ValueError("Unknown table type")


# In SI units
constants = OrderedDict()
# Speed of light
constants['c'] = 299792458.
# Avogrados constant
constants['A'] = 6.02214076e23
# Boltzmanns constant
constants['k'] = 1.380649e-23
# Plancks constant
constants['h'] = 6.62607015e-34


# Here we have all the units
units = OrderedDict()

def get_table(name):
    global table_suffix
    return globals()[f"{name}_{table_suffix}"]()


def mass_legacy():
    mass = OrderedDict()
    mass['g'] = '1.e-3'
    mass['kg'] = '1.e0'
    mass['amu'] = '1.66054e-27'
    return mass

def mass_codata2018():
    mass = OrderedDict()
    mass['g'] = '1.e-3'
    mass['kg'] = '1.'
    mass['amu'] = '1.66053906660e-27'
    mass['da'] = 11.9999999958e-3 / constants['A']
    return mass


def length_legacy():
    length = OrderedDict()
    length['m'] = '1.e0'
    length['cm'] = '1.e-2'
    length['nm'] = '1.e-9'
    length['pm'] = '1.e-12'
    length['Ang'] = '1.e-10'
    length['Bohr'] = '0.529177e-10'
    return length

def length_codata2018():
    length = OrderedDict()
    length['m'] = '1.'
    length['cm'] = '1.e-2'
    length['nm'] = '1.e-9'
    length['pm'] = '1.e-12'
    length['Ang'] = '1.e-10'
    length['Bohr'] = '0.529177210903e-10'
    return length


def energy_legacy():
    energy = OrderedDict()
    energy['J'] = '1.e0'
    energy['kJ'] = '1.e3'
    energy['erg'] = '1.e-7'
    energy['meV'] = '1.60219e-22'
    energy['eV'] = '1.60219e-19'
    energy['mRy'] = '2.17991e-21'
    energy['Ry'] = '2.17991e-18'
    energy['mHa'] = '4.35982e-21'
    energy['mHartree'] = '4.35982e-21'
    energy['Ha'] = '4.35982e-18'
    energy['Hartree'] = '4.35982e-18'
    energy['K'] = '1.38066e-23'
    energy['Kelvin'] = '1.38066e-23'
    energy['kcal/mol'] = '6.94780e-21'
    energy['kJ/mol'] = '1.6606e-21'
    energy['Hz'] = '6.6262e-34'
    energy['THz'] = '6.6262e-22'
    energy['cm-1'] = '1.986e-23'
    energy['cm^-1'] = '1.986e-23'
    energy['cm**-1'] = '1.986e-23'
    return energy

def energy_codata2018():
    energy = OrderedDict()
    energy['J'] = '1.'
    energy['kJ'] = '1.e3'
    energy['erg'] = '1.e-7'
    energy['meV'] = '1.602176634e-22'
    energy['eV'] = '1.602176634e-19'
    energy['mRy'] = '2.1798723611035e-21'
    energy['Ry'] = '2.1798723611035e-18'
    energy['mHa'] = '4.3597447222071e-21'
    energy['mHa'] = '4.3597447222071e-21'
    energy['Ha'] = '4.3597447222071e-18'
    energy['Hartree'] = '4.3597447222071e-18'
    energy['K'] = '1.380649e-23'
    energy['Kelvin'] = '1.380649e-23'
    energy['kJ/mol'] = '1.6605390671738467e-21'
    #THIS IS NOT FROM CODATA
    # This is kJ/mol * 4.184
    energy['kcal/mol'] = float(energy['kJ/mol']) * 4.184
    energy['Hz'] = '6.62607015e-34'
    energy['THz'] = '6.62607015e-22'
    energy['cm-1'] = '1.986445857e-23'
    energy['cm^-1'] = '1.986445857e-23'
    energy['cm**-1'] = '1.986445857e-23'
    return energy


def time_legacy():
    time = OrderedDict()
    time['s'] = '1.e0'
    time['ns'] = '1.e-9'
    time['ps'] = '1.e-12'
    time['fs'] = '1.e-15'
    time['min'] = '60.e0'
    time['mins'] = '60.e0'
    time['hour'] = '3600.e0'
    time['hours'] = '3600.e0'
    time['day'] = '86400.e0'
    time['days'] = '86400.e0'
    return time

def time_codata2018():
    time = OrderedDict()
    time['s'] = '1.'
    time['ns'] = '1.e-9'
    time['ps'] = '1.e-12'
    time['fs'] = '1.e-15'
    time['min'] = '60.e0'
    time['mins'] = '60.e0'
    time['hour'] = '3600.e0'
    time['hours'] = '3600.e0'
    time['day'] = '86400.e0'
    time['days'] = '86400.e0'
    return time


def force_legacy():
    force = OrderedDict()
    force['N'] = '1.e0'
    force['eV/Ang'] = '1.60219e-9'
    force['Ry/Bohr'] = '4.11943e-8'
    force['Ha/Bohr'] = '8.23886e-08'
    return force

def force_codata2018():
    global units
    energy = units["energy"]
    length = units["length"]
    force = OrderedDict()
    force['N'] = '1.'
    force['eV/Ang'] = float(energy["eV"]) / float(length["Ang"])
    force['Ry/Bohr'] = float(energy["Ry"]) / float(length["Bohr"])
    force['Ha/Bohr'] = float(energy["Ha"]) / float(length["Bohr"])
    return force


def pressure_legacy():
    global units
    pressure = OrderedDict()
    pressure['Pa'] = '1.e0'
    pressure['GPa'] = '1.e9'
    pressure['atm'] = '1.01325e5'
    pressure['bar'] = '1.e5'
    pressure['kbar'] = '1.e8'
    pressure['Mbar'] = '1.e11'
    pressure['eV/Ang**3'] = '1.60219e11'
    pressure['eV/Ang^3'] = '1.60219e11'
    pressure['Ry/Bohr**3'] = '1.47108e13'
    pressure['Ry/Bohr^3'] = '1.47108e13'
    pressure['Ha/Bohr^3'] = '2.94216e13'
    pressure['Ha/Bohr**3'] = '2.94216e13'
    return pressure

def pressure_codata2018():
    global units
    pressure = OrderedDict()
    energy = units["energy"]
    length = units["length"]
    pressure['Pa'] = '1.'
    pressure['GPa'] = '1.e9'
    pressure['atm'] = '1.01325e5'
    pressure['bar'] = '1.e5'
    pressure['kbar'] = '1.e8'
    pressure['Mbar'] = '1.e11'
    pressure['eV/Ang**3'] = float(energy["eV"]) / float(length["Ang"]) ** 3
    pressure['eV/Ang^3'] = pressure['eV/Ang**3']
    pressure['Ry/Bohr**3'] = float(energy["Ry"]) / float(length["Bohr"]) ** 3
    pressure['Ry/Bohr^3'] = pressure['Ry/Bohr**3']
    pressure['Ha/Bohr**3'] = float(energy["Ha"]) / float(length["Bohr"]) ** 3
    pressure['Ha/Bohr^3'] = pressure['Ha/Bohr**3']
    return pressure


def surftens_legacy():
    surface_tension = OrderedDict()
    # These are all define equivalently, always! :)
    surface_tension['N/m'] = '1.e0'
    surface_tension['mN/m'] = '1.e3'
    surface_tension['dyn/cm'] = '1.e3'
    surface_tension['erg/cm**2'] = '1.e3'
    return surface_tension

def surftens_codata2018():
    return surftens_legacy()


def charge_legacy():
    charge = OrderedDict()
    charge['c'] = '1.e0'
    charge['e'] = '1.602177e-19'
    return charge

def charge_codata2018():
    charge = OrderedDict()
    charge['c'] = '1.'
    charge['e'] = '1.602176634e-19'
    return charge


def dipole_legacy():
    dipole = OrderedDict()
    dipole['c*m'] = '1.e0'
    dipole['D'] = '3.33564e-30'
    dipole['Debye'] = '3.33564e-30'
    dipole['e*Bohr'] = '8.47835e-30'
    dipole['e*Ang'] = '1.602177e-29'
    return dipole

def dipole_codata2018():
    global units, constants
    dipole = OrderedDict()
    charge = units["charge"]
    length = units["length"]
    dipole['c*m'] = '1.'
    dipole['e*Bohr'] = float(charge['e']) * float(length['Bohr'])
    dipole['e*Ang'] = float(charge['e']) * float(length['Ang'])
    # D = 1e-18 statC * cm
    #    statC = 0.1 / c with c = speed of light
    dipole['D'] = 0.1 / constants['c'] * float(length['cm']) * 1e-18
    dipole['Debye'] = dipole['D']
    return dipole


def mominert_legacy():
    mominert = OrderedDict()
    mominert['kg*m**2'] = '1.e0'
    mominert['Ry*fs**2'] = '2.17991e-48'
    return mominert

def mominert_codata2018():
    global units
    energy = units["energy"]
    time = units["time"]
    mominert = OrderedDict()
    mominert['kg*m**2'] = '1.'
    mominert['Ry*fs**2'] = float(energy['Ry']) * float(time['fs']) ** 2
    return mominert


def efield_legacy():
    efield = OrderedDict()
    efield['V/m'] = '1.e0'
    efield['V/cm'] = '1.e2'
    efield['V/um'] = '1.e6'
    efield['V/nm'] = '1.e9'
    efield['V/Ang'] = '1.e10'
    efield['V/Bohr'] = '1.8897268e10'
    efield['Ry/Bohr/e'] = '2.5711273e11'
    efield['Ha/Bohr/e'] = '5.1422546e11'
    efield['Har/Bohr/e'] = '5.1422546e11'
    return efield

def efield_codata2018():
    global units
    length = units["length"]
    energy = units["energy"]
    charge = units["charge"]
    efield = OrderedDict()
    efield['V/m'] = '1.'
    efield['V/cm'] = '1.e2'
    efield['V/um'] = '1.e6'
    efield['V/nm'] = '1.e9'
    efield['V/Ang'] = '1.e10'
    efield['eV/Ang/e'] = '1.e10' # just for completeness sake
    efield['V/Bohr'] = 1. / float(length['Bohr'])
    efield['Ry/Bohr/e'] = float(energy['Ry']) / float(length['Bohr']) / float(charge["e"])
    efield['Ha/Bohr/e'] = float(energy['Ha']) / float(length['Bohr']) / float(charge["e"])
    efield['Har/Bohr/e'] = efield['Ha/Bohr/e']
    return efield


def angle_legacy():
    angle = OrderedDict()
    angle['deg'] = '1.e0'
    angle['rad'] = '5.72957795e1'
    return angle

def angle_codata2018():
    angle = OrderedDict()
    angle['deg'] = '1.'
    angle['rad'] = 180 / math.pi
    return angle


def torque_legacy():
    torque = OrderedDict()
    torque['meV/deg'] = '1.0e-3'
    torque['meV/rad'] = '1.745533e-5'
    torque['eV/deg'] = '1.0e0'
    torque['eV/rad'] = '1.745533e-2'
    torque['mRy/deg'] = '13.6058e-3'
    torque['mRy/rad'] = '0.237466e-3'
    torque['Ry/deg'] = '13.6058e0'
    torque['Ry/rad'] = '0.237466e0'
    return torque

def torque_codata2018():
    global units
    energy = units["energy"]
    angle = units["angle"]
    torque = OrderedDict()
    # NOTE
    # This is changed drastically since we now default to the SI unit system
    # This is to remain consistent with the rest of the code.
    torque['N*m'] = '1.'
    torque['meV/deg'] = energy['meV']
    torque['meV/rad'] = float(torque['meV/deg']) / float(angle['rad'])
    torque['eV/deg'] = energy['eV']
    torque['eV/rad'] = float(torque['eV/deg']) / float(angle['rad'])
    torque['mRy/deg'] = energy['mRy']
    torque['mRy/rad'] = float(torque['mRy/deg']) / float(angle['rad'])
    torque['Ry/deg'] = energy['Ry']
    torque['Ry/rad'] = float(torque['Ry/deg']) / float(angle['rad'])
    torque['Ha/deg'] = energy['Ha']
    torque['Ha/rad'] = float(torque['Ha/deg']) / float(angle['rad'])
    return torque


def velocity_legacy():
    velocity = OrderedDict()
    velocity["m/s"] = "1.d0"
    velocity["cm/s"] = "1.d-2"
    velocity["mm/s"] = "1.d-3"
    velocity["nm/s"] = "1.d-9"
    velocity["Ang/ns"] = "1.d-1"
    velocity["Ang/ps"] = "1.d2"
    velocity["Ang/fs"] = "1.d5"
    velocity["Bohr/ns"] = "0.529177d-1"
    velocity["Bohr/ps"] = "0.529177d2"
    velocity["Bohr/fs"] = "0.529177d5"
    return velocity

def velocity_codata2018():
    global units
    length = units["length"]
    time = units["time"]
    velocity = OrderedDict()
    velocity["m/s"] = "1.d0"
    velocity["cm/s"] = "1.d-2"
    velocity["mm/s"] = "1.d-3"
    velocity["nm/s"] = "1.d-9"
    velocity["Ang/ns"] = "1.d-1"
    velocity["Ang/ps"] = "1.d2"
    velocity["Ang/fs"] = "1.d5"
    velocity["Bohr/ns"] = float(length["Bohr"]) / float(time["ns"])
    velocity["Bohr/ps"] = float(length["Bohr"]) / float(time["ps"])
    velocity["Bohr/fs"] = float(length["Bohr"]) / float(time["fs"])
    return velocity


def freq_legacy():
    global units
    length = units["length"]

    freq = OrderedDict()
    freq["1/s"] = "1.d0"
    freq["1/ns"] = "1.d9"
    freq["1/ps"] = "1.d12"
    freq["1/Ang"] = "1.d10"
    freq["1/Bohr"] = 1/float(length["Bohr"])
    freq["1/fs"] = "1.d15"
    for f in list(freq.keys()):
        freq[f"{f[2:]}**-1"] = freq[f]
        freq[f"{f[2:]}^-1"] = freq[f]
        freq[f"{f[2:]}-1"] = freq[f]
    return freq

freq_codata2018 = freq_legacy


def byte_legacy():
    byte = OrderedDict()
    byte["byte"] = 1
    byte["kB"] = "1.d3"
    byte["kiB"] = 1024
    byte["MB"] = "1.d6"
    byte["MiB"] = 1024**2
    byte["GB"] = "1.d9"
    byte["GiB"] = 1024**3
    byte["TB"] = "1.d12"
    byte["TiB"] = 1024**4
    byte["PB"] = "1.d15"
    byte["PiB"] = 1024**5
    return byte

byte_codata2018 = byte_legacy


def bfield_legacy():
    bfield = OrderedDict()
    bfield["Tesla"] = '1.'
    bfield["G"] = '1.e-4'
    return bfield

def bfield_codata2018():
    return bfield_legacy()


for t in ("mass", "length",
          "energy", "time", "force",
          "pressure", "surftens", "charge", "dipole",
          "mominert", "efield", "angle", "torque",
          "bfield", "byte",
          "velocity", "freq",
        ):
    units[t] = get_table(t)


def check_ambiguity(units):
    for field1 in units:
        unit1 = units[field1]
        for name1 in unit1:
            i = 0
            lst = []
            for field2 in units:
                unit2 = units[field2]
                for name2 in unit2:
                    if name1.lower() == name2.lower():
                        i += 1
                        if i > 1:
                            lst.append((field1, name1, field2, name2))
            if i != 1:
                print(f"Unit specification is requried for {lst}")

check_ambiguity(units)

def max_length(units):
    fl = 1
    ul = 1
    for field in units:
        fl = max(fl, len(field))
        for unit in units[field]:
            ul = max(ul, len(unit))
    return fl, ul

def total_units(units):
    l = 0
    for field in units:
        l = l + len(units[field])
    return l


dimm_l, name_l = max_length(units)
max_units = 10


if args.out_format == "fortran":

    if args.lower:
        def p(unit):
            return unit.lower()
    else:
        def p(unit):
            return unit
    # Now write everything
    ind = ' ' * 2
    print(ind + '! Data generated from table: {}'.format(args.table))
    print(ind + 'integer, parameter :: nu = {}'.format(total_units(units)))
    print(ind + 'character({}), save :: dimm(nu)'.format(dimm_l))
    print(ind + 'character({}), save :: name(nu)'.format(name_l))
    print(ind + 'real(dp), save :: unit(nu)')

    fmt_s = "'{0:" + str(dimm_l) + "s}', '{1:" + str(name_l) + "s}', {2:s}"
    # Double precision has (up to) 17 significant digits.
    fmt_f = "'{0:" + str(dimm_l) + "s}', '{1:" + str(name_l) + "s}', {2:<.17e}_dp"

    def get_dict(field, unit, val):
        """ Convert to dict """
        dimm_f = "'{0:" + str(dimm_l) + "s}'"
        unit_f = "'{0:" + str(name_l) + "s}'"
        if isinstance(val, (float, int)):
            val_f = "{0:<.17e}_dp"
        else:
            val = val.replace("e", "d")
            val_f = "{0:s}"
        return {"dimm": dimm_f.format(field),
                "unit": unit_f.format(unit),
                "val": val_f.format(val)}

    N = 0
    for field in units:
        for i, unit in enumerate(units[field]):
            N += 1
            print(ind + 'data dimm({i})/{dimm}/, name({i})/{unit}/, unit({i})/{val}/'
                  .format(i=N, **get_dict(field, p(unit), units[field][unit])))

        print()

elif args.out_format == "latex-table":

    fmt = f"{{0:{dimm_l}s}} & {{1:{name_l}s}} \\\\"

    for name, table in units.items():
        for unit, value in table.items():
            print(fmt.format(name, unit.replace('^', '\^')))
