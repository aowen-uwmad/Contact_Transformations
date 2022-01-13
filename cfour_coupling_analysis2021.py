#!/home/aowen4/python3/bin/python3

# Updated coupling analysis for predicting GaJ, etc. using CFOUR, coded in python3,
# using results of automated derivation of Hamiltonian reduction.

import sys
import numpy as np

if len(sys.argv) == 1:
    print("""
    
        |   Usage: cfour_coupling_analysis.py arguments filename
        |
        |   Optional arguments:
        |       -c[utoff] N (considers only vibrational states below N cm^-1)
        |       -d[ifference] M (considers only vibrational states with frequencies differing by less than M cm^-1)
        |
        |   This script parses CFOUR anharmonic vibrational analysis results to do high order vibrot coupling
        |   predictions and will write the comma-delimited data to a "_coupling.csv" file.
        
    """)
    quit()

freq_cutoff = None
max_freq_diff = None
file_indices = range(1, len(sys.argv))
for flag in ['-c', '-cutoff', '-d', '-difference']:
    if flag in sys.argv:
        flag_pos = sys.argv.index(flag)
        arg_pos = flag_pos + 1
        if flag in ['-c', '-cutoff']:
            freq_cutoff = float(sys.argv[arg_pos])
        elif flag in ['-d', '-difference']:
            max_freq_diff = float(sys.argv[arg_pos])
        file_indices = [x for x in file_indices if x != flag_pos and x != arg_pos]
file_list = [sys.argv[x] for x in file_indices]

if freq_cutoff is None:
    freq_cutoff = 200.0  # cm^-1
if max_freq_diff is None:
    max_freq_diff = 50.0  # cm^-1 

test = str(input('...Did you remember to use `qsub -I -l nodes=1:ppn=1`? (y/n)...')).lower()
if test not in ['y', 'yes']:
    print('Try again after entering the qsub command.')
    quit()

if len(file_list) == 0:
    print('\n*** Filename not provided, output will be saved as `coupling-analysis_coupling.csv` ***\n')
    out_file_name = 'coupling_analysis'
else:
    out_file_name = file_list[0]

# file imports
try:
    with open('anharm.out', 'r') as anharm_in_raw:
        anharm_in = anharm_in_raw.read()
    with open('proc/corioliszeta', 'r') as zetas_in_raw:
        zetas_in = zetas_in_raw.read().split('Coriolis')
    with open('proc/cubic', 'r') as cubic_in_raw:
        cubic_in = cubic_in_raw.read().split('\n')
    with open('didQ', 'r') as didQ_in_raw:
        didQ_in = didQ_in_raw.read().split('\n')
except IOError:
    print("""
    !!! Could not find the requisite files. !!!
        Check that the following files are present:
            anharm.out
            didQ
            proc/
                corioliszeta
                cubic
    
    """)
    quit()

# --------------------- #
# Fundamental constants #
# --------------------- #
h_not_bar = 6.626 * 10 ** (-34)  # kg*m**2/s
h_bar = h_not_bar / (2 * np.pi)
speed_of_light = 2.998 * 10 ** 8  # m/s

# -------------------------------------- #
# Importing & formatting requisite data. #
# -------------------------------------- #

# parsing anharm.out
rot_in = anharm_in.split('Be, B0 AND B-B0 SHIFTS FOR SINGLY EXCITED VIBRATIONAL STATES (CM-1)')[1].split(
    'Vibrationally averaged dipole moment')[0].split('\n')[4].split()[1:]
rot_formed = [float(x) for x in rot_in]
vib_in = anharm_in.split('Rotationally projected vibrational frequencies')[-1].split(
    'Finite temperature thermodynamic corrections')[0].split('\n')[7:-3]
# first 6 "vibrations" in CFOUR are the rotation and translations with frequencies of 0. Using dummy energies for those.
fun_vibs_formed = [1000000.0] * 6
fun_vibs_formed.extend([float(vibration.split()[1]) for vibration in vib_in])
n_modes = len(fun_vibs_formed)

# parsing corioliszeta
x_zetas_in = zetas_in[1].split('\n')[1:-1]
y_zetas_in = zetas_in[2].split('\n')[1:-1]
z_zetas_in = zetas_in[3].split('\n')[1:-1]
x_zetas = [[[int(x.split()[0]) - 1, int(x.split()[1]) - 1], float(x.split()[2])] for x in x_zetas_in]
y_zetas = [[[int(x.split()[0]) - 1, int(x.split()[1]) - 1], float(x.split()[2])] for x in y_zetas_in]
z_zetas = [[[int(x.split()[0]) - 1, int(x.split()[1]) - 1], float(x.split()[2])] for x in z_zetas_in]
# Starting with null data array of proper size
x_zetas_out = [[0 for x in range(0, n_modes)] for y in range(0, n_modes)]
y_zetas_out = [[0 for x in range(0, n_modes)] for y in range(0, n_modes)]
z_zetas_out = [[0 for x in range(0, n_modes)] for y in range(0, n_modes)]
for x in range(0, len(x_zetas)):
    x_zetas_out[x_zetas[x][0][0]][x_zetas[x][0][1]] = x_zetas[x][1]
    x_zetas_out[x_zetas[x][0][1]][x_zetas[x][0][0]] = x_zetas[x][1] * (-1)
for x in range(0, len(y_zetas)):
    y_zetas_out[y_zetas[x][0][0]][y_zetas[x][0][1]] = y_zetas[x][1]
    y_zetas_out[y_zetas[x][0][1]][y_zetas[x][0][0]] = y_zetas[x][1] * (-1)
for x in range(0, len(z_zetas)):
    z_zetas_out[z_zetas[x][0][0]][z_zetas[x][0][1]] = z_zetas[x][1]
    z_zetas_out[z_zetas[x][0][1]][z_zetas[x][0][0]] = z_zetas[x][1] * (-1)
zetas_formed = [x_zetas_out, y_zetas_out, z_zetas_out]

# parsing cubic
cubic = [
    [[int(x.split()[0]) - 1, int(x.split()[1]) - 1, int(x.split()[2]) - 1], float(x.split()[3])] for x in cubic_in if
    x != '']
# Starting with null data array of proper size
cubic_out = [[[0 for x in range(0, n_modes)] for y in range(0, n_modes)] for z in range(0, n_modes)]
for x in range(0, len(cubic)):
    cubic_out[cubic[x][0][0]][cubic[x][0][1]][cubic[x][0][2]] = cubic[x][1]
    cubic_out[cubic[x][0][0]][cubic[x][0][2]][cubic[x][0][1]] = cubic[x][1]
    cubic_out[cubic[x][0][1]][cubic[x][0][0]][cubic[x][0][2]] = cubic[x][1]
    cubic_out[cubic[x][0][1]][cubic[x][0][2]][cubic[x][0][0]] = cubic[x][1]
    cubic_out[cubic[x][0][2]][cubic[x][0][0]][cubic[x][0][1]] = cubic[x][1]
    cubic_out[cubic[x][0][2]][cubic[x][0][1]][cubic[x][0][0]] = cubic[x][1]
cubic_fc_formed = cubic_out

# parsing didq
didQ = [[[int(didQ_in[x].split()[0]) - 1, int(didQ_in[x].split()[1]) - 1, int(didQ_in[x].split()[2]) - 1],
         float(didQ_in[x].split()[3])] for x in range(0, len(didQ_in)) if didQ_in[x] != '']
# Making full inertia derivatives data array
didQ_out = [[[0 for x in range(0, 3)] for y in range(0, 3)] for z in range(0, n_modes)]
for x in range(0, len(didQ)):
    didQ_out[didQ[x][0][2]][didQ[x][0][0]][didQ[x][0][1]] = didQ[x][1]
inertia_der_formed = didQ_out


# ---------------------- #
# Simplifying data calls #
# ---------------------- #
# vibrational energies in wavenumbers
def w(vib_index):
    # wavenumbers
    return fun_vibs_formed[vib_index]


# vibrational energies in hertz
def v(vib_index):
    # wavenumbers
    old_units = fun_vibs_formed[vib_index]
    # hertz
    new_units = old_units * (100 / 1) * speed_of_light
    return new_units


# coriolis zetas
def zeta(rot_index, vib_index1, vib_index2):
    return zetas_formed[rot_index][vib_index1][vib_index2]  # unitless


# rotational constants
b = {i: rot_formed[i] for i in [0, 1, 2]}  # wavenumbers
B = {i: (rot_formed[i] * 100 * speed_of_light) for i in [0, 1, 2]}  # hertz


# cubic force constants
def k3(vib_index1, vib_index2, vib_index3):
    return cubic_fc_formed[vib_index1][vib_index2][vib_index3] * 100 * speed_of_light  # hertz


# inertial derivatives
def aD(rot_index1, rot_index2, vib_index):
    value = inertia_der_formed[vib_index][rot_index1][rot_index2]  # amu**(1/2)*Bohr
    return value*(1.660529*10**(-27)/1)**(1/2)*(5.291*10**(-11))  # kg**(1/2)*m


# --------------------- #
# Fundamental functions #
# --------------------- #
# moments of inertia
def i0(rot_index):
    return h_not_bar / (8 * (np.pi ** 2) * B[rot_index])  # kg*m**2


def omega(vib_index1, vib_index2, vib_index3):
    value = (v(vib_index1) + v(vib_index2) + v(vib_index3)) * \
            (-v(vib_index1) + v(vib_index2) + v(vib_index3)) * \
            (v(vib_index1) - v(vib_index2) + v(vib_index3)) * \
            (v(vib_index1) + v(vib_index2) - v(vib_index3))
    return value  # hertz**4


def baan(rot_index1, rot_index2, vib_index1):
    # hertz
    value = ((-1) * (h_bar ** 3) / (2 * (h_not_bar ** ( 3 / 2)))) * (
            aD(rot_index1, rot_index2, vib_index1) / (i0(rot_index1) * i0(rot_index2) * (v(vib_index1) ** (1 / 2))))
    return value


def baann(rot_index1, rot_index2, vib_index1, vib_index2):
    # hertz
    value = (3 / 8) * sum([((1 / B[rotIn3]) * (
                baan(rot_index1, rotIn3, vib_index1) * baan(rotIn3, rot_index2, vib_index2) +
                baan(rot_index1, rotIn3, vib_index2) * baan(rotIn3, rot_index2, vib_index1)
        )) for rotIn3 in [0, 1, 2]])
    return value


def ee(rot_index1, rot_index2, rot_index3):
    # Levi-Cevita symbol (unitless)
    if rot_index1 != rot_index2 and rot_index1 != rot_index3 and rot_index2 != rot_index3:
        if (rot_index1 == 0 and rot_index2 == 1 and rot_index3 == 2) or \
                (rot_index1 == 1 and rot_index2 == 2 and rot_index3 == 0) or \
                (rot_index1 == 2 and rot_index2 == 0 and rot_index3 == 1):
            value = 1
        else:
            value = -1
    else:
        value = 0
    return value


def tau(rot_index1, rot_index2, rot_index3, rot_index4):
    value = -2 * sum(
        [(baan(rot_index1, rot_index2, vib_index)
          * baan(rot_index3, rot_index4, vib_index)
          / v(vib_index)) for vib_index in range(6, n_modes)])
    return value  # hertz


# --------------------------------------- #
# Vibration-rotation interaction (alphas) #
# --------------------------------------- #
def alpha_vib(rot_index, vib_index):
    value = (-2) * pow(B[rot_index], 2) * pow(v(vib_index), -1) * (
            sum([(3 * pow(aD(rot_index, rot_index2, vib_index), 2) / (4 * i0(rot_index2))) 
                 for rot_index2 in [0, 1, 2]])
            + sum([(pow(zeta(rot_index, vib_index, vib_index2), 2) 
                    * (3 * pow(v(vib_index), 2) + pow(v(vib_index2), 2)) 
                    / (pow(v(vib_index), 2) - pow(v(vib_index2), 2))) 
                   for vib_index2 in range(6, n_modes) if vib_index2 != vib_index])
            + sum([(k3(vib_index, vib_index, vib_index2)
                    * aD(rot_index, rot_index, vib_index2)
                    * v(vib_index) * pow(v(vib_index2), (-3 / 2)))
                   for vib_index2 in range(6, n_modes)]) * np.pi * pow(h_not_bar, (-1 / 2)))
    return value


def alpha_sum(rot_index, quanta):
    value = sum([(alpha_vib(rot_index, vibIn) * (quanta[vibIn] + (1 / 2))) for vibIn in range(6, n_modes)])
    return value


def deperturbed_alpha_vib(rot_index, vib_index, coriolis_pairs, fermi_pairs):
    value = alpha_vib(rot_index, vib_index)
    for coriolis_pair in coriolis_pairs:
        c0 = coriolis_pair[0]
        c1 = coriolis_pair[1]
        if vib_index == c0:
            if c0 != c1:
                value = value - (-(2 * (B[rot_index] ** 2) / v(c0))
                                 * zeta(rot_index, c0, c1) ** 2
                                 * (3 * v(c0) ** 2 + v(c1) ** 2)
                                 * ((v(c0) ** 2 - v(c1) ** 2) ** (-1)))
            value = value + (-(2 * (B[rot_index] ** 2) / v(c0))
                             * (-zeta(rot_index, c0, c1) ** 2)
                             * (B[rot_index] / v(c0))
                             * (v(c0) - v(c1)) ** 2
                             / (v(c1) * (v(c0) + v(c1))))
        elif vib_index == c1:
            if c0 != c1:
                value = value - (-(2 * (B[rot_index] ** 2) / v(c1))
                                 * zeta(rot_index, c1, c0) ** 2
                                 * (3 * v(c1) ** 2 + v(c0) ** 2)
                                 * ((v(c1) ** 2 - v(c0) ** 2) ** (-1)))
            value = value + (-(2 * (B[rot_index] ** 2) / v(c1))
                             * (-zeta(rot_index, c1, c0) ** 2)
                             * (B[rot_index] / v(c1))
                             * (v(c1) - v(c0)) ** 2
                             / (v(c0) * (v(c1) + v(c0))))
    for fermi_pair in fermi_pairs:
        f0 = fermi_pair[0]
        f1 = fermi_pair[1]
        if vib_index == f0:
            value = value - (-(2 * B[rot_index] ** 2 / v(f0))
                             * np.pi * h_not_bar ** (-1 / 2)
                             * k3(f0, f0, f1)
                             * aD(rot_index, rot_index, f1)
                             * v(f0) * v(f1) ** (-3 / 2))
        elif vib_index == f1:
            value = value - (-(2 * B[rot_index] ** 2 / v(f1))
                             * np.pi * h_not_bar ** (-1 / 2)
                             * k3(f1, f1, f0)
                             * aD(rot_index, rot_index, f0)
                             * v(f1) * v(f0) ** (-3 / 2))
    return value


def deperturbed_alpha_sum(rot_index, quanta, coriolis_pairs, fermi_pairs):
    value = sum([(deperturbed_alpha_vib(rot_index, vib_index, coriolis_pairs, fermi_pairs)
                  * (quanta[vib_index] + (1 / 2)))
                 for vib_index in range(6, n_modes)])
    return value


# ----------------------------------- #
# Couping constant definitions (t) #
# ----------------------------------- #
# *** Note: t is being used here in a generic sense, w.r.t how the Hamiltonian can be expressed as
#       H = Sum([t[n]*J[n] for n in [x,y,z]]) + Sum([t[n,m]*J[n].J[m] for n in [x,y,z] for m in [x,y,z]]) + ...
#     The exact nature of each order of t is dependent on what Hamiltonian is being written, specifically
#     whether the terms are vibrationally on-diagonal (i.e. centrifugal distortion) or off-diagonal (i.e.
#     coriolis coupling).

# CORIOLIS COUPLING DEFINITIONS
# First order (t[i]), h21
def h21(rot_index, vib_index1, vib_index2, vib_quanta1, vib_quanta2):
    value = (B[rot_index]
             * zeta(rot_index, vib_index1, vib_index2)
             * ((v(vib_index1) + v(vib_index2)) / ((v(vib_index1) * v(vib_index2)) ** (1 / 2)))
             * (((vib_quanta1 + 1) * vib_quanta2) ** (1 / 2)))
    return value  # hertz


# Second order (t[i,j]), h22
def h22_sum_A(rot_index1, rot_index2, vib_index1, vib_index2):
    value = sum([(zeta(rot_index1, vib_index1, vib_index3)
                  * zeta(rot_index2, vib_index2, vib_index3)
                  * (v(vib_index1) ** 2 + v(vib_index3) ** 2)
                  * ((v(vib_index1) ** 2 - v(vib_index3) ** 2) ** (-1))
                  + zeta(rot_index1, vib_index2, vib_index3)
                  * zeta(rot_index2, vib_index1, vib_index3)
                  * (v(vib_index2) ** 2 + v(vib_index3) ** 2)
                  * ((v(vib_index2) ** 2 - v(vib_index3) ** 2) ** (-1)))
                 for vib_index3 in range(6, n_modes) if (vib_index3 != vib_index1 and vib_index3 != vib_index2)])
    return value


def h22_sum_B(rot_index1, rot_index2, vib_index1, vib_index2):
    value = sum([(zeta(rot_index1, vib_index1, vib_index3)
                  * zeta(rot_index2, vib_index2, vib_index3)
                  * ((v(vib_index1) ** 2 - v(vib_index3) ** 2) ** (-1))
                  + zeta(rot_index1, vib_index2, vib_index3)
                  * zeta(rot_index2, vib_index1, vib_index3)
                  * ((v(vib_index2) ** 2 - v(vib_index3) ** 2) ** (-1)))
                 for vib_index3 in range(6, n_modes) if(vib_index3 != vib_index1 and vib_index3 != vib_index2)])
    return value


def h22(rot_index1, rot_index2, vib_index1, vib_index2, vib_quanta1, vib_quanta2):
    term1 = baann(rot_index1, rot_index2, vib_index1, vib_index2)
    term2 = -(1 / 4) * sum([(k3(vib_index1, vib_index2, vib_index3)
                             * baan(rot_index1, rot_index2, vib_index3)
                             * (v(vib_index3)
                                * (v(vib_index1) ** 2 + v(vib_index2) ** 2 - v(vib_index3) ** 2)
                                / omega(vib_index1, vib_index2, vib_index3)
                                + (1 / v(vib_index3))))
                            for vib_index3 in range(6, n_modes)])
    term3 = (B[rot_index1]
             * B[rot_index2]
             * ((v(vib_index1) * v(vib_index2)) ** (-1 / 2))
             * h22_sum_A(rot_index1, rot_index2, vib_index1, vib_index2))
    term4 = -(1 / 4) * sum([(k3(vib_index1, vib_index2, vib_index3)
                             * baan(rot_index1, rot_index2, vib_index3)
                             * v(vib_index1)
                             * v(vib_index2)
                             * v(vib_index3)
                             / omega(vib_index1, vib_index2, vib_index3))
                            for vib_index3 in range(6, n_modes)])
    term5 = (2
             * B[rot_index1]
             * B[rot_index2]
             * ((v(vib_index1) * v(vib_index2)) ** (1 / 2))
             * h22_sum_B(rot_index1, rot_index2, vib_index1, vib_index2))
    value = (term1 + term2 + term3 + term4 + term5) * (1 / 2) * (((vib_quanta1 + 1) * vib_quanta2) ** (1 / 2))
    return value


# Third order (t[i,j,k]), h23
def h23sumA(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2):
    value = sum([(((v(vib_index1) * v(vib_index3)) ** (1 / 2))
                  * ((v(vib_index1) ** 2 - v(vib_index3) ** 2) ** (-1))
                  * sum([((B[rot_index2] / B[rot_index4])
                          * (baan(rot_index1, rot_index4, vib_index1)
                             * zeta(rot_index2, vib_index2, vib_index3)
                             * baan(rot_index4, rot_index3, vib_index3)
                             + baan(rot_index3, rot_index4, vib_index1)
                             * zeta(rot_index2, vib_index2, vib_index3)
                             * baan(rot_index4, rot_index1, vib_index3)
                             - baan(rot_index1, rot_index4, vib_index2)
                             * zeta(rot_index2, vib_index1, vib_index3)
                             * baan(rot_index4, rot_index3, vib_index3)
                             - baan(rot_index3, rot_index4, vib_index2)
                             * zeta(rot_index2, vib_index1, vib_index3)
                             * baan(rot_index4, rot_index1, vib_index3)))
                         for rot_index4 in [0, 1, 2]]))
                 for vib_index3 in range(6, n_modes) if vib_index3 != vib_index1])
    return value


def h23sumB(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2):
    value = sum([
        sum([((2 * v(vib_index1) ** 2 + v(vib_index3) ** 2 + v(vib_index4) ** 2)
              * (zeta(rot_index1, vib_index1, vib_index3)
                 * zeta(rot_index3, vib_index2, vib_index4)
                 - zeta(rot_index1, vib_index2, vib_index3)
                 * zeta(rot_index3, vib_index1, vib_index4))
              * zeta(rot_index2, vib_index3, vib_index4)
              / ((v(vib_index1) ** 2 - v(vib_index3) ** 2)
                 * (v(vib_index1) ** 2 - v(vib_index4) ** 2)))
             for vib_index4 in range(6, n_modes) if vib_index4 != vib_index1])
        for vib_index3 in range(6, n_modes) if vib_index3 != vib_index1])
    return value


def h23sumC(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2):
    value = sum([((v(vib_index1) ** 2 + v(vib_index3) ** 2)
                  * ((v(vib_index1) ** 2 - v(vib_index3) ** 2) ** (-2))
                  * (zeta(rot_index1, vib_index1, vib_index1)
                     * zeta(rot_index3, vib_index2, vib_index3)
                     - zeta(rot_index1, vib_index2, vib_index1)
                     * zeta(rot_index3, vib_index1, vib_index3)
                     + zeta(rot_index3, vib_index1, vib_index1)
                     * zeta(rot_index1, vib_index2, vib_index3)
                     - zeta(rot_index3, vib_index2, vib_index1)
                     * zeta(rot_index1, vib_index1, vib_index3))
                  * zeta(rot_index2, vib_index1, vib_index3))
                 for vib_index3 in range(6, n_modes) if vib_index3 != vib_index1])
    return value


def h23sumD(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2):
    value = sum([((v(vib_index1) ** 2 + v(vib_index3) ** 2)
                  * ((v(vib_index1) ** 2 - v(vib_index3) ** 2) ** (-2))
                  * sum([(((B[rot_index1] - B[rot_index3])
                           * ee(rot_index1, rot_index3, rot_index4)
                           * B[rot_index4]
                           * B[rot_index2]
                           * ((zeta(rot_index4, vib_index1, vib_index3)
                               * zeta(rot_index2, vib_index2, vib_index3))
                              - (zeta(rot_index4, vib_index2, vib_index3)
                                 * zeta(rot_index2, vib_index1, vib_index3)))))
                         for rot_index4 in [0, 1, 2]]))
                 for vib_index3 in range(6, n_modes) if vib_index3 != vib_index1])
    return value


def h23sumE(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2):
    value = sum([
        sum([(((v(vib_index1) * v(vib_index4)) ** (1 / 2))
              * ((v(vib_index1) * v(vib_index1)) - (v(vib_index4) * v(vib_index4))) ** (-1)
              * ((v(vib_index3)) ** (-1))
              * baan(rot_index1, rot_index3, vib_index3)
              * B[rot_index2] *
              ((k3(vib_index3, vib_index4, vib_index1)
                * zeta(rot_index2, vib_index4, vib_index2))
               - (k3(vib_index3, vib_index4, vib_index2)
                  * zeta(rot_index2, vib_index4, vib_index1))))
             for vib_index4 in range(6, n_modes) if vib_index4 != vib_index1])
        for vib_index3 in range(6, n_modes)])
    return value


def h23(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2, vib_quanta1, vib_quanta2):
    h23A = (-1) * sum([(tau(rot_index1, rot_index3, rot_index2, rot_index4)
                        * zeta(rot_index4, vib_index1, vib_index2))
                       for rot_index4 in [0, 1, 2]])
    h23B = sum([((v(vib_index1) ** (1 / 2))
                 * (v(vib_index3) ** (-3 / 2))
                 * baan(rot_index1, rot_index3, vib_index3)
                 * sum([(baan(rot_index2, rot_index4, vib_index1)
                         * zeta(rot_index4, vib_index3, vib_index2)
                         - baan(rot_index2, rot_index4, vib_index2)
                         * zeta(rot_index4, vib_index3, vib_index1))
                        for rot_index4 in [0, 1, 2]]))
                for vib_index3 in range(6, n_modes)])
    h23C = ((v(vib_index1) ** (-1))
            * sum([
                sum([((baan(rot_index1, rot_index4, vib_index1)
                       * baan(rot_index3, rot_index5, vib_index2)
                       - baan(rot_index1, rot_index4, vib_index2)
                       * baan(rot_index3, rot_index5, vib_index1))
                      * ee(rot_index4, rot_index5, rot_index2))
                     for rot_index5 in [0, 1, 2]])
                for rot_index4 in [0, 1, 2]]))
    h23D = (3 / 2) * h23sumA(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2)
    h23E = (-4
            * B[rot_index1]
            * B[rot_index2]
            * B[rot_index3]
            * h23sumB(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2))
    h23F = (4
            * B[rot_index1]
            * B[rot_index2]
            * B[rot_index3]
            * h23sumC(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2))
    h23G = 4 * h23sumD(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2)
    h23H = 2 * h23sumE(rot_index1, rot_index2, rot_index3, vib_index1, vib_index2)
    value = ((h23A + h23B + h23C + h23D + h23E + h23F + h23G + h23H)
             * (1 / 2)
             * (((vib_quanta1 + 1) * vib_quanta2) ** (1 / 2)))
    # return [h23A,h23B,h23C,h23D,h23E,h23F,h23G,h23H,value]
    return value


# -------------------------------------------------- #
# Identifying candidate pairs of vibrational states. #
# -------------------------------------------------- #
def custom_sort(finalist):
    value = sorted([sorted([finalist[0][0:2], finalist[0][2:4]]), sorted([finalist[1][0:2], finalist[1][2:4]])])
    return value


max_k = int(freq_cutoff // min(fun_vibs_formed))
candidates = []
for kk in range(0, max_k + 1):
    for ll in range(0, max_k + 1):
        for mode_k in range(6, n_modes):
            for mode_l in range(6, n_modes):
                if mode_k != mode_l:
                    state_freq = kk * w(mode_k) + ll * w(mode_l)
                    if 0 < state_freq < freq_cutoff:
                        candidates.append([[kk, mode_k, ll, mode_l], state_freq])

finalists = []
for x in range(0, len(candidates)):
    for y in range(0, len(candidates)):
        ex = candidates[x][1]
        ey = candidates[y][1]
        energy_diff = abs(ex - ey)
        if energy_diff < max_freq_diff:
            if candidates[x][0][1] == candidates[y][0][1]:
                if candidates[x][0][3] == candidates[y][0][3]:
                    if candidates[x][0][0] + 1 == candidates[y][0][0]:
                        if candidates[x][0][2] - 1 == candidates[y][0][2]:
                            finalists.append([candidates[x][0], candidates[y][0], [ex, ey]])

finalists2 = sorted([[x, custom_sort(x)] for x in finalists], key=lambda element: element[1])
uniques = [finalists2[0][0]]
for x in range(1, len(finalists2)):
    if finalists2[x - 1][1] != finalists2[x][1]:
        uniques.append(finalists2[x][0])

uniques.sort(key=lambda element: sorted(element[2]))

column_labels = ','.join(['Mode 1', 'Mode 2', 'Mode 1 (cm^-1)', 'Mode 2 (cm^-1)',
                          'A', 'B', 'C',  # Ground deperturbed sum of alphas.
                          'A', 'B', 'C',  # Mode 1 deperturbed sum of alphas.
                          'A', 'B', 'C',  # Mode 2 deperturbed sum of alphas.
                          'Ga(A)', 'Gb(A)', 'Gc(A)',  # With higher order corrections
                          'Ga(S)', 'Gb(S)', 'Gc(S)',
                          'Fbc(A)', 'Fac(A)', 'Fab(A)',
                          'Fbc(S)', 'Fac(S)', 'Fab(S)',
                          'GaJ(A)', 'GbJ(A)', '-i*GcJ(A)',
                          'GaK(A)', 'GbK(A)', '-i*GcK(A)',
                          'GaJ(S)', 'GbJ(S)', '-i*GcJ(S)',
                          'GaK(S)', 'G+K(S)', '-i*G-K(S)',
                          'A', 'B', 'C',  # Ground perturbed sum of alphas.
                          'A', 'B', 'C',  # Mode 1 perturbed sum of alphas.
                          'A', 'B', 'C',  # Mode 2 perturbed sum of alphas.
                          'Ga', 'Gb', 'Gc'  # Without higher order corrections.
                          ])

header = """
Using CFOUR output - assuming IIIr representation
- frequency-order naming of vibrational modes
Maximum vibrational frequency = {freq_cutoff:.2f} cm^-1
Maximum frequency difference = {max_freq_diff:.2f} cm^-1

# Sign of Ga relative to Fbc depends on order of vibration specified in fit file:
#   Mode 1 <-> x0011 and Mode 2 <-> x0022
#   If vice versa then multiply only Ga Gb and Gc terms by -1
# Sign of Fbc Gaj Gak etc. are independent of order of vibration
# (S) => symmetric reduction
# (A) => asymmetric reduction
# Using convention of Bv=Be-alphasum

Units in MHz unless specified otherwise
""".format(freq_cutoff=freq_cutoff, max_freq_diff=max_freq_diff) + ',' * 4 + \
         'Ground Deperturbed Sum of Alphas' + ',' * 3 + \
         'Mode 1 Deperturbed Sum of Alphas' + ',' * 3 + \
         'Mode 2 Deperturbed Sum of Alphas' + ',' * 3 + \
         'With All Higher Order Corrections' + ',' * 24 + \
         'Ground Perturbed Sum of Alphas' + ',' * 3 + \
         'Mode 1 Perturbed Sum of Alphas' + ',' * 3 + \
         'Mode 2 Perturbed Sum of Alphas' + ',' * 3 + \
         'Without Higher Order Corrections' + ',' * 3 + \
         '\n' + column_labels + '\n'

rows = []
#for unique in uniques:
for unique in [uniques[0]]:
    [[kk, mode_k, ll, mode_l], [kk2, mode_k2, ll2, mode_l2], [freq_1, freq_2]] = unique

    if kk > 0 and ll == 0:
        mode_1 = '{}m{}'.format(kk, mode_k + 1)
    elif kk == 0 and ll > 0:
        mode_1 = '{}m{}'.format(ll, mode_l + 1)
    else:
        mode_1 = '{}m{}+{}m{}'.format(kk, mode_k+1, ll, mode_l+1)

    if kk2 > 0 and ll2 == 0:
        mode_2 = '{}m{}'.format(kk2, mode_k2 + 1)
    elif kk2 == 0 and ll2 > 0:
        mode_2 = '{}m{}'.format(ll2, mode_l2 + 1)
    else:
        mode_2 = '{}m{}+{}m{}'.format(kk2, mode_k2, ll2, mode_l2)

    mode_1_freq = '{:.2f}'.format(freq_1)
    mode_2_freq = '{:.2f}'.format(freq_2)
    quanta_ground = [0 for i in range(0, n_modes)]
    quanta1 = [0 for i in range(0, n_modes)]
    quanta1[mode_k] = kk
    quanta1[mode_l] = ll
    quanta2 = [0 for i in range(0, n_modes)]
    quanta2[mode_k2] = kk2
    quanta2[mode_l2] = ll2

    # ------------------ #
    # Alphas #
    # ------------------ #

    deperturbed_alphas = []
    perturbed_alphas = []
    for quanta, modes in [[quanta_ground, [mode_k, mode_l]],
                          [quanta1, [mode_k, mode_l]],
                          [quanta2, [mode_k2, mode_l2]]]:
        for axis in [0, 1, 2]:
            deperturbed_alphas.append(deperturbed_alpha_sum(axis, quanta, [modes], []))
            perturbed_alphas.append(alpha_sum(axis, quanta))

    # ------------------ #
    # Defining the "t"'s #
    # ------------------ #
    t = {}
    for rot_index1 in [0, 1, 2]:
        t[rot_index1] = h21(rot_index1, mode_k, mode_l, kk, ll)
        for rot_index2 in [0, 1, 2]:
            t[rot_index1, rot_index2] = h22(rot_index1, rot_index2, mode_k, mode_l, kk, ll)
            for rot_index3 in [0, 1, 2]:
                t[rot_index1, rot_index2, rot_index3] = h23(rot_index1, rot_index2, rot_index3, mode_k, mode_l, kk, ll)

    # -------------------------------- #
    # Automated Derivation definitions #
    # -------------------------------- #
    # STANDARD FORM
    h = {}
    h_imaginary = {}
    # from export_txt function. Lower order terms include corrections from higher order terms.
    h[3, 0, 0] = t[0, 0, 0] / 2
    h[2, 1, 0] = t[0, 0, 1] / 2 + t[0, 1, 0] / 2 + t[1, 0, 0] / 2
    h[2, 0, 1] = t[0, 0, 2] / 2 + t[0, 2, 0] / 2 + t[2, 0, 0] / 2
    h[1, 2, 0] = t[0, 1, 1] / 2 + t[1, 0, 1] / 2 + t[1, 1, 0] / 2
    h[1, 1, 1] = t[0, 1, 2] / 2 + t[0, 2, 1] / 2 + t[1, 0, 2] / 2 + t[1, 2, 0] / 2 + t[2, 0, 1] / 2 + t[2, 1, 0] / 2
    h[1, 0, 2] = t[0, 2, 2] / 2 + t[2, 0, 2] / 2 + t[2, 2, 0] / 2
    h[0, 3, 0] = t[1, 1, 1] / 2
    h[0, 2, 1] = t[1, 1, 2] / 2 + t[1, 2, 1] / 2 + t[2, 1, 1] / 2
    h[0, 1, 2] = t[1, 2, 2] / 2 + t[2, 1, 2] / 2 + t[2, 2, 1] / 2
    h[0, 0, 3] = t[2, 2, 2] / 2

    h[2, 0, 0] = (t[0, 0] / 2)
    h_imaginary[2, 0, 0] = (- t[0, 1, 2] / 4 + t[0, 2, 1] / 4
                            - t[1, 0, 2] / 4 - t[1, 2, 0] / 4
                            + t[2, 0, 1] / 4 + t[2, 1, 0] / 4)
    h[1, 1, 0] = (t[0, 1] / 2 + t[1, 0] / 2)
    h_imaginary[1, 1, 0] = (t[0, 0, 2] / 2 - t[1, 1, 2] / 2
                            - t[2, 0, 0] / 2 + t[2, 1, 1] / 2)
    h[1, 0, 1] = (t[0, 2] / 2 + t[2, 0] / 2)
    h_imaginary[1, 0, 1] = (- t[0, 0, 1] / 2 + t[1, 0, 0] / 2
                            - t[1, 2, 2] / 2 + t[2, 2, 1] / 2)
    h[0, 2, 0] = (t[1, 1] / 2)
    h_imaginary[0, 2, 0] = (+ t[0, 1, 2] / 4 + t[0, 2, 1] / 4
                            + t[1, 0, 2] / 4 - t[1, 2, 0] / 4
                            - t[2, 0, 1] / 4 - t[2, 1, 0] / 4)
    h[0, 1, 1] = (t[1, 2] / 2 + t[2, 1] / 2)
    h_imaginary[0, 1, 1] = (- t[0, 1, 1] / 2 + t[0, 2, 2] / 2
                            + t[1, 1, 0] / 2 - t[2, 2, 0] / 2)
    h[0, 0, 2] = (t[2, 2] / 2)
    h_imaginary[0, 0, 2] = (- t[0, 1, 2] / 4 - t[0, 2, 1] / 4
                            + t[1, 0, 2] / 4 + t[1, 2, 0] / 4
                            - t[2, 0, 1] / 4 + t[2, 1, 0] / 4)

    h[1, 0, 0] = (t[0] / 2 - t[1, 0, 1] / 4 - t[2, 0, 2] / 4)
    h_imaginary[1, 0, 0] = - t[1, 2] / 4 + t[2, 1] / 4
    h[0, 1, 0] = (t[1] / 2 - t[0, 1, 0] / 4 - t[2, 1, 2] / 4)
    h_imaginary[0, 1, 0] = + t[0, 2] / 4 - t[2, 0] / 4
    h[0, 0, 1] = (t[2] / 2 - t[0, 2, 0] / 4 - t[1, 2, 1] / 4)
    h_imaginary[0, 0, 1] = - t[0, 1] / 4 + t[1, 0] / 4

    bad_imaginary = []
    for key, value in h_imaginary.items():
        if h[key] > 0:
            if abs(value/h[key]) > 0.001:
                bad_imaginary.append(key)
    if len(bad_imaginary) > 0:
        print('Warning: The imaginary portions of the following Standard Form terms are significantly larger than 0:\n')
        for bad in bad_imaginary:
            print('h[{key}] (real) = {real}  vs  h[{key}] (imag) = {imag}\n'.format(
                key=bad, real=h[bad], imag=h_imaginary[bad]))

    # ASYMMETRIC REDUCTION
    ca = {}
    ca[2, 1, 0, 0] = (2 * (B[0] - B[2]) * h[0, 2, 1] / (B[0] + B[1] - 2 * B[2])
                      + 2 * (B[1] - B[2]) * h[2, 0, 1] / (B[0] + B[1] - 2 * B[2]))
    ca[2, 0, 1, 0] = h[3, 0, 0]
    ca[2, 0, 0, 1] = -h[0, 3, 0]  # imaginary
    ca[0, 3, 0, 0] = (2 * (-B[0] + B[2]) * h[0, 2, 1] / (B[0] + B[1] - 2 * B[2])
                      + 2 * (-B[1] + B[2]) * h[2, 0, 1] / (B[0] + B[1] - 2 * B[2])
                      + 2 * h[0, 0, 3])
    ca[0, 2, 1, 0] = (h[1, 0, 2] / 2
                      + (B[0] - B[2]) * h[1, 2, 0] / (2 * (B[0] - B[1]))
                      + (-B[0] + B[1] / 2 + B[2] / 2) * h[3, 0, 0] / (B[0] - B[1]))
    ca[0, 2, 0, 1] = (-h[0, 1, 2] / 2
                      + (B[1] - B[2]) * h[2, 1, 0] / (2 * (B[0] - B[1]))
                      + (B[0] / 2 - B[1] + B[2] / 2) * h[0, 3, 0] / (B[0] - B[1]))  # imaginary
    ca_third_order = [
        ca[2, 1, 0, 0],  # G_a^J
        ca[2, 0, 1, 0],  # G_b^J
        ca[2, 0, 0, 1],  # G_c^J, imaginary
        ca[0, 3, 0, 0],  # G_a^K
        ca[0, 2, 1, 0],  # G_b^K
        ca[0, 2, 0, 1]  # G_c^K, imaginary
    ]

    ca[2, 0, 0, 0] = h[0, 2, 0] + h[2, 0, 0]
    ca[0, 2, 0, 0] = 2 * h[0, 0, 2] - h[0, 2, 0] - h[2, 0, 0]
    ca[0, 0, 2, 0] = -h[0, 2, 0] / 2 + h[2, 0, 0] / 2
    ca_second_order = [
        ca[2, 0, 0, 0],  # F_bc
        ca[0, 2, 0, 0],  # F_ac
        ca[0, 0, 2, 0]  # F_ab
    ]

    ca_noHigherOrderCorrections = {key: value for key, value in ca.items()}
    ca[0, 1, 0, 0] = ((-B[0] + B[1]) * h[0, 2, 1] / (2 * (B[0] + B[1] - 2 * B[2]))
                      + (B[0] - B[1]) * h[2, 0, 1] / (2 * (B[0] + B[1] - 2 * B[2]))
                      + 2 * h[0, 0, 1])
    ca_noHigherOrderCorrections[0, 1, 0, 0] = 2 * h[0, 0, 1]
    ca[0, 0, 1, 0] = (h[1, 0, 0]
                      + (-B[1] + B[2]) * h[1, 2, 0] / (4 * (B[0] - B[1]))
                      + (B[1] - B[2]) * h[3, 0, 0] / (4 * (B[0] - B[1])))
    ca_noHigherOrderCorrections[0, 0, 1, 0] = h[1, 0, 0]
    ca[0, 0, 0, 1] = ((-B[0] + B[2]) * h[2, 1, 0] / (4 * (B[0] - B[1]))
                      - h[0, 1, 0]
                      + (B[0] - B[2]) * h[0, 3, 0] / (4 * (B[0] - B[1])))  # imaginary
    ca_noHigherOrderCorrections[0, 0, 0, 1] = - h[0, 1, 0]  # imaginary
    ca_first_order = [
        ca[0, 1, 0, 0],  # G_a
        ca[0, 0, 1, 0],  # G_b
        ca[0, 0, 0, 1]  # G_c, imaginary
    ]
    ca_first_order_noHigherOrderCorrections = [
        ca_noHigherOrderCorrections[0, 1, 0, 0],  # G_a
        ca_noHigherOrderCorrections[0, 0, 1, 0],  # G_b
        ca_noHigherOrderCorrections[0, 0, 0, 1]  # G_c, imaginary
    ]

    # SYMMETRIC REDUCTION
    cs = {}
    cs[2, 1, 0, 0] = (2 * (B[0] - B[2]) * h[0, 2, 1] / (B[0] + B[1] - 2 * B[2])
                      + 2 * (B[1] - B[2]) * h[2, 0, 1] / (B[0] + B[1] - 2 * B[2]))
    cs[2, 0, 1, 0] = ((-B[0] + B[1]) * h[1, 0, 2] / (-5 * B[0] + B[1] + 4 * B[2])
                      + (-B[0] + B[2]) * h[1, 2, 0] / (-5 * B[0] + B[1] + 4 * B[2])
                      + 3 * (-B[0] + B[2]) * h[3, 0, 0] / (-5 * B[0] + B[1] + 4 * B[2]))
    cs[2, 0, 0, 1] = ((-B[0] + B[1]) * h[0, 1, 2] / (B[0] - 5 * B[1] + 4 * B[2])
                      + 3 * (B[1] - B[2]) * h[0, 3, 0] / (B[0] - 5 * B[1] + 4 * B[2])
                      + (B[1] - B[2]) * h[2, 1, 0] / (B[0] - 5 * B[1] + 4 * B[2]))  # imaginary
    cs[0, 3, 0, 0] = (2 * (-B[0] + B[2]) * h[0, 2, 1] / (B[0] + B[1] - 2 * B[2])
                      + 2 * (-B[1] + B[2]) * h[2, 0, 1] / (B[0] + B[1] - 2 * B[2])
                      + 2 * h[0, 0, 3])
    cs[0, 0, 3, 0] = ((B[0] - B[1]) * h[1, 0, 2] / (-5 * B[0] + B[1] + 4 * B[2])
                      + (B[0] - B[2]) * h[1, 2, 0] / (-5 * B[0] + B[1] + 4 * B[2])
                      + (-2 * B[0] + B[1] + B[2]) * h[3, 0, 0] / (-5 * B[0] + B[1] + 4 * B[2]))
    cs[0, 0, 0, 3] = ((-B[0] + B[1]) * h[0, 1, 2] / (B[0] - 5 * B[1] + 4 * B[2])
                      + (B[1] - B[2]) * h[2, 1, 0] / (B[0] - 5 * B[1] + 4 * B[2])
                      + (B[0] - 2 * B[1] + B[2]) * h[0, 3, 0] / (B[0] - 5 * B[1] + 4 * B[2]))  # imaginary
    cs_third_order = [
        cs[2, 1, 0, 0],  # G_a^J
        cs[2, 0, 1, 0],  # G_b^J
        cs[2, 0, 0, 1],  # G_c^J, imaginary
        cs[0, 3, 0, 0],  # G_a^K
        cs[0, 0, 3, 0],  # G_+^K
        cs[0, 0, 0, 3]  # G_-^K, imaginary
    ]

    cs[2, 0, 0, 0] = h[0, 2, 0] + h[2, 0, 0]
    cs[0, 2, 0, 0] = 2 * h[0, 0, 2] - h[0, 2, 0] - h[2, 0, 0]
    cs[0, 0, 2, 0] = -h[0, 2, 0] / 2 + h[2, 0, 0] / 2
    cs_second_order = [
        cs[2, 0, 0, 0],  # F_bc
        cs[0, 2, 0, 0],  # F_ac
        cs[0, 0, 2, 0]  # F_ab
    ]

    cs_noHigherOrderCorrections = {key: value for key, value in cs.items()}
    cs[0, 1, 0, 0] = ((-B[0] + B[1]) * h[0, 2, 1] / (2 * (B[0] + B[1] - 2 * B[2]))
                      + (B[0] - B[1]) * h[2, 0, 1] / (2 * (B[0] + B[1] - 2 * B[2]))
                      + 2 * h[0, 0, 1])
    cs_noHigherOrderCorrections[0, 1, 0, 0] = 2 * h[0, 0, 1]
    cs[0, 0, 1, 0] = (h[1, 0, 0]
                      + (-2 * B[0] + B[1] + B[2]) * h[1, 2, 0] / (4 * (-5 * B[0] + B[1] + 4 * B[2]))
                      + (-B[0] - B[1] + 2 * B[2]) * h[1, 0, 2] / (2 * (-5 * B[0] + B[1] + 4 * B[2]))
                      + (4 * B[0] + B[1] - 5 * B[2]) * h[3, 0, 0] / (4 * (-5 * B[0] + B[1] + 4 * B[2])))
    cs_noHigherOrderCorrections[0, 0, 1, 0] = h[1, 0, 0]
    cs[0, 0, 0, 1] = ((-B[0] - 4 * B[1] + 5 * B[2]) * h[0, 3, 0] / (4 * (B[0] - 5 * B[1] + 4 * B[2]))
                      + (-B[0] + 2 * B[1] - B[2]) * h[2, 1, 0] / (4 * (B[0] - 5 * B[1] + 4 * B[2]))
                      - h[0, 1, 0]
                      + (B[0] + B[1] - 2 * B[2]) * h[0, 1, 2] / (2 * (B[0] - 5 * B[1] + 4 * B[2])))  # imaginary
    cs_noHigherOrderCorrections[0, 0, 0, 1] = - h[0, 1, 0]  # imaginary
    cs_first_order = [
        cs[0, 1, 0, 0],  # G_a
        cs[0, 0, 1, 0],  # G_b
        cs[0, 0, 0, 1]  # G_c, imaginary
    ]
    cs_first_order_noHigherOrderCorrections = [
        cs_noHigherOrderCorrections[0, 1, 0, 0],  # G_a
        cs_noHigherOrderCorrections[0, 0, 1, 0],  # G_b
        cs_noHigherOrderCorrections[0, 0, 0, 1]  # G_c, imaginary
    ]

    # ------------------- #
    # Summarizing results #
    # ------------------- #
    data = [
        *deperturbed_alphas,
        *ca_first_order,
        *cs_first_order,
        *ca_second_order,
        *cs_second_order,
        *ca_third_order,
        *cs_third_order,
        *perturbed_alphas,
        *ca_first_order_noHigherOrderCorrections
    ]

    data = [hertz / 1000000 for hertz in data]  # converting to MHz
    data_str = ['{:.6e}'.format(i) for i in data]
    row = ','.join([mode_1, mode_2, mode_1_freq, mode_2_freq])+','+','.join(data_str)
    rows.append(row)


# ------------ #
# Final output #
# ------------ #
all_rows_str = '\n'.join(rows)
write_string = header + all_rows_str + '\n'
filename = out_file_name.split('.')[0]
outfile = 'coupling2021-of-'+filename+'.csv'
with open(outfile, 'w') as raw_out:
    raw_out.write(write_string)

