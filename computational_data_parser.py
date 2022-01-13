from pathlib import Path
import numpy as np


def gaussian16_anharmonic_import(gaussian_output_file: Path,
                                 h_not_bar=6.626*10**(-34),
                                 speed_of_light=2.998*10**8):
    if not gaussian_output_file.is_file():
        raise FileNotFoundError('Could not find {}'.format(gaussian_output_file))

    # TODO: create anharmonic import function. See previous import script as template.from

    # return h_not_bar, speed_of_light, h_bar, n_modes, w, v, zeta, b0, B0, k3, aD, i0, omega, baan, baann, ee, tau


def cfour_vpt2_import(anharm_out_file: Path,
                      didQ_file: Path,
                      coriolis_zeta_file: Path,
                      cubic_file: Path,
                      h_not_bar=6.626*10**(-34),
                      speed_of_light=2.998*10**8):
    files_present = True
    for in_file in [anharm_out_file, didQ_file, coriolis_zeta_file, cubic_file]:
        if not in_file.is_file():
            files_present = False
            print('Could not find {}'.format(in_file))
    if not files_present:
        raise FileNotFoundError

    h_bar = h_not_bar / (2 * np.pi)

    with open(anharm_out_file, 'r') as anharm_in_raw:
        anharm_in = anharm_in_raw.read()
    with open(didQ_file, 'r') as didQ_in_raw:
        didQ_in = didQ_in_raw.read().split('\n')
    with open(coriolis_zeta_file, 'r') as zetas_in_raw:
        zetas_in = zetas_in_raw.read().split('Coriolis')
    with open(cubic_file, 'r') as cubic_in_raw:
        cubic_in = cubic_in_raw.read().split('\n')

    # parsing anharm.out
    try:
        rot_in = anharm_in.split('Be, B0 AND B-B0 SHIFTS FOR SINGLY EXCITED VIBRATIONAL STATES (CM-1)')[1].split(
            'Vibrationally averaged dipole moment')[0].split('\n')[4].split()[1:]
        rot_formed = [float(x) for x in rot_in]
        vib_in = anharm_in.split('Rotationally projected vibrational frequencies')[-1].split(
            'Finite temperature thermodynamic corrections')[0].split('\n')[7:-3]
        # first 6 "vibrations" in CFOUR are the rotation and translations with frequencies of 0. Using dummy energies
        # for those.
        fun_vibs_formed = [1000000.0] * 6
        fun_vibs_formed.extend([float(vibration.split()[1]) for vibration in vib_in])
        n_modes = len(fun_vibs_formed)
    except IndexError:
        raise ValueError('Unable to parse the anharm.out file located at {}'.format(anharm_out_file))

    # parsing corioliszeta
    try:
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
    except IndexError:
        raise ValueError('Unable to parse the cubic file located at {}'.format(coriolis_zeta_file))

    # parsing cubic
    try:
        cubic = [[[int(x.split()[0]) - 1, int(x.split()[1]) - 1, int(x.split()[2]) - 1], float(x.split()[3])]
                 for x in cubic_in if x != '']
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
    except IndexError:
        raise ValueError('Unable to parse the cubic file located at {}'.format(cubic_file))

    # parsing didq
    try:
        didQ = [[[int(didQ_in[x].split()[0]) - 1, int(didQ_in[x].split()[1]) - 1, int(didQ_in[x].split()[2]) - 1],
                 float(didQ_in[x].split()[3])] for x in range(0, len(didQ_in)) if didQ_in[x] != '']
        # Making full inertia derivatives data array
        didQ_out = [[[0 for x in range(0, 3)] for y in range(0, 3)] for z in range(0, n_modes)]
        for x in range(0, len(didQ)):
            didQ_out[didQ[x][0][2]][didQ[x][0][0]][didQ[x][0][1]] = didQ[x][1]
        inertia_der_formed = didQ_out
    except IndexError:
        raise ValueError('Unable to parse the didQ file located at {}'.format(didQ_file))

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
    def b0(i):
        return rot_formed[i]  # wavenumbers

    def B0(i):
        return rot_formed[i] * 100 * speed_of_light  # hertz

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
        return h_not_bar / (8 * (np.pi ** 2) * B0(rot_index))  # kg*m**2

    def omega(vib_index1, vib_index2, vib_index3):
        value = (v(vib_index1) + v(vib_index2) + v(vib_index3)) * \
                (-v(vib_index1) + v(vib_index2) + v(vib_index3)) * \
                (v(vib_index1) - v(vib_index2) + v(vib_index3)) * \
                (v(vib_index1) + v(vib_index2) - v(vib_index3))
        return value  # hertz**4

    def baan(rot_index1, rot_index2, vib_index1):
        # hertz
        value = ((-1) * (h_bar ** 3) / (2 * (h_not_bar ** (3 / 2)))) * (
                aD(rot_index1, rot_index2, vib_index1) / (i0(rot_index1) * i0(rot_index2) * (v(vib_index1) ** (1 / 2))))
        return value

    def baann(rot_index1, rot_index2, vib_index1, vib_index2):
        # hertz
        value = (3 / 8) * sum([((1 / B0(rotIn3)) * (
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

    return h_not_bar, speed_of_light, h_bar, n_modes, w, v, zeta, b0, B0, k3, aD, i0, omega, baan, baann, ee, tau


def anharmonic_data_import(files, program: str, h_not_bar=6.626*10**(-34), speed_of_light=2.998*10**8):
    if program.lower() in ['cfour', 'cfour18', 'c', 'c4', 'cccc']:
        if isinstance(files, (list, tuple)) and len(files) == 4:
            list_of_path_objects = []
            for file in files:
                if isinstance(file, Path):
                    list_of_path_objects.append(file)
                elif isinstance(file, str):
                    list_of_path_objects.append(Path(file))
                else:
                    raise TypeError('Unrecognized input: {}'.format(file))
            anharm_out_file, didQ_file, coriolis_zeta_file, cubic_file = list_of_path_objects
            return cfour_vpt2_import(anharm_out_file,
                                     didQ_file,
                                     coriolis_zeta_file,
                                     cubic_file,
                                     h_not_bar=h_not_bar,
                                     speed_of_light=speed_of_light)
        else:
            raise TypeError('Improper list of CFOUR input files.')
    elif program.lower() in ['gaussian', 'gaussian16', 'g16']:
        if isinstance(files, (list, tuple)) and len(files) == 1:
            gaussian_file = files[0]
        elif isinstance(files, (Path, str)):
            gaussian_file = files
        else:
            raise TypeError('Improper input of Gaussian file.')
        if not isinstance(gaussian_file, Path):
            if isinstance(gaussian_file, str):
                gaussian_file = Path(gaussian_file)
            else:
                raise TypeError('Unrecognized input: {}'.format(gaussian_file))
        return gaussian16_anharmonic_import(gaussian_file, h_not_bar=h_not_bar, speed_of_light=speed_of_light)
