from main import *
from pathlib import Path
from computational_data_parser import anharmonic_data_import
import numpy as np
from sympy import lambdify, LeviCivita, KroneckerDelta, Function


resonance_threshold = 25  # wavenumbers (CM-1)


(h_not_bar, speed_of_light, h_bar, n_modes, w0, v0,
 zeta, b0, B0, k3, aD, i0, omega, baan, baann, ee_f, tau) = anharmonic_data_import(
    [Path('benzonitrile_cfour_outputs/anharm.out'),
        Path('benzonitrile_cfour_outputs/didQ'),
        Path('benzonitrile_cfour_outputs/corioliszeta'),
        Path('benzonitrile_cfour_outputs/cubic')],
    program='CFOUR')


resonance_threshold_hz = resonance_threshold * 100 * speed_of_light  # Hertz (Hz)

old_term_definitions = {
    H(2, 0): (Term([qop(v(0)), qop(v(1))], [], A(2, 0)(v(0), v(1)), [v(0), v(1)], [])
              + Term([pop(v(0)), pop(v(1))], [], A(2, 0)(v(0), v(1)), [v(0), v(1)], [])),
    H(3, 0): Term([qop(v(0)), qop(v(1)), qop(v(2))], [], A(3, 0)(v(0), v(1), v(2)), [v(0), v(1), v(2)], []),
    H(4, 0): (Term([qop(v(0)), qop(v(1)), qop(v(2)), qop(v(3))],
                   [],
                   A(4, 0)(v(0), v(1), v(2), v(3)),
                   [v(0), v(1), v(2), v(3)],
                   [])
              + Term([qop(v(0)), pop(v(1)), qop(v(2)), pop(v(3))],
                     [],
                     B(4, 0)(v(0), v(1), v(2), v(3)),
                     [v(0), v(1), v(2), v(3)],
                     [])
              ),
    H(2, 1): Term([qop(v(0)), pop(v(1))], [jop(r(0))], A(2, 1)(v(0), v(1), r(0)), [v(0), v(1)], [r(0)]),
    H(3, 1): Term([qop(v(0)), pop(v(1)), qop(v(2))],
                  [jop(r(0))],
                  A(3, 1)(v(0), v(1), v(2), r(0)),
                  [v(0), v(1), v(2)],
                  [r(0)]),
    H(0, 2): Term([], [jop(r(0)), jop(r(1))], A(0, 2)(r(0), r(1)), [], [r(0), r(1)]),
    H(1, 2): Term([qop(v(0))], [jop(r(0)), jop(r(1))], A(1, 2)(v(0), r(0), r(1)), [v(0)], [r(0), r(1)]),
    H(2, 2): Term([qop(v(0)), qop(v(1))],
                  [jop(r(0)), jop(r(1))],
                  A(2, 2)(v(0), v(1), r(0), r(1)),
                  [v(0), v(1)],
                  [r(0), r(1)]),
    H(3, 2): Term([qop(v(0)), qop(v(1)), qop(v(2))],
                  [jop(r(0)), jop(r(1))],
                  A(3, 2)(v(0), v(1), v(2), r(0), r(1)),
                  [v(0), v(1), v(2)],
                  [r(0), r(1)])
}


def caan(ra, rb, va):
    value = (-1)*baan(ra, rb, va)/v0(va)
    return value


class A20(Function):
    n_args = 2

    @classmethod
    def eval(cls, va, vb):
        if va.is_Integer and vb.is_Integer:
            if va == vb:
                return v0(va) / 2
            else:
                return 0


class A30(Function):
    n_args = 3

    @classmethod
    def eval(cls, va, vb, vc):
        if all([va.is_Integer, vb.is_Integer, vc.is_Integer]):
            return k3(va, vb, vc)/6


class A21(Function):
    n_args = 3

    @classmethod
    def eval(cls, va, vb, ra):
        if all([va.is_Integer, vb.is_Integer, ra.is_Integer]):
            value = (-2) * (v0(vb) / v0(va)) ** (1 / 2) * B0(ra) * zeta(ra, va, vb)
            return value


class A31(Function):
    n_args = 4

    @classmethod
    def eval(cls, va, vb, vc, ra):
        if all(x.is_Integer for x in [va, vb, vc, ra]):
            value = (v0(vb) / (v0(va) * v0(vc))) ** (1 / 2) * sum([
                (v0(va) ** (3 / 2) * caan(ra, rb, va) * zeta(rb, vc, vb) + v0(vc) ** (3 / 2) * caan(ra, rb, vc) * zeta(
                    rb, va, vb))
                for rb in [0, 1, 2]])
            return value


class A02(Function):
    n_args = 2

    @classmethod
    def eval(cls, ra, rb):
        if ra.is_Integer and rb.is_Integer:
            if ra == rb:
                return B0(ra)
            else:
                return 0


class A12(Function):
    n_args = 3

    @classmethod
    def eval(cls, va, ra, rb):
        if all(x.is_Integer for x in [va, ra, rb]):
            return (-1) * v0(va) * caan(ra, rb, va)


class A22(Function):
    n_args = 4

    @classmethod
    def eval(cls, va, vb, ra, rb):
        if all(x.is_Integer for x in [va, vb, ra, rb]):
            value = (3 / 8) * v0(va) * v0(vb) * sum([
                ((caan(ra, rc, va) * caan(rc, rb, vb) + caan(ra, rc, vb) * caan(rc, rb, va)) / B0(rc))
                for rc in [0, 1, 2]])
            return value


class DC(Function):
    n_args = [i for i in range(2, 20)]  # CHANGE THIS for vib orders greater than ten.

    @classmethod
    def eval(cls, *args):
        args_list = [*args]
        counter = 0
        pairs = []
        while counter+1 < len(args_list):
            new_pair = [args_list[counter], args_list[counter+1]]
            counter += 2
            pairs.append(new_pair)
        denominator_sum = 0
        for pair in pairs:
            sign = pair[0]
            index = list(pair[1].args)[0]
            denominator_sum += sign * v0(index)
        if denominator_sum <= resonance_threshold_hz:
            return 0
        else:
            return 1 / denominator_sum


term_definitions = {
    H(2, 0): (Term([qop(v(0)), qop(v(1))], [], A20(v(0), v(1)), [v(0), v(1)], [])
              + Term([pop(v(0)), pop(v(1))], [], A20(v(0), v(1)), [v(0), v(1)], [])),
    H(3, 0): Term([qop(v(0)), qop(v(1)), qop(v(2))], [], A30(v(0), v(1), v(2)), [v(0), v(1), v(2)], []),
    # H(4, 0): (Term([qop(v(0)), qop(v(1)), qop(v(2)), qop(v(3))],
    #                [],
    #                A40(v(0), v(1), v(2), v(3)),
    #                [v(0), v(1), v(2), v(3)],
    #                [])
    #           + Term([qop(v(0)), pop(v(1)), qop(v(2)), pop(v(3))],
    #                  [],
    #                  B40(v(0), v(1), v(2), v(3)),
    #                  [v(0), v(1), v(2), v(3)],
    #                  [])
    #           ),
    H(2, 1): Term([qop(v(0)), pop(v(1))], [jop(r(0))], A21(v(0), v(1), r(0)), [v(0), v(1)], [r(0)]),
    H(3, 1): Term([qop(v(0)), pop(v(1)), qop(v(2))],
                  [jop(r(0))],
                  A31(v(0), v(1), v(2), r(0)),
                  [v(0), v(1), v(2)],
                  [r(0)]),
    H(0, 2): Term([], [jop(r(0)), jop(r(1))], A02(r(0), r(1)), [], [r(0), r(1)]),
    H(1, 2): Term([qop(v(0))], [jop(r(0)), jop(r(1))], A12(v(0), r(0), r(1)), [v(0)], [r(0), r(1)]),
    H(2, 2): Term([qop(v(0)), qop(v(1))],
                  [jop(r(0)), jop(r(1))],
                  A22(v(0), v(1), r(0), r(1)),
                  [v(0), v(1)],
                  [r(0), r(1)]),
    # H(3, 2): Term([qop(v(0)), qop(v(1)), qop(v(2))],
    #               [jop(r(0)), jop(r(1))],
    #               A32(v(0), v(1), v(2), r(0), r(1)),
    #               [v(0), v(1), v(2)],
    #               [r(0), r(1)])
}


ht22, ht22_t, ht22_s = find_target_and_definitions(2, 2, term_definitions, full_simplify=True)

ht22_A, ht22_B, ht22_C, ht22_D = ht22


def ht22_A_function(va: int, vb: int, ra: int, rb: int):
    vib_op_indices = [x.args[0] for x in ht22_A.vib_op]  # n = 2
    rot_op_indices = [x.args[0] for x in ht22_A.rot_op]  # n = 2
    vib_indices = ht22_A.vib_indices  # n = 5
    rot_indices = ht22_A.rot_indices  # n = 5
    vib_sum_indices = [x for x in vib_indices if x not in vib_op_indices]  # n = 3
    rot_sum_indices = [x for x in rot_indices if x not in rot_op_indices]  # n = 3
    new_vib_indices = []
    counter = 0
    for index in vib_indices:
        if index in vib_op_indices:
            new_vib_indices.append([va, vb][counter])
            counter += 1
        else:
            new_vib_indices.append(index)
    new_rot_indices = []
    counter = 0
    for index in rot_indices:
        if index in rot_op_indices:
            new_rot_indices.append([ra, rb][counter])
            counter += 1
        else:
            new_rot_indices.append(index)
    new_ht22_A = ht22_A
    new_ht22_A = new_ht22_A.changeIndices(new_vib_indices, new_rot_indices)
    coefficient_expression = new_ht22_A.coefficient.replace(D, DC)
    coefficient_expression = coefficient_expression.replace(ee, LeviCivita)
    new_coefficient = 0
    for vs0 in range(n_modes):
        for vs1 in range(n_modes):
            for vs2 in range(n_modes):
                for rs0 in range(3):
                    for rs1 in range(3):
                        for rs2 in range(3):
                            sub_rules = {
                                vib_sum_indices[0]: vs0,
                                vib_sum_indices[1]: vs1,
                                vib_sum_indices[2]: vs2,
                                rot_sum_indices[0]: vs0,
                                rot_sum_indices[1]: vs1,
                                rot_sum_indices[2]: vs2
                            }
                            new_coefficient += coefficient_expression.subs(sub_rules, simultaneous=True)
    new_ht22_A.coefficient = new_coefficient
    return new_ht22_A



