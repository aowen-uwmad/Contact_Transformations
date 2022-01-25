from main import *
from pathlib import Path
from computational_data_parser import anharmonic_data_import
import numpy as np
from functools import lru_cache
from tqdm import tqdm


resonance_threshold = 25  # wavenumbers (CM-1)


(h_not_bar, speed_of_light, h_bar, n_modes, w0, v0,
 zeta, b0, B0, k3, aD, i0, omega, baan, baann, ee_f, tau) = anharmonic_data_import(
    [Path('benzonitrile_cfour_outputs/anharm.out'),
        Path('benzonitrile_cfour_outputs/didQ'),
        Path('benzonitrile_cfour_outputs/corioliszeta'),
        Path('benzonitrile_cfour_outputs/cubic')],
    program='CFOUR')


resonance_threshold_hz = resonance_threshold * 100 * speed_of_light  # Hertz (Hz)


cache_max_size = 2 ** 16


@lru_cache(maxsize=cache_max_size)
def caan(ra, rb, va):
    value = (-1)*baan(ra, rb, va)/v0(va)
    return value


@lru_cache(maxsize=cache_max_size)
def A20(va, vb):
    if va == vb:
        return v0(va) / 2
    else:
        return 0


@lru_cache(maxsize=cache_max_size)
def A30(va, vb, vc):
    return k3(va, vb, vc)/6


@lru_cache(maxsize=cache_max_size)
def A21(va, vb, ra):
    value = (-2) * (v0(vb) / v0(va)) ** (1 / 2) * B0(ra) * zeta(ra, va, vb)
    return value


@lru_cache(maxsize=cache_max_size)
def A31(va, vb, vc, ra):
    value = (v0(vb) / (v0(va) * v0(vc))) ** (1 / 2) * sum([
        (v0(va) ** (3 / 2) * caan(ra, rb, va) * zeta(rb, vc, vb) + v0(vc) ** (3 / 2) * caan(ra, rb, vc) * zeta(
            rb, va, vb))
        for rb in [0, 1, 2]])
    return value


@lru_cache(maxsize=cache_max_size)
def A02(ra, rb):
    if ra == rb:
        return B0(ra)
    else:
        return 0


@lru_cache(maxsize=cache_max_size)
def A12(va, ra, rb):
    return (-1) * v0(va) * caan(ra, rb, va)


@lru_cache(maxsize=cache_max_size)
def A22(va, vb, ra, rb):
    value = (3 / 8) * v0(va) * v0(vb) * sum([
        ((caan(ra, rc, va) * caan(rc, rb, vb) + caan(ra, rc, vb) * caan(rc, rb, va)) / B0(rc))
        for rc in [0, 1, 2]])
    return value


@lru_cache(maxsize=cache_max_size)
def DC(*args):
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




