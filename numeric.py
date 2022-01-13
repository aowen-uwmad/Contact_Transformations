from main import *
from pathlib import Path
from computational_data_parser import anharmonic_data_import
import numpy as np

# term_definitions = {
#     H(2, 0): (Term([qop(v(0)), qop(v(1))], [], A(2, 0)(v(0), v(1)), [v(0), v(1)], [])
#               + Term([pop(v(0)), pop(v(1))], [], A(2, 0)(v(0), v(1)), [v(0), v(1)], [])),
#     H(3, 0): Term([qop(v(0)), qop(v(1)), qop(v(2))], [], A(3, 0)(v(0), v(1), v(2)), [v(0), v(1), v(2)], []),
#     H(4, 0): (Term([qop(v(0)), qop(v(1)), qop(v(2)), qop(v(3))],
#                    [],
#                    A(4, 0)(v(0), v(1), v(2), v(3)),
#                    [v(0), v(1), v(2), v(3)],
#                    [])
#               + Term([qop(v(0)), pop(v(1)), qop(v(2)), pop(v(3))],
#                      [],
#                      B(4, 0)(v(0), v(1), v(2), v(3)),
#                      [v(0), v(1), v(2), v(3)],
#                      [])
#               ),
#     H(2, 1): Term([qop(v(0)), pop(v(1))], [jop(r(0))], A(2, 1)(v(0), v(1), r(0)), [v(0), v(1)], [r(0)]),
#     H(3, 1): Term([qop(v(0)), pop(v(1)), qop(v(2))],
#                   [jop(r(0))],
#                   A(3, 1)(v(0), v(1), v(2), r(0)),
#                   [v(0), v(1), v(2)],
#                   [r(0)]),
#     H(0, 2): Term([], [jop(r(0)), jop(r(1))], A(0, 2)(r(0), r(1)), [], [r(0), r(1)]),
#     H(1, 2): Term([qop(v(0))], [jop(r(0)), jop(r(1))], A(1, 2)(v(0), r(0), r(1)), [v(0)], [r(0), r(1)]),
#     H(2, 2): Term([qop(v(0)), qop(v(1))],
#                   [jop(r(0)), jop(r(1))],
#                   A(2, 2)(v(0), v(1), r(0), r(1)),
#                   [v(0), v(1)],
#                   [r(0), r(1)]),
#     H(3, 2): Term([qop(v(0)), qop(v(1)), qop(v(2))],
#                   [jop(r(0)), jop(r(1))],
#                   A(3, 2)(v(0), v(1), v(2), r(0), r(1)),
#                   [v(0), v(1), v(2)],
#                   [r(0), r(1)])
# }
#
# ht22, ht22_t, ht22_s = find_target_and_definitions(2, 2, term_definitions,)

h_not_bar, speed_of_light, h_bar, n_modes, w, v, \
zeta, b0, B0, k3, aD, i0, omega, baan, baann, ee, tau = anharmonic_data_import(
    [Path('benzonitrile_cfour_outputs/anharm.out'),
        Path('benzonitrile_cfour_outputs/didQ'),
        Path('benzonitrile_cfour_outputs/corioliszeta'),
        Path('benzonitrile_cfour_outputs/cubic')],
    program='CFOUR')


