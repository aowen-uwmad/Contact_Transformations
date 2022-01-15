from pathlib import Path
from computational_data_parser import anharmonic_data_import
import numpy as np
from sympy import LeviCivita, Function, I, preorder_traversal, Symbol, symbols, Rational
from math import factorial


def printDictionary(dictionary: dict):
    for key, value in dictionary.items():
        print('{}:  {}'.format(key, value))


v0, v1, v2, v3, v4, v5, v6, v7, v8, v9 = symbols('v0, v1, v2, v3, v4, v5, v6, v7, v8, v9')
vibration_indices = [v0, v1, v2, v3, v4, v5, v6, v7, v8, v9]

r0, r1, r2, r3, r4, r5, r6, r7, r8, r9 = symbols('r0, r1, r2, r3, r4, r5, r6, r7, r8, r9')
rotation_indices = [r0, r1, r2, r3, r4, r5, r6, r7, r8, r9]

jop = Function('jop')  # j rotational operator
pop = Function('pop')  # p (momentum) vibrational operator
qop = Function('qop')  # q (position) vibrational operator
lop = Function('lop')  # L (ladder) vibrational operator
sigma = Function('sigma')  # sign (+/-)
omega = Function('omega')  # harmonic frequency
VC = Function('VC')  # Vibrational Commutator
RC = Function('RC')  # Rotational Commutator
nm = Symbol('nm')  # Number of Vibrational Modes


class H(Function):
    n_args = 2

    @classmethod
    def eval(cls, m, n):
        if m.is_Integer and n.is_Integer:
            return Symbol('H{}{}'.format(m, n))


class S(Function):
    n_args = 2

    @classmethod
    def eval(cls, m, n):
        if m.is_Integer and n.is_Integer:
            return Symbol('S{}{}'.format(m, n))


class GenericTerm:
    def __init__(self, coefficient, m: int, n: int, term_type: str):
        if term_type in ['H', 'S']:
            self.type = term_type
        else:
            raise AttributeError('Unrecognized type {}.'.format(term_type))
        self.coefficient = coefficient
        if m < 0 or n < 0 or (m+n-2) < 0:
            raise AttributeError('Invalid orders of operators {} and {}'.format(m, n))
        else:
            self.vib_order = m
            self.rot_order = n
            if m == 0 and n == 2:
                self.order = 1
            else:
                self.order = m + n - 2
        if self.type == 'H':
            self.symbol = H(self.vib_order, self.rot_order)
        elif self.type == 'S':
            self.symbol = S(self.vib_order, self.rot_order)

    def __repr__(self):
        if self.coefficient == 1:
            return str(self.symbol)
        else:
            return '({})*{}'.format(self.coefficient, self.symbol)

    def __add__(self, other):
        if isinstance(other, GenericTerm):
            if self.combinesWith(other):
                self.coefficient = self.coefficient + other.coefficient
                return self
            else:
                return GenericExpression([self, other])
        elif isinstance(other, GenericExpression):
            items_list = other.items
            items_list.append(self)
            return GenericExpression(items_list)
        elif isinstance(other, GenericCommutator):
            return GenericExpression([self, other])
        else:
            raise NotImplementedError

    def __sub__(self, other):
        return self + (other*(-1))

    def __mul__(self, other):
        if isinstance(other, (GenericTerm, GenericExpression, GenericCommutator)):
            raise NotImplementedError
        else:
            try:
                new = self
                new.coefficient = new.coefficient*other
                return new
            except:
                raise TypeError

    def __rmul__(self, other):
        if isinstance(other, (GenericTerm, GenericExpression, GenericCommutator)):
            raise NotImplementedError
        else:
            try:
                new = self
                new.coefficient = new.coefficient*other
                return new
            except:
                raise TypeError

    def __eq__(self, other):
        if self.combinesWith(other) and self.coefficient == other.coefficient:
            return True
        else:
            return False

    def sympy(self):
        return self.coefficient * self.symbol

    def combinesWith(self, other):
        if isinstance(other, GenericTerm):
            if self.type == other.type and self.vib_order == other.vib_order and self.rot_order == other.rot_order:
                return True
            else:
                return False
        else:
            return False

    def vibCommutator(self, other):
        if isinstance(other, (GenericTerm, GenericCommutator)):
            if self.vib_order == 0 or other.vib_order == 0:
                return 0
            else:
                return GenericCommutator(1, self, other, 'V')
        elif isinstance(other, GenericExpression):
            item_list = [self.vibCommutator(item) for item in other.items]
            return GenericExpression(item_list)
        else:
            raise NotImplementedError

    def rotCommutator(self, other):
        if isinstance(other, (GenericTerm, GenericCommutator)):
            if self.rot_order == 0 or other.rot_order == 0:
                return 0
            else:
                return GenericCommutator(1, self, other, 'R')
        elif isinstance(other, GenericExpression):
            item_list = [self.rotCommutator(item) for item in other.items]
            return GenericExpression(item_list)
        else:
            raise NotImplementedError

    def select(self, m: int, n: int):
        if self.vib_order == m and self.rot_order == n:
            return self
        else:
            return 0

    def toEquation(self, definitions: dict):
        if self.symbol not in definitions.keys():
            raise KeyError('{} is not defined in the provided dictionary.'.format(self.symbol))
        return definitions[self.symbol] * self.coefficient


class GenericExpression:
    def __init__(self, items_list):
        self.items = []
        for item in items_list:
            if isinstance(item, GenericExpression):
                for sub_item in item.items:
                    self.items.append(sub_item)
            elif isinstance(item, (GenericTerm, GenericCommutator)):
                if item.order >= 0 and not ([item.vib_order, item.rot_order] in [[0, 0], [1, 0], [0, 1], [1, 1]]):
                    self.items.append(item)
            elif item == 0:
                pass
            else:
                raise ValueError('Unrecognized object encountered: {}'.format(item))

        unique_items = []
        for item in self.items:
            if not any(item.combinesWith(unique) for unique in unique_items):
                unique_items.append(item)

        if len(unique_items) != len(self.items):
            combined_items = []
            for unique in unique_items:
                combines_with_items = []
                for item in self.items:
                    if item.combinesWith(unique):
                        combines_with_items.append(item)
                item_sum = combines_with_items[0]
                for item in combines_with_items[1:]:
                    item_sum += item
                combined_items.append(item_sum)
            self.items = combined_items

        final_items = []
        for item in self.items:
            if item.coefficient != 0:
                final_items.append(item)

    def __repr__(self):
        return str(self.items)

    def __add__(self, other):
        if isinstance(other, (GenericTerm, GenericCommutator)):
            return GenericExpression([*self.items, other])
        elif isinstance(other, GenericExpression):
            return GenericExpression([*self.items, *other.items])
        elif other == 0:
            return self
        else:
            raise TypeError('Cannot add {} to GenericExpression object'.format(other))

    def __sub__(self, other):
        return self + (other*(-1))

    def __mul__(self, other):
        if isinstance(other, (GenericTerm, GenericExpression, GenericCommutator)):
            raise NotImplementedError
        else:
            try:
                value = GenericExpression([item * other for item in self.items])
                return value
            except:
                raise TypeError('Cannot multiply GenericExpression object by {}'.format(other))

    def __rmul__(self, other):
        return self*other

    def __getitem__(self, item):
        return self.items[item]

    def __eq__(self, other):
        if isinstance(other, GenericExpression):
            raise NotImplementedError
        else:
            return False

    def sympy(self):
        n_items = len(self.items)
        item_sum = self.items[0].sympy()
        if n_items == 0:
            return 0
        elif n_items == 1:
            return item_sum
        else:
            for item in self.items[1:]:
                item_sum = item_sum + item.sympy()
            return item_sum

    def select(self, m: int, n: int):
        return GenericExpression([item.select(m, n) for item in self.items])

    def toEquation(self, definitions: dict):
        return Expression([item.toEquation(definitions) for item in self.items])


class GenericCommutator:
    def __init__(self, coefficient, term1, term2, com_type: str):
        if com_type == 'V' or com_type == 'R':
            self.type = com_type
        else:
            raise TypeError
        if isinstance(term1, (GenericTerm, GenericCommutator)) and isinstance(term2, (GenericTerm, GenericCommutator)):
            self.coefficient = coefficient*term1.coefficient*term2.coefficient
            self.term1 = term1
            self.term2 = term2
            self.term1.coefficient = 1
            self.term2.coefficient = 1
            if self.type == 'V':
                self.vib_order = term1.vib_order + term2.vib_order - 2
                self.rot_order = term1.rot_order + term2.rot_order
            elif self.type == 'R':
                self.vib_order = term1.vib_order + term2.vib_order
                self.rot_order = term1.rot_order + term2.rot_order - 1
            self.order = self.vib_order + self.rot_order - 2
        elif isinstance(term1, GenericExpression) or isinstance(term2, GenericExpression):
            raise TypeError('Use .vibCommutator(GenericExpression) or .rotCommutator(GenericExpression) to evaluate '
                            'the commutator of an expression')
        else:
            raise NotImplementedError

    def __repr__(self):
        if self.coefficient == 1:
            return '[{},{}]{}'.format(self.term1, self.term2, self.type)
        else:
            return '({})*[{},{}]{}'.format(self.coefficient, self.term1, self.term2, self.type)

    def __add__(self, other):
        if isinstance(other, GenericCommutator):
            if self.combinesWith(other):
                new = self
                new.coefficient = new.coefficient + other.coefficient
                return new
            else:
                return GenericExpression([self, other])
        elif isinstance(other, (GenericTerm, GenericExpression)):
            return GenericExpression([self, other])
        else:
            raise TypeError

    def __sub__(self, other):
        return self + (other * (-1))

    def __mul__(self, other):
        if isinstance(other, (GenericTerm, GenericCommutator, GenericExpression)):
            raise NotImplementedError
        else:
            try:
                new = self
                new.coefficient = new.coefficient * other
                return new
            except:
                raise TypeError

    def __rmul__(self, other):
        if isinstance(other, (GenericTerm, GenericCommutator, GenericExpression)):
            raise NotImplementedError
        else:
            try:
                new = self
                new.coefficient = new.coefficient * other
                return new
            except:
                raise TypeError

    def __eq__(self, other):
        if self.combinesWith(other) and self.coefficient == other.coefficient:
            return True
        else:
            return False

    def sympy(self):
        if self.type == 'V':
            return self.coefficient * VC(self.term1.sympy(), self.term2.sympy())
        elif self.type == 'R':
            return self.coefficient * RC(self.term1.sympy(), self.term2.sympy())

    def combinesWith(self, other):
        if isinstance(other, GenericCommutator):
            if self.term1 == other.term1 and self.term2 == other.term2 and self.type == other.type:
                return True
            else:
                return False
        else:
            return False

    def select(self, m: int, n: int):
        if self.vib_order == m and self.rot_order == n:
            return self
        else:
            return 0

    def toEquation(self, definitions: dict):
        new_term1 = self.term1.toEquation(definitions)
        new_term2 = self.term2.toEquation(definitions)
        if self.type == 'V':
            return new_term1.vibCommutator(new_term2)
        elif self.type == 'R':
            return new_term1.rotCommutator(new_term2)
        else:
            raise AttributeError


def H_Group(i: int):
    if i == 0:
        value = GenericTerm(1, 2, 0, 'H')
    elif i == 1:
        value = (GenericTerm(1, 0, 2, 'H')
                 + GenericTerm(1, 3, 0, 'H')
                 + GenericTerm(1, 2, 1, 'H')
                 + GenericTerm(1, 1, 2, 'H'))
    elif i > 1:
        value = GenericTerm(1, i + 2, 0, 'H') + GenericTerm(1, i + 1, 1, 'H') + GenericTerm(1, i, 2, 'H')
    else:
        raise ValueError('Order of {} is not allowed.'.format(i))
    return value


def targetExpression(i: int, j: int, do_print=False):
    def transformedH(k, level):
        if level == 0:
            return H_Group(k)
        elif level > 0:
            if k == 0:
                return H_Group(0)
            elif k > 0:
                value = transformedH(k, level-1)
                s_term = transforms[level-1]
                for l in range(0, k):
                    coefficient = I**(k - l)*Rational(1, (factorial(k - l)))
                    commutator = GenericExpression([
                        s_term.vibCommutator(transformedH(l, level-1)),
                        s_term.rotCommutator(transformedH(l, level-1))
                    ])
                    counter = 1
                    while counter < k - l:
                        commutator = GenericExpression([
                            s_term.vibCommutator(commutator),
                            s_term.rotCommutator(commutator)
                        ])
                        counter += 1
                    value = GenericExpression([commutator * coefficient, value])
                return value
            else:
                raise AttributeError
        else:
            raise AttributeError

    max_order = i + j - 2
    transforms = []
    if i > 0:
        m_range = i + 1
        n_range = j + 1
    elif i == 0:
        m_range = 2
        n_range = int(j / 2) + 1
    else:
        raise ValueError
    for m in range(1, m_range):
        n_start = max([3-m, 0])
        for n in range(n_start, n_range):
            transforms.append(GenericTerm(1, m, n, 'S'))
    n_transforms = len(transforms)
    up_to_max_order = GenericExpression([transformedH(k, n_transforms) for k in range(0, max_order + 1)])
    select_orders = up_to_max_order.select(i, j)
    defining_equations = []
    for level in range(0, n_transforms):
        transform = transforms[level]
        a = transform.vib_order
        b = transform.rot_order
        term_to_block_diagonalize = GenericExpression(
            [transformedH(k, level) for k in range(0, max_order + 1)]
        ).select(a, b)
        defining_equations.append([transform, term_to_block_diagonalize])

    if do_print:
        print(select_orders.sympy())
        printDefiningEquations(defining_equations)
    return select_orders, defining_equations


def printDefiningEquations(defining_equations):
    for item in defining_equations:
        s_term = item[0]
        other_part = item[1]
        vib_order = s_term.vib_order
        rot_order = s_term.rot_order
        print('H_t({},{}) = ({}) + I*VC({}, H20)'.format(vib_order, rot_order, other_part.sympy(), s_term))


class Term:

