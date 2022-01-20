from pathlib import Path
from computational_data_parser import anharmonic_data_import
import numpy as np
from sympy import LeviCivita, Function, I, preorder_traversal, Symbol, symbols, Rational, summation, Mul
from sympy.physics.quantum.operator import Operator
from math import factorial
from itertools import product


def printDictionary(dictionary: dict):
    for key, value in dictionary.items():
        print('{}:  {}'.format(key, value))


class v(Function):
    n_args = 1

    @classmethod
    def eval(cls, index):
        if index.is_Integer:
            return Symbol('v{}'.format(index))


class r(Function):
    n_args = 1

    @classmethod
    def eval(cls, index):
        if index.is_Integer:
            return Symbol('r{}'.format(index))


class ee(LeviCivita):
    pass


vibration_indices = [v(i) for i in range(10)]
rotation_indices = [r(i) for i in range(10)]

jop = Function('jop', is_commutative=False)  # j rotational operator
pop = Function('pop', is_commutative=False)  # p (momentum) vibrational operator
qop = Function('qop', is_commutative=False)  # q (position) vibrational operator
lop = Function('lop', is_commutative=False)  # L (ladder) vibrational operator
sigma = Function('sigma')  # sign (+/-)
omega = Function('omega')  # harmonic frequency
VC = Function('VC')  # Vibrational Commutator
RC = Function('RC')  # Rotational Commutator
D = Function('D')  # Denominator
nm = Symbol('nm')  # Number of Vibrational Modes


class H(Function):
    n_args = 2

    @classmethod
    def eval(cls, m, n):
        if m.is_Integer and n.is_Integer:
            return Symbol('H{}{}'.format(m, n))


class h(Function):
    n_args = (2, 3)

    @classmethod
    def eval(cls, *args):
        if all(x.is_Integer for x in args):
            if len(args) == 2:
                return Function('h{}{}'.format(*args))
            elif len(args) == 3:
                return Function('h{}{}_{}'.format(*args))


class S(Function):
    n_args = 2

    @classmethod
    def eval(cls, m, n):
        if m.is_Integer and n.is_Integer:
            return Symbol('S{}{}'.format(m, n))


class s(Function):
    n_args = (2, 3)

    @classmethod
    def eval(cls, *args):
        if all(x.is_Integer for x in args):
            if len(args) == 2:
                return Function('s{}{}'.format(*args))
            elif len(args) == 3:
                return Function('s{}{}_{}'.format(*args))


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
    def __init__(self, vib_op_list: list, rot_op_list: list, coefficient, vib_indices: list, rot_indices: list):
        self._n_vib_op = len(vib_op_list)
        self._n_rot_op = len(rot_op_list)
        self._vib_op = vib_op_list
        self._rot_op = rot_op_list
        self._coefficient = coefficient
        self._vib_indices = vib_indices
        self._rot_indices = rot_indices

    @property
    def vib_op(self):
        return self._vib_op

    @property
    def rot_op(self):
        return self._rot_op

    @property
    def coefficient(self):
        return self._coefficient

    @coefficient.setter
    def coefficient(self, value):
        self._coefficient = value

    @property
    def n_vib_op(self):
        return self._n_vib_op

    @property
    def n_rot_op(self):
        return self._n_rot_op

    @property
    def vib_indices(self):
        return self._vib_indices

    @property
    def rot_indices(self):
        return self._rot_indices

    def __repr__(self):
        operators = []
        for vib_op in self.vib_op:
            operators.append(str(vib_op))
        for rot_op in self.rot_op:
            operators.append(str(rot_op))
        return "Sum(({})*{})".format(self.coefficient, "*".join(operators))

    def __add__(self, other):
        if isinstance(other, Term):
            if self.willCombineWith(other):
                return self.combineWith(other)
            else:
                return Expression([self, other])
        elif isinstance(other, Expression):
            return Expression([self, *other.items])
        elif other == 0:
            return self
        else:
            raise TypeError('Cannot add {} to Term object.'.format(other))

    def __sub__(self, other):
        return self + (other*(-1))

    def __mul__(self, other):
        if isinstance(other, Expression):
            return Expression([self*term for term in other.items])
        elif isinstance(other, Term):
            n_left_vib_indices = len(self.vib_indices)
            n_left_rot_indices = len(self.rot_indices)
            n_right_vib_indices = len(other.vib_indices)
            n_right_rot_indices = len(other.rot_indices)

            left_vib_indices = [vibration_indices[i] for i in range(0, n_left_vib_indices)]
            right_vib_indices = [vibration_indices[i] for i in range(n_left_vib_indices, n_left_vib_indices+n_right_vib_indices)]
            left_rot_indices = [rotation_indices[i] for i in range(0, n_left_rot_indices)]
            right_rot_indices = [rotation_indices[i] for i in range(n_left_rot_indices, n_left_rot_indices+n_right_rot_indices)]

            new_left = self.changeIndices(left_vib_indices, left_rot_indices)
            new_right = other.changeIndices(right_vib_indices, right_rot_indices)

            final_vib_indices = left_vib_indices + right_vib_indices
            final_rot_indices = left_rot_indices + right_rot_indices
            final_coefficient = new_left.coefficient * new_right.coefficient
            final_vib_op = new_left.vib_op + new_right.vib_op
            final_rot_op = new_left.rot_op + new_right.rot_op

            return Term(final_vib_op, final_rot_op, final_coefficient, final_vib_indices, final_rot_indices)
        else:
            try:
                new_coefficient = self.coefficient*other
                return Term(self.vib_op, self.rot_op, new_coefficient, self.vib_indices, self.rot_indices)
            except:
                raise TypeError('Cannot multiply coefficient by {}.'.format(other))

    def __rmul__(self, other):
        raise NotImplementedError

    def __eq__(self, other):
        if isinstance(other, Term):
            if all([self.vib_op == other.vib_op,
                    self.rot_op == other.rot_op,
                    self.vib_indices == other.vib_indices,
                    self.rot_indices == other.rot_indices,
                    self.coefficient == other.coefficient]):
                return True
            else:
                return False
        else:
            return False

    def __len__(self):
        return 1

    def willCombineWith(self, other):
        if isinstance(other, Term):
            if self.n_vib_op != other.n_vib_op:
                return False
            elif self.n_rot_op != other.n_rot_op:
                return False
            for i in range(0, self.n_vib_op):
                if self.vib_op[i].func != other.vib_op[i].func:
                    return False
            for i in range(0, self.n_rot_op):
                if self.rot_op[i].func != other.rot_op[i].func:
                    return False
            return True
        else:
            return False

    def changeVibIndices(self, new_vib_indices_list):
        if len(new_vib_indices_list) != len(self.vib_indices):
            raise ValueError('Unequal length of vib indices lists.')
        substitution_rules = {}
        for i in range(0, len(self.vib_indices)):
            substitution_rules[self.vib_indices[i]] = new_vib_indices_list[i]
        new_vib_op = [x.subs(substitution_rules, simultaneous=True) for x in self.vib_op]
        new_coefficient = self.coefficient.subs(substitution_rules, simultaneous=True)
        new_vib_indices = []
        for x in new_vib_indices_list:
            if x not in new_vib_indices:
                new_vib_indices.append(x)
        return Term(new_vib_op, self.rot_op, new_coefficient, new_vib_indices, self.rot_indices)

    def changeRotIndices(self, new_rot_indices_list):
        if len(new_rot_indices_list) != len(self.rot_indices):
            raise ValueError('Unequal length of rot indices lists.')
        substitution_rules = {}
        for i in range(0, len(self.rot_indices)):
            substitution_rules[self.rot_indices[i]] = new_rot_indices_list[i]
        new_rot_op = [x.subs(substitution_rules, simultaneous=True) for x in self.rot_op]
        new_coefficient = self.coefficient.subs(substitution_rules, simultaneous=True)
        return Term(self.vib_op, new_rot_op, new_coefficient, self.vib_indices, new_rot_indices_list)

    def changeIndices(self, new_vib_indices_list, new_rot_indices_list):
        new_vib_term = self.changeVibIndices(new_vib_indices_list)
        final_term = new_vib_term.changeRotIndices(new_rot_indices_list)
        return final_term

    def combineWith(self, other):
        if self.willCombineWith(other):
            if len(self.vib_indices) <= len(other.vib_indices):
                smaller_vib_term = self
                larger_vib_term = other
            else:
                smaller_vib_term = other
                larger_vib_term = self
            diff_n_vib_indices = len(larger_vib_term.vib_indices) - len(smaller_vib_term.vib_indices)
            combined_vib_indices = larger_vib_term.vib_indices
            combined_vib_op = larger_vib_term.vib_op
            if smaller_vib_term.vib_indices == larger_vib_term.vib_indices:
                new_smaller_vib_term = smaller_vib_term
            else:
                # Changing smaller term to match indices of larger term.  That way the additional indices of the larger
                # term are unaffected.  If same, then doesn't matter.
                vib_substitutions = {}
                for i in range(0, smaller_vib_term.n_vib_op):
                    smaller_vib_index = [x for x in preorder_traversal(smaller_vib_term.vib_op[i])][1]
                    larger_vib_index = [x for x in preorder_traversal(larger_vib_term.vib_op[i])][1]
                    vib_substitutions[smaller_vib_index] = larger_vib_index
                if len(vib_substitutions.keys()) < len(smaller_vib_term.vib_indices):
                    has_rules = [i for i in vib_substitutions.keys()]
                    needs_rules = [i for i in smaller_vib_term.vib_indices if i not in has_rules]
                    used_indices = [value for key, value in vib_substitutions.items()]
                    unused_indices = [i for i in larger_vib_term.vib_indices if i not in used_indices]
                    for i in range(0, len(needs_rules)):
                        vib_substitutions[needs_rules[i]] = unused_indices[i]
                new_vib_indices = [x.subs(vib_substitutions, simultaneous=True) for x in smaller_vib_term.vib_indices]
                new_smaller_vib_term = smaller_vib_term.changeVibIndices(new_vib_indices)
                if diff_n_vib_indices > 0:
                    new_smaller_vib_term.coefficient = (1/nm)**diff_n_vib_indices*new_smaller_vib_term.coefficient
            if len(new_smaller_vib_term.rot_indices) <= len(larger_vib_term.rot_indices):
                smaller_rot_term = new_smaller_vib_term
                larger_rot_term = larger_vib_term
            else:
                smaller_rot_term = larger_vib_term
                larger_rot_term = new_smaller_vib_term
            diff_n_rot_indices = len(larger_rot_term.rot_indices) - len(smaller_rot_term.rot_indices)
            combined_rot_indices = larger_rot_term.rot_indices
            combined_rot_op = larger_rot_term.rot_op
            if smaller_rot_term.rot_indices == larger_rot_term.rot_indices:
                new_smaller_rot_term = smaller_rot_term
            else:
                rot_substitutions = {}
                for i in range(0, smaller_rot_term.n_rot_op):
                    smaller_rot_index = [x for x in preorder_traversal(smaller_rot_term.rot_op[i])][1]
                    larger_rot_index = [x for x in preorder_traversal(larger_rot_term.rot_op[i])][1]
                    rot_substitutions[smaller_rot_index] = larger_rot_index
                if len(rot_substitutions.keys()) < len(smaller_rot_term.rot_indices):
                    has_rules = [i for i in rot_substitutions.keys()]
                    needs_rules = [i for i in smaller_rot_term.rot_indices if i not in has_rules]
                    used_indices = [value for key, value in rot_substitutions.items()]
                    unused_indices = [i for i in larger_rot_term.rot_indices if i not in used_indices]
                    for i in range(0, len(needs_rules)):
                        rot_substitutions[needs_rules[i]] = unused_indices[i]
                new_rot_indices = [x.subs(rot_substitutions, simultaneous=True) for x in smaller_rot_term.rot_indices]
                new_smaller_rot_term = smaller_rot_term.changeRotIndices(new_rot_indices)
                if diff_n_rot_indices > 0:
                    new_smaller_rot_term.coefficient = Rational(1, 3)**diff_n_rot_indices*new_smaller_rot_term.coefficient
            combined_coefficient = larger_rot_term.coefficient + new_smaller_rot_term.coefficient
            combined_term = Term(combined_vib_op,
                                 combined_rot_op,
                                 combined_coefficient,
                                 combined_vib_indices,
                                 combined_rot_indices)
            return combined_term
        else:
            return Expression([self, other])

    def vibCommutator(self, other):
        if self.n_vib_op == 0:
            return 0
        elif isinstance(other, Expression):
            return Expression([self.vibCommutator(term) for term in other.items])
        elif isinstance(other, Term):
            if other.n_vib_op == 0:
                return 0
            else:
                n_left_vib_indices = len(self.vib_indices)
                n_left_rot_indices = len(self.rot_indices)
                n_right_vib_indices = len(other.vib_indices)
                n_right_rot_indices = len(other.rot_indices)

                left_vib_indices = [v(i) for i in range(0, n_left_vib_indices)]
                right_vib_indices = [v(i) for i in range(n_left_vib_indices, n_left_vib_indices+n_right_vib_indices)]
                left_rot_indices = [r(i) for i in range(0, n_left_rot_indices)]
                right_rot_indices = [r(i) for i in range(n_left_rot_indices, n_left_rot_indices+n_right_rot_indices)]

                new_left = self.changeIndices(left_vib_indices, left_rot_indices)
                new_right = other.changeIndices(right_vib_indices, right_rot_indices)

                final_vib_indices = left_vib_indices + right_vib_indices
                final_rot_indices = left_rot_indices + right_rot_indices
                final_coefficient = Rational(1, 2) * new_left.coefficient * new_right.coefficient

                vib_op_list = pure_vibration_commutator(new_left.vib_op, new_right.vib_op)
                rot_op_list = [new_left.rot_op+new_right.rot_op, new_right.rot_op+new_left.rot_op]

                final_terms = []
                if vib_op_list == 0:
                    return 0
                else:
                    for item in vib_op_list:
                        if item == 0:
                            pass
                        else:
                            vib_op = item[0]
                            kronecker_delta_rules = item[1]
                            multiplier = item[2]
                            for rot_op in rot_op_list:
                                new_term = Term(vib_op,
                                                rot_op,
                                                multiplier*final_coefficient,
                                                final_vib_indices,
                                                final_rot_indices)
                                new_vib_indices = [x.subs(kronecker_delta_rules,
                                                          simultaneous=True) for x in final_vib_indices]
                                new_term = new_term.changeVibIndices(new_vib_indices)
                                final_terms.append(new_term)

                return Expression(final_terms)
        else:
            return 0

    def rotCommutator(self, other):
        if self.n_rot_op == 0:
            return 0
        elif isinstance(other, Expression):
            return Expression([self.rotCommutator(term) for term in other.items])
        elif isinstance(other, Term):
            if other.n_rot_op == 0:
                return 0
            else:
                n_left_vib_indices = len(self.vib_indices)
                n_left_rot_indices = len(self.rot_indices)
                n_right_vib_indices = len(other.vib_indices)
                n_right_rot_indices = len(other.rot_indices)

                left_vib_indices = [v(i) for i in range(0, n_left_vib_indices)]
                right_vib_indices = [v(i) for i in range(n_left_vib_indices, n_left_vib_indices+n_right_vib_indices)]
                left_rot_indices = [r(i) for i in range(0, n_left_rot_indices)]
                right_rot_indices = [r(i) for i in range(n_left_rot_indices, n_left_rot_indices+n_right_rot_indices)]

                new_left = self.changeIndices(left_vib_indices, left_rot_indices)
                new_right = other.changeIndices(right_vib_indices, right_rot_indices)

                final_vib_indices = left_vib_indices + right_vib_indices
                final_rot_indices = left_rot_indices + right_rot_indices
                final_coefficient = Rational(1, 2) * new_left.coefficient * new_right.coefficient

                vib_op_list = [new_left.vib_op + new_right.vib_op, new_right.vib_op + new_left.vib_op]
                rot_op_list = pure_rotation_commutator(new_left.rot_op,
                                                       new_right.rot_op,
                                                       new_left.rot_indices,
                                                       new_right.rot_indices)

                final_terms = []
                if rot_op_list == 0:
                    return 0
                else:
                    for item in rot_op_list:
                        if item == 0:
                            pass
                        else:
                            rot_op = item[0]
                            new_indices = item[1]
                            multiplier = item[2]
                            for vib_op in vib_op_list:
                                new_term = Term(vib_op,
                                                rot_op,
                                                multiplier*final_coefficient,
                                                final_vib_indices,
                                                final_rot_indices + list(new_indices))
                                final_terms.append(new_term)
                final_expression = Expression(final_terms)
                for term in final_expression:
                    ee_atoms = []
                    for atom in list(preorder_traversal(term.coefficient)):
                        try:
                            if atom.func == ee:
                                ee_atoms.append(atom)
                        except:
                            pass
                    ee_replace_rules = {}
                    for atom in ee_atoms:
                        index1, index2, index3 = atom.args
                        if index1 != index2 and index1 != index3 and index2 != index3:
                            pass
                        else:
                            ee_replace_rules[atom] = 0
                    new_coefficient = term.coefficient.subs(ee_replace_rules, simultaneous=True)
                    term.coefficient = new_coefficient
                return final_expression
        else:
            return 0

    def toLadder(self):
        new_coefficient = self.coefficient
        ladder_operators = []
        for vib_op in self.vib_op:
            index = list(preorder_traversal(vib_op))[1]
            ladder_operators.append(lop(index, sigma(index)))
            if vib_op.func == qop:
                new_coefficient = new_coefficient*Rational(1, 2)
            elif vib_op.func == pop:
                new_coefficient = new_coefficient*Rational(1, 2)*I*sigma(index)
        return LadderTerm(ladder_operators, self.rot_op, new_coefficient, self.vib_indices, self.rot_indices)

    def printProperties(self):
        print('vib_op (n={}): {}'.format(self.n_vib_op, self.vib_op))
        print('rot_op (n={}): {}'.format(self.n_rot_op, self.rot_op))
        print('vib_indices: {}'.format(self.vib_indices))
        print('rot_indices: {}'.format(self.rot_indices))
        print('coefficient: {}\n'.format(self.coefficient))

    def sympy(self, number_of_vib_modes=nm):
        summation_indices = []
        for vib_index in self.vib_indices:
            if not isinstance(vib_index, int):
                summation_indices.append((vib_index, 0, number_of_vib_modes-1))
        for rot_index in self.rot_indices:
            if not isinstance(rot_index, int):
                summation_indices.append((rot_index, 0, 2))
        operators = []
        for vib_op in self.vib_op:
            operators.append(vib_op)
        for rot_op in self.rot_op:
            operators.append(rot_op)
        op_part = Mul(*operators)
        if len(summation_indices) > 0:
            result = summation(self.coefficient*op_part, *summation_indices)
        else:
            result = self.coefficient*op_part
        return result

    def vibrational_matrix_element(self, state1, state2, basis):
        if not isinstance(state1, VibState):
            raise TypeError
        elif not isinstance(state2, VibState):
            raise TypeError

        operator_types = []
        for operator in self.vib_op:
            if operator.func == qop:
                operator_types.append('q')
            elif operator.func == pop:
                operator_types.append('p')

        possible_indices_ranges = [range(starting_mode, n_modes)]*self.n_vib_op
        possible_indices = list(product(*possible_indices_ranges))

        sympy_result = 0
        for indices_list in possible_indices:
            operators_list = list(zip(operator_types, indices_list))
            result = evaluate_multi_vib_op(operators_list, state1, state2, basis)
            if len(result) > 0:
                for entry in result:
                    changed_term = self.changeVibIndices(indices_list)
                    new_term = Term([],
                                    changed_term.rot_op,
                                    changed_term.coefficient,
                                    [i for i in changed_term.vib_indices if not isinstance(i, int)],
                                    changed_term.rot_indices)
                    sympy_result += entry[0]*new_term.sympy(n_modes)

        return sympy_result


def pure_vibration_commutator(left: list, right: list):
    if len(left) == 0:
        return 0
    elif len(left) == 1:
        if len(right) == 0:
            return 0
        elif len(right) == 1:
            a = left[0]
            b = right[0]
            a_index = list(preorder_traversal(a))[1]
            b_index = list(preorder_traversal(b))[1]
            if a.func == b.func:
                return 0
            elif a.func == qop and b.func == pop:
                kronecker_delta_rules = {a_index: b_index}
                multiplier = I
                return [[[], kronecker_delta_rules, multiplier]]
            elif a.func == pop and b.func == qop:
                kronecker_delta_rules = {a_index: b_index}
                multiplier = (-1)*I
                return [[[], kronecker_delta_rules, multiplier]]
            else:
                raise ValueError
        elif len(right) > 1:
            a = left[0]
            b = right[0]
            c = right[1:]
            first_commutator = pure_vibration_commutator([a], c)
            second_commutator = pure_vibration_commutator([a], [b])
            final_list = []
            if first_commutator == 0:
                pass
            else:
                for item in first_commutator:
                    if item == 0:
                        pass
                    else:
                        vib_op = item[0]
                        kd_rules = item[1]
                        mult = item[2]
                        new_op = [b]
                        for op in vib_op:
                            new_op.append(op)
                        final_list.append([new_op, kd_rules, mult])
            if second_commutator == 0:
                pass
            else:
                for item in second_commutator:
                    if item == 0:
                        pass
                    else:
                        vib_op = item[0]
                        kd_rules = item[1]
                        mult = item[2]
                        for op in c:
                            vib_op.append(op)
                        final_list.append([vib_op, kd_rules, mult])
            if len(final_list) == 0:
                return 0
            else:
                return final_list
        else:
            raise ValueError
    elif len(left) > 1:
        a = left[0]
        b = left[1:]
        c = right
        first_commutator = pure_vibration_commutator(b, c)
        second_commutator = pure_vibration_commutator([a], c)
        final_list = []
        if first_commutator == 0:
            pass
        else:
            for item in first_commutator:
                if item == 0:
                    pass
                else:
                    vib_op = item[0]
                    kd_rules = item[1]
                    mult = item[2]
                    new_op = [a]
                    for op in vib_op:
                        new_op.append(op)
                    final_list.append([new_op, kd_rules, mult])
        if second_commutator == 0:
            pass
        else:
            for item in second_commutator:
                if item == 0:
                    pass
                else:
                    vib_op = item[0]
                    kd_rules = item[1]
                    mult = item[2]
                    for op in b:
                        vib_op.append(op)
                    final_list.append([vib_op, kd_rules, mult])
        if len(final_list) == 0:
            return 0
        else:
            return final_list
    else:
        raise ValueError


def pure_rotation_commutator(left: list, right: list, left_indices: list, right_indices: list):
    if len(left) == 0:
        return 0
    elif len(left) == 1:
        if len(right) == 0:
            return 0
        elif len(right) == 1:
            a = left[0]
            b = right[0]
            a_index = list(preorder_traversal(a))[1]
            b_index = list(preorder_traversal(b))[1]
            if a_index == b_index:
                return 0
            else:
                used_indices = left_indices + right_indices
                new_index = None
                counter = 0
                while new_index is None:
                    if r(counter) in used_indices:
                        counter += 1
                    else:
                        new_index = r(counter)
                multiplier = (-1)*I*ee(a_index, b_index, new_index)
                rot_op = jop(new_index)
                return [[[rot_op], [new_index], multiplier]]
        elif len(right) > 1:
            a = left[0]
            b = right[0]
            c = right[1:]
            first_commutator = pure_rotation_commutator([a], c, left_indices, right_indices)
            second_commutator = pure_rotation_commutator([a], [b], left_indices, right_indices)
            final_list = []
            if first_commutator == 0:
                pass
            else:
                for item in first_commutator:
                    if item == 0:
                        pass
                    else:
                        rot_op = item[0]
                        new_indices = item[1]
                        mult = item[2]
                        new_op = [b]
                        for op in rot_op:
                            new_op.append(op)
                        final_list.append([new_op, new_indices, mult])
            if second_commutator == 0:
                pass
            else:
                for item in second_commutator:
                    if item == 0:
                        pass
                    else:
                        rot_op = item[0]
                        new_indices = item[1]
                        mult = item[2]
                        for op in c:
                            rot_op.append(op)
                        final_list.append([rot_op, new_indices, mult])
            if len(final_list) == 0:
                return 0
            else:
                return final_list
        else:
            raise ValueError
    elif len(left) > 1:
        a = left[0]
        b = left[1:]
        c = right
        first_commutator = pure_rotation_commutator(b, c, left_indices, right_indices)
        second_commutator = pure_rotation_commutator([a], c, left_indices, right_indices)
        final_list = []
        if first_commutator == 0:
            pass
        else:
            for item in first_commutator:
                if item == 0:
                    pass
                else:
                    rot_op = item[0]
                    new_indices = item[1]
                    mult = item[2]
                    new_op = [a]
                    for op in rot_op:
                        new_op.append(op)
                    final_list.append([new_op, new_indices, mult])
        if second_commutator == 0:
            pass
        else:
            for item in second_commutator:
                if item == 0:
                    pass
                else:
                    rot_op = item[0]
                    new_indices = item[1]
                    mult = item[2]
                    for op in b:
                        rot_op.append(op)
                    final_list.append([rot_op, new_indices, mult])
        if len(final_list) == 0:
            return 0
        else:
            return final_list
    else:
        raise ValueError


class Expression:
    def __init__(self, items_list):
        self.items = []
        for item in items_list:
            if isinstance(item, Expression):
                for sub_item in item.items:
                    self.items.append(sub_item)
            elif isinstance(item, Term):
                self.items.append(item)
            elif item == 0:
                pass
            else:
                raise ValueError('Unrecognized object encountered: {}'.format(item))

        combined_terms = []
        for term1 in self.items:
            if any(term1.willCombineWith(term) for term in combined_terms):
                pass
            else:
                combines_with_term1 = []
                for term2 in self.items:
                    if term1.willCombineWith(term2):
                        combines_with_term1.append(term2)
                combined_sum = combines_with_term1[0]
                if len(combines_with_term1) > 1:
                    for term in combines_with_term1[1:]:
                        combined_sum += term
                combined_terms.append(combined_sum)

        final_terms = []
        for term in combined_terms:
            if term.coefficient != 0:
                final_terms.append(term)

        self.items = final_terms

    def __repr__(self):
        return "{}".format([item for item in self.items])

    def __add__(self, other):
        if isinstance(other, Term):
            return Expression([*self.items, other])
        elif isinstance(other, Expression):
            return Expression([*self.items, *other.items])
        elif other == 0:
            return self
        else:
            raise TypeError('Cannot add {} to GenericExpression object'.format(other))

    def __sub__(self, other):
        return self + (other*(-1))

    def __mul__(self, other):
        return Expression([item * other for item in self.items])

    def __rmul__(self, other):
        raise NotImplementedError

    def __getitem__(self, item):
        return self.items[item]

    def __eq__(self, other):
        if isinstance(other, Expression):
            raise NotImplementedError('Need to define expression_sort function first.')
        else:
            return False

    def __len__(self):
        return len([item for item in self])

    def vibCommutator(self, other):
        return Expression([term.vibCommutator(other) for term in self.items])

    def rotCommutator(self, other):
        return Expression([term.rotCommutator(other) for term in self.items])

    def toLadder(self):
        return LadderExpression([term.toLadder() for term in self.items])

    def sympy(self):
        total = 0
        for item in self.items:
            total += item.sympy()
        return total

    def vibrational_matrix_element(self, state1, state2, basis):
        return sum([item.vibrational_matrix_element(state1, state2, basis) for item in self.items])


class LadderTerm:
    def __init__(self, vib_op_list: list, rot_op_list: list, coefficient, vib_indices: list, rot_indices: list):
        self._n_vib_op = len(vib_op_list)
        self._n_rot_op = len(rot_op_list)
        self._vib_indices = vib_indices
        self._rot_indices = rot_indices
        self._vib_op = vib_op_list
        self._rot_op = rot_op_list
        self._coefficient = coefficient

    @property
    def vib_op(self):
        return self._vib_op

    @property
    def rot_op(self):
        return self._rot_op

    @property
    def coefficient(self):
        return self._coefficient

    @coefficient.setter
    def coefficient(self, value):
        self._coefficient = value

    @property
    def n_vib_op(self):
        return self._n_vib_op

    @property
    def n_rot_op(self):
        return self._n_rot_op

    @property
    def vib_indices(self):
        return self._vib_indices

    @property
    def rot_indices(self):
        return self._rot_indices

    def __repr__(self):
        operators = []
        for vib_op in self.vib_op:
            operators.append(str(vib_op))
        for rot_op in self.rot_op:
            operators.append(str(rot_op))
        return "Sum(({})*{})".format(self.coefficient, "*".join(operators))

    def __add__(self, other):
        if isinstance(other, LadderTerm):
            if self.willCombineWith(other):
                return self.combineWith(other)
            else:
                return LadderExpression([self, other])
        elif isinstance(other, LadderExpression):
            return LadderExpression([self, *other.items])
        elif other == 0:
            return self
        else:
            raise TypeError('Cannot add {} to Term object.'.format(other))

    def __sub__(self, other):
        return self + (other*(-1))

    def __mul__(self, other):
        if isinstance(other, LadderExpression):
            return LadderExpression([self * term for term in other.items])
        elif isinstance(other, LadderTerm):
            n_left_vib_indices = len(self.vib_indices)
            n_left_rot_indices = len(self.rot_indices)
            n_right_vib_indices = len(other.vib_indices)
            n_right_rot_indices = len(other.rot_indices)

            left_vib_indices = [v(i) for i in range(0, n_left_vib_indices)]
            right_vib_indices = [v(i) for i in range(n_left_vib_indices, n_left_vib_indices+n_right_vib_indices)]
            left_rot_indices = [r(i) for i in range(0, n_left_rot_indices)]
            right_rot_indices = [r(i) for i in range(n_left_rot_indices, n_left_rot_indices+n_right_rot_indices)]

            new_left = self.changeIndices(left_vib_indices, left_rot_indices)
            new_right = other.changeIndices(right_vib_indices, right_rot_indices)

            final_vib_indices = left_vib_indices + right_vib_indices
            final_rot_indices = left_rot_indices + right_rot_indices
            final_coefficient = new_left.coefficient * new_right.coefficient
            final_vib_op = new_left.vib_op + new_right.vib_op
            final_rot_op = new_left.rot_op + new_right.rot_op

            return LadderTerm(final_vib_op, final_rot_op, final_coefficient, final_vib_indices, final_rot_indices)
        else:
            try:
                new_coefficient = self.coefficient*other
                return LadderTerm(self.vib_op, self.rot_op, new_coefficient, self.vib_indices, self.rot_indices)
            except:
                raise TypeError('Cannot multiply coefficient by {}.'.format(other))

    def __rmul__(self, other):
        raise NotImplementedError

    def __eq__(self, other):
        if isinstance(other, LadderTerm):
            if all([self.vib_op == other.vib_op,
                    self.rot_op == other.rot_op,
                    self.vib_indices == other.vib_indices,
                    self.rot_indices == other.rot_indices,
                    self.coefficient == other.coefficient]):
                return True
            else:
                return False
        else:
            return False

    def __len__(self):
        return 1

    def willCombineWith(self, other):
        if isinstance(other, LadderTerm):
            if self.n_vib_op != other.n_vib_op:
                return False
            for i in range(0, self.n_vib_op):
                if self.vib_op[i].func != other.vib_op[i].func:
                    return False
            if self.n_rot_op != other.n_rot_op:
                return False
            for i in range(0, self.n_rot_op):
                if self.rot_op[i].func != other.rot_op[i].func:
                    return False
            return True
        else:
            return False

    def changeVibIndices(self, new_vib_indices_list):
        if len(new_vib_indices_list) != len(self.vib_indices):
            raise ValueError('Unequal length of vib indices lists.')
        substitution_rules = {}
        for i in range(0, len(self.vib_indices)):
            substitution_rules[self.vib_indices[i]] = new_vib_indices_list[i]
        new_vib_op = [x.subs(substitution_rules, simultaneous=True) for x in self.vib_op]
        new_coefficient = self.coefficient.subs(substitution_rules, simultaneous=True)
        new_vib_indices = []
        for x in new_vib_indices_list:
            if x not in new_vib_indices:
                new_vib_indices.append(x)
        return LadderTerm(new_vib_op, self.rot_op, new_coefficient, new_vib_indices, self.rot_indices)

    def changeRotIndices(self, new_rot_indices_list):
        if len(new_rot_indices_list) != len(self.rot_indices):
            raise ValueError('Unequal length of rot indices lists.')
        substitution_rules = {}
        for i in range(0, len(self.rot_indices)):
            substitution_rules[self.rot_indices[i]] = new_rot_indices_list[i]
        new_rot_op = [x.subs(substitution_rules, simultaneous=True) for x in self.rot_op]
        new_coefficient = self.coefficient.subs(substitution_rules, simultaneous=True)
        return LadderTerm(self.vib_op, new_rot_op, new_coefficient, self.vib_indices, new_rot_indices_list)

    def changeIndices(self, new_vib_indices_list, new_rot_indices_list):
        new_vib_term = self.changeVibIndices(new_vib_indices_list)
        final_term = new_vib_term.changeRotIndices(new_rot_indices_list)
        return final_term

    def combineWith(self, other):
        if self.willCombineWith(other):
            if len(self.vib_indices) <= len(other.vib_indices):
                smaller_vib_term = self
                larger_vib_term = other
            else:
                smaller_vib_term = other
                larger_vib_term = self
            diff_n_vib_indices = len(larger_vib_term.vib_indices) - len(smaller_vib_term.vib_indices)
            # Changing smaller term to match indices of larger term.  That way the additional indices of the larger
            # term are unaffected.  If same, then doesn't matter.
            combined_vib_indices = larger_vib_term.vib_indices
            combined_vib_op = larger_vib_term.vib_op
            if smaller_vib_term.vib_indices == larger_vib_term.vib_indices:
                new_smaller_vib_term = smaller_vib_term
            else:
                vib_substitutions = {}
                for i in range(0, smaller_vib_term.n_vib_op):
                    smaller_vib_index = [x for x in preorder_traversal(smaller_vib_term.vib_op[i])][1]
                    larger_vib_index = [x for x in preorder_traversal(larger_vib_term.vib_op[i])][1]
                    vib_substitutions[smaller_vib_index] = larger_vib_index
                if len(vib_substitutions.keys()) < len(smaller_vib_term.vib_indices):
                    has_rules = [i for i in vib_substitutions.keys()]
                    needs_rules = [i for i in smaller_vib_term.vib_indices if i not in has_rules]
                    used_indices = [value for key, value in vib_substitutions.items()]
                    unused_indices = [i for i in larger_vib_term.vib_indices if i not in used_indices]
                    for i in range(0, len(needs_rules)):
                        vib_substitutions[needs_rules[i]] = unused_indices[i]
                new_vib_indices = [x.subs(vib_substitutions, simultaneous=True) for x in smaller_vib_term.vib_indices]
                new_smaller_vib_term = smaller_vib_term.changeVibIndices(new_vib_indices)
                if diff_n_vib_indices > 0:
                    new_smaller_vib_term.coefficient = (1/nm)**diff_n_vib_indices*new_smaller_vib_term.coefficient
            if len(new_smaller_vib_term.rot_indices) <= len(larger_vib_term.rot_indices):
                smaller_rot_term = new_smaller_vib_term
                larger_rot_term = larger_vib_term
            else:
                smaller_rot_term = larger_vib_term
                larger_rot_term = new_smaller_vib_term
            diff_n_rot_indices = len(larger_rot_term.rot_indices) - len(smaller_rot_term.rot_indices)
            combined_rot_indices = larger_rot_term.rot_indices
            combined_rot_op = larger_rot_term.rot_op
            rot_substitutions = {}
            for i in range(0, smaller_rot_term.n_rot_op):
                smaller_rot_index = [x for x in preorder_traversal(smaller_rot_term.rot_op[i])][1]
                larger_rot_index = [x for x in preorder_traversal(larger_rot_term.rot_op[i])][1]
                rot_substitutions[smaller_rot_index] = larger_rot_index
            new_rot_indices = [x.subs(rot_substitutions, simultaneous=True) for x in smaller_rot_term.rot_indices]
            new_smaller_rot_term = smaller_rot_term.changeRotIndices(new_rot_indices)
            if diff_n_rot_indices > 0:
                new_smaller_rot_term.coefficient = Rational(1, 3)**diff_n_rot_indices*new_smaller_rot_term.coefficient
            combined_coefficient = larger_rot_term.coefficient + new_smaller_rot_term.coefficient
            combined_term = LadderTerm(combined_vib_op,
                                       combined_rot_op,
                                       combined_coefficient,
                                       combined_vib_indices,
                                       combined_rot_indices)
            return combined_term
        else:
            return LadderExpression([self, other])

    def toOperator(self):
        transform_list = [[self.coefficient, self.vib_op]]
        for index in self.vib_indices:
            new_transform_list = []
            for pair in transform_list:
                coefficient = pair[0]
                ops = pair[1]
                new_transform_list.append([coefficient.subs({sigma(index): 1}, simultaneous=True),
                                           [op.subs({sigma(index): 1}, simultaneous=True) for op in ops]])
                new_transform_list.append([coefficient.subs({sigma(index): -1}, simultaneous=True),
                                           [op.subs({sigma(index): -1}, simultaneous=True) for op in ops]])
            transform_list = [x for x in new_transform_list]
        final_terms = []
        for pair in transform_list:
            coefficient = pair[0]
            ladder_ops = pair[1]
            ops_and_mult = []
            for op in ladder_ops:
                index, sign = list(preorder_traversal(op))[1:3]
                q_part = [1, [qop(index)]]
                p_part = [-I * sign, [pop(index)]]
                ops_and_mult.append([q_part, p_part])
            new_pairs = ops_and_mult[0]
            for op_mult in ops_and_mult[1:]:
                new_new_pairs = []
                q_part = op_mult[0]
                p_part = op_mult[1]
                for new_pair in new_pairs:
                    new_new_pairs.append([new_pair[0] * q_part[0], [*new_pair[1], *q_part[1]]])
                    new_new_pairs.append([new_pair[0] * p_part[0], [*new_pair[1], *p_part[1]]])
                new_pairs = [x for x in new_new_pairs]
            for new_pair in new_pairs:
                mult = new_pair[0]
                ops = new_pair[1]
                final_coefficient = coefficient*mult
                new_term = Term(ops, self.rot_op, final_coefficient, self.vib_indices, self.rot_indices)
                final_terms.append(new_term)
        return Expression(final_terms)

    def printProperties(self):
        print('vib_op (n={}): {}'.format(self.n_vib_op, self.vib_op))
        print('rot_op (n={}): {}'.format(self.n_rot_op, self.rot_op))
        print('vib_indices: {}'.format(self.vib_indices))
        print('rot_indices: {}'.format(self.rot_indices))
        print('coefficient: {}\n'.format(self.coefficient))

    def sympy(self):
        summation_indices = []
        for vib_index in self.vib_indices:
            summation_indices.append((vib_index, 0, nm))
        for rot_index in self.rot_indices:
            summation_indices.append((rot_index, 0, 2))
        operators = []
        for vib_op in self.vib_op:
            operators.append(vib_op)
        for rot_op in self.rot_op:
            operators.append(rot_op)
        op_part = Mul(*operators)
        sum_part = summation(self.coefficient * op_part, *summation_indices)
        return sum_part


class LadderExpression:
    def __init__(self, items_list):
        self.items = []
        for item in items_list:
            if isinstance(item, LadderExpression):
                for sub_item in item.items:
                    self.items.append(sub_item)
            elif isinstance(item, LadderTerm):
                self.items.append(item)
            elif item == 0:
                pass
            else:
                raise ValueError('Unrecognized object encountered: {}'.format(item))

        combined_terms = []
        for term1 in self.items:
            if any(term1.willCombineWith(term) for term in combined_terms):
                pass
            else:
                combines_with_term1 = []
                for term2 in self.items:
                    if term1.willCombineWith(term2):
                        combines_with_term1.append(term2)
                combined_sum = combines_with_term1[0]
                if len(combines_with_term1) > 1:
                    for term in combines_with_term1[1:]:
                        combined_sum += term
                combined_terms.append(combined_sum)

        final_terms = []
        for term in combined_terms:
            if term.coefficient != 0:
                final_terms.append(term)

        self.items = final_terms

    def __repr__(self):
        return "{}".format([item for item in self.items])

    def __add__(self, other):
        if isinstance(other, LadderTerm):
            return LadderExpression([*self.items, other])
        elif isinstance(other, LadderExpression):
            return LadderExpression([*self.items, *other.items])
        elif other == 0:
            return self
        else:
            raise TypeError('Cannot add {} to GenericExpression object'.format(other))

    def __sub__(self, other):
        return self + (other*(-1))

    def __mul__(self, other):
        return LadderExpression([item * other for item in self.items])

    def __rmul__(self, other):
        raise NotImplementedError

    def __getitem__(self, item):
        return self.items[item]

    def __eq__(self, other):
        if isinstance(other, LadderExpression):
            raise NotImplementedError('Need to define expression_sort function first.')
        else:
            return False

    def __len__(self):
        return len([item for item in self])

    def toOperator(self):
        return Expression([term.toOperator() for term in self.items])

    def sympy(self):
        total = 0
        for item in self.items:
            total += item.sympy()
        return total


def vib_indices_sorter(indices_list):
    str_list = [str(x) for x in indices_list]
    sorted_list = sorted(str_list)
    if str_list == sorted_list:
        return list(indices_list)
    else:
        new_indices_list = []
        for i in range(0, len(sorted_list)):
            correct_str = sorted_list[i]
            for j in range(0, len(str_list)):
                candidate = str_list[j]
                if candidate == correct_str:
                    new_indices_list.append(indices_list[j])
                    break
        return list(new_indices_list)


def rot_indices_sorter(indices_list):
    str_list = [str(x) for x in indices_list]
    sorted_list = sorted(str_list)
    if str_list == sorted_list:
        return list(indices_list)
    else:
        new_indices_list = []
        for i in range(0, len(sorted_list)):
            correct_str = sorted_list[i]
            for j in range(0, len(str_list)):
                candidate = str_list[j]
                if candidate == correct_str:
                    new_indices_list.append(indices_list[j])
                    break
        return list(new_indices_list)


(h_not_bar, speed_of_light, h_bar, n_modes, w0, v0,
 zeta, b0, B0, k3, aD, i0, capital_omega, baan, baann, ee_f, tau) = anharmonic_data_import(
    [Path('benzonitrile_cfour_outputs/anharm.out'),
        Path('benzonitrile_cfour_outputs/didQ'),
        Path('benzonitrile_cfour_outputs/corioliszeta'),
        Path('benzonitrile_cfour_outputs/cubic')],
    program='CFOUR')


starting_mode = 6  # if CFOUR
# starting_mode = 1  # if Gaussian?

resonance_threshold = 25  # wavenumbers (CM-1)
resonance_threshold_hz = resonance_threshold * 100 * speed_of_light  # Hertz (Hz)

max_state_energy = 1000  # wavenumbers


def describe_normal_modes():
    print("i    wavenumbers")
    for i in range(starting_mode, n_modes):
        print("{}:  {}".format(i, w0(i)))


class VibState:
    def __init__(self, quanta_dict=None):
        state = {i: 0 for i in range(starting_mode, n_modes)}
        if quanta_dict is None:
            self._quanta = state
        elif isinstance(quanta_dict, dict):
            for mode, quanta in quanta_dict.items():
                if quanta < 0:
                    raise ValueError
                state[mode] = quanta
            self._quanta = state
        else:
            raise TypeError

    @property
    def quanta(self):
        return self._quanta

    def __repr__(self):
        quanta_list = ["{}v{}".format(self.quanta[i], i) for i in range(starting_mode, n_modes) if self.quanta[i] > 0]
        if len(quanta_list) == 0:
            quanta_string = "0"
        else:
            quanta_string = " ".join(quanta_list)
        return "|{}>".format(quanta_string)

    def __eq__(self, other):
        if isinstance(other, VibState):
            for key in self.quanta.keys():
                if self.quanta[key] != other.quanta[key]:
                    return False
            return True
        else:
            return False

    def energy(self, MHz=False):
        if MHz:
            result = sum([v0(i)*self.quanta[i] for i in range(starting_mode, n_modes)])/1000000
        else:
            result = sum([w0(i)*self.quanta[i] for i in range(starting_mode, n_modes)])
        return result


def vib_basis(max_energy=max_state_energy):
    max_quanta = [int(max_energy // w0(i)) for i in range(starting_mode, n_modes)]
    ranges = [range(0, i+1) for i in max_quanta]
    basis_quanta_list = product(*ranges)
    vib_state_list = []
    for quanta in basis_quanta_list:
        quanta_dict = {i+6: quanta[i] for i in range(len(quanta))}
        state_energy = sum([quanta_dict[i]*w0(i) for i in range(starting_mode, n_modes)])
        if state_energy <= max_energy:
            vib_state_list.append(VibState(quanta_dict))
    result = sorted(vib_state_list, key=lambda i: i.energy())
    return result


def evaluate_qop(mode: int, state1: VibState, state2: VibState):
    for i in range(starting_mode, n_modes):
        if i != mode:
            if state1.quanta[i] != state2.quanta[i]:
                return 0
    state1_mode_quanta = state1.quanta[mode]
    state2_mode_quanta = state2.quanta[mode]
    if state2_mode_quanta == state1_mode_quanta+1:
        return np.sqrt((state1_mode_quanta+1)/2)
    elif state2_mode_quanta == state1_mode_quanta-1:
        return np.sqrt((state1_mode_quanta/2))
    else:
        return 0


def evaluate_pop(mode: int, state1: VibState, state2: VibState):
    for i in range(starting_mode, n_modes):
        if i != mode:
            if state1.quanta[i] != state2.quanta[i]:
                return 0
    state1_mode_quanta = state1.quanta[mode]
    state2_mode_quanta = state2.quanta[mode]
    if state2_mode_quanta == state1_mode_quanta+1:
        return np.sqrt((state1_mode_quanta+1)/2) / 1j
    elif state2_mode_quanta == state1_mode_quanta-1:
        return np.sqrt((state1_mode_quanta/2)) / (-1j)
    else:
        return 0


def evaluate_multi_vib_op(operator_list: list, state1: VibState, state2: VibState, basis: list):
    n_operators = len(operator_list)
    combinations_list = [[1, [state1]]]
    counter = 1
    for operator_type, mode in operator_list:
        new_combinations_list = []
        for coefficient, states in combinations_list:
            left_state = states[-1]
            if counter == n_operators:
                right_states_list = [state2]
            else:
                right_states_list = [i for i in basis]
            for right_state in right_states_list:
                if operator_type == 'q':
                    evaluation = evaluate_qop(mode, left_state, right_state)
                elif operator_type == 'p':
                    evaluation = evaluate_pop(mode, left_state, right_state)
                else:
                    raise TypeError
                if evaluation != 0:
                    new_coefficient = coefficient*evaluation
                    new_states = [i for i in states] + [right_state]
                    new_combinations_list.append([new_coefficient, new_states])
        combinations_list = [i for i in new_combinations_list]
        counter += 1
    return combinations_list


def extract_jop_coefficients(sympy_expression, n_rot_op):
    possible_rot_operators = [[jop(0), jop(1), jop(2)]]*n_rot_op
    if len(possible_rot_operators) == 0:
        # !!!Only guaranteed to work on expanded expressions.!!!
        return sum([
            i for i in sympy_expression.args
            if not any(j.func == jop for j in preorder_traversal(i)
                       if isinstance(j, Function))])
    elif len(possible_rot_operators) == 1:
        possible_rot_operators_list = [[i] for i in possible_rot_operators[0]]
    elif len(possible_rot_operators) > 1:
        possible_rot_operators_list = list(product(*possible_rot_operators))
    else:
        raise ValueError
    shape_tuple = tuple(3 for i in range(n_rot_op))
    extraction_array = np.zeros(shape=shape_tuple, dtype='complex')
    for operators in possible_rot_operators_list:
        operator_product = Mul(*operators)
        indices = [i.args[0] for i in operators]
        coefficient = sympy_expression.as_coefficient(operator_product)
        if coefficient is not None:
            extraction_array.itemset(tuple(indices), coefficient)
    return extraction_array


def caan(ra, rb, va):
    value = (-1)*baan(ra, rb, va)/v0(va)
    return value


class rot_con(Function):
    n_args = 1

    @classmethod
    def eval(cls, i):
        if i.is_Integer:
            return B0(i)


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
                return rot_con(ra)
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


ground = VibState()
nu22 = VibState({6: 1})
nu33 = VibState({7: 1})

mvb = vib_basis()


h21 = GenericTerm(1, 2, 1, 'H').toEquation(term_definitions)
nu22_h21_nu33 = h21.vibrational_matrix_element(nu22, nu33, mvb)

