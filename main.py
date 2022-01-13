from sympy import I, symbols, Function, preorder_traversal, signsimp, Add, Mul, Pow, sympify, Symbol, Rational
from sympy.core.function import UndefinedFunction
from itertools import permutations
from math import factorial


def printDictionary(dictionary: dict):
    for key, value in dictionary.items():
        print('{}:  {}'.format(key, value))


jop = Function('jop')  # j rotational operator
pop = Function('pop')  # p (momentum) vibrational operator
qop = Function('qop')  # q (position) vibrational operator
ee = Function('ee')  # Levi-Cevita symbol
lop = Function('lop')  # L (ladder) vibrational operator
sigma = Function('sigma')  # sign (+/-)
omega = Function('omega')  # harmonic frequency
VC = Function('VC')  # Vibrational Commutator
RC = Function('RC')  # Rotational Commutator
D = Function('D')  # Denominator
nvm = Symbol('nvm')  # Number of Vibrational Modes

used_symbols = {i: set() for i in ['v', 'r', 'A', 'B', 'H', 'S']}


def v(i: int):
    symbol = symbols('v{}'.format(i))
    used_symbols['v'].add(symbol)
    return symbol


def r(i: int):
    symbol = symbols('r{}'.format(i))
    used_symbols['r'].add(symbol)
    return symbol


def functionGenerator(base_string: str):
    def make_function(m: int, n: int):
        function = Function('{}{}{}'.format(base_string, m, n))
        try:
            used_symbols[base_string].add(function)
        except KeyError:
            used_symbols[base_string] = set(function)
        return function

    return make_function


A = functionGenerator('A')
B = functionGenerator('B')


def symbolGenerator(base_string: str):
    def make_symbol(m: int, n: int):
        symbol = symbols('{}{}{}'.format(base_string, m, n))
        try:
            used_symbols[base_string].add(symbol)
        except KeyError:
            used_symbols[base_string] = set(symbol)
        return symbol

    return make_symbol


H = symbolGenerator('H')
S = symbolGenerator('S')


def gH(m: int, n: int):
    return GenericTerm(1, m, n, 'H')


def gS(m: int, n: int):
    return GenericTerm(1, m, n, 'S')


def printUsedSymbols(used_symbols_dict):
    for key, value in used_symbols_dict.items():
        ordered_list = sorted(value, key=lambda i: str(i))
        print('{}: {}'.format(key, ordered_list))


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


def targetExpression(i: int, j: int):
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
    for m in range(1, i+1):
        n_start = max([3-m, 0])
        for n in range(n_start, j+1):
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

    return select_orders, defining_equations


def printDefiningEquations(defining_equations):
    for item in defining_equations:
        s_term = item[0]
        other_part = item[1]
        vib_order = s_term.vib_order
        rot_order = s_term.rot_order
        print('H_t({},{}) = {} + I*VC({}, H20)'.format(vib_order, rot_order, other_part.sympy(), s_term))


class Term:
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
                    new_smaller_vib_term.coefficient = (1/nvm)**diff_n_vib_indices*new_smaller_vib_term.coefficient
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
                    new_smaller_vib_term.coefficient = (1/nvm)**diff_n_vib_indices*new_smaller_vib_term.coefficient
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


def transform_solution(defining_expression, base_name: Symbol):
    if isinstance(defining_expression, Term):
        ladder_expression = Expression([defining_expression]).toLadder()
    elif isinstance(defining_expression, Expression):
        ladder_expression = defining_expression.toLadder()
    else:
        raise TypeError

    final_ladder_terms = []

    for i in range(0, len(ladder_expression)):
        ladder_term = ladder_expression[i]
        reversed_ladder_term = ladder_term.changeVibIndices(list(reversed(ladder_term.vib_indices)))
        forward_coefficient = ladder_term.coefficient
        reversed_coefficient = reversed_ladder_term.coefficient
        operator_indices = [list(preorder_traversal(x))[1] for x in ladder_term.vib_op]
        denominator = D(*[ii for jj in [[sigma(index), omega(index)] for index in operator_indices] for ii in jj])
        new_coefficient = -I*Rational(1, 2)*(forward_coefficient+reversed_coefficient)/denominator
        new_term = LadderTerm(ladder_term.vib_op,
                              ladder_term.rot_op,
                              new_coefficient,
                              ladder_term.vib_indices,
                              ladder_term.rot_indices)
        final_ladder_terms.append(new_term)
    final_ladder_expression = LadderExpression(final_ladder_terms)
    new_operator_expression = final_ladder_expression.toOperator()
    subbed_terms = []
    substitution_definitions = {}
    for term in new_operator_expression:
        term.coefficient = signsimp(term.coefficient)
    new_operator_expression = Expression(new_operator_expression.items)
    for i in range(0, len(new_operator_expression)):
        term = new_operator_expression[i]
        sub_name = Function(str(base_name)+'_'+str(i))
        full_coefficient = term.coefficient
        vib_indices = []
        rot_indices = []
        for x in preorder_traversal(full_coefficient):
            if isinstance(x, Symbol):
                if str(x)[0] == 'v':
                    if x not in vib_indices:
                        vib_indices.append(x)
                elif str(x)[0] == 'r':
                    if x not in rot_indices:
                        rot_indices.append(x)
        indices = vib_indices_sorter(vib_indices) + rot_indices_sorter(rot_indices)
        new_coefficient = sub_name(*indices)
        substitution_definitions[new_coefficient] = full_coefficient
        new_term = Term(term.vib_op,
                        term.rot_op,
                        new_coefficient,
                        term.vib_indices,
                        term.rot_indices)
        subbed_terms.append(new_term)
    final_expression = Expression(subbed_terms)
    return final_expression, substitution_definitions


def find_transforms(defining_equations_list: list, definitions: dict):
    print('The following transforms will be defined: \n    {}'.format([item[0] for item in defining_equations_list]))
    with_new_definitions = {**definitions}
    sub_definitions = {}
    for transform_term, defining_equation in defining_equations_list:
        if (not isinstance(transform_term, GenericTerm)) or\
                (not isinstance(defining_equation, (GenericTerm, GenericExpression))):
            raise TypeError
        transform_name = transform_term.symbol
        print('{} is being defined ...'.format(transform_name))
        full_defining_equation = defining_equation.toEquation(with_new_definitions)
        transform_definition, transform_sub_definition = transform_solution(
            full_defining_equation, transform_name)
        with_new_definitions[transform_name] = transform_definition
        sub_definitions = {**sub_definitions, **transform_sub_definition}
        print(' ... Done.')

    return with_new_definitions, sub_definitions


def find_target_and_definitions(m: int, n: int, term_definitions: dict, full_simplify=False):
    target, target_defining_equations = targetExpression(m, n)
    with_transform_definitions, transform_sub_definitions = find_transforms(target_defining_equations, term_definitions)
    final_target = target.toEquation(with_transform_definitions)
    if full_simplify:
        subbed_terms = []
        for i in range(len(final_target)):
            term = final_target[i]
            base_name = Function('H{}{}_{}'.format(m, n, i))
            full_coefficient = term.coefficient
            vib_indices = list(set([x for x in preorder_traversal(full_coefficient) if x in used_symbols['v']]))
            rot_indices = list(set([x for x in preorder_traversal(full_coefficient) if x in used_symbols['r']]))
            indices = vib_indices_sorter(vib_indices) + rot_indices_sorter(rot_indices)
            new_coefficient = base_name(*indices)
            transform_sub_definitions[new_coefficient] = full_coefficient
            new_term = Term(term.vib_op,
                            term.rot_op,
                            new_coefficient,
                            term.vib_indices,
                            term.rot_indices)
            subbed_terms.append(new_term)
        final_subbed_target = Expression(subbed_terms)
        return final_subbed_target, with_transform_definitions, transform_sub_definitions
    else:
        return final_target, with_transform_definitions, transform_sub_definitions

