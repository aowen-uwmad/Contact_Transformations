from math import factorial
import numpy as np
from sympy import I, symbols, Function, preorder_traversal, Basic, Expr
from sympy.core.function import UndefinedFunction

vibrational_indices_str = ' '.join(['v{:02d}'.format(i) for i in range(0, 30)])
v = symbols(vibrational_indices_str)

rotational_indices_str = ' '.join(['r{:02d}'.format(i) for i in range(0, 30)])
r = symbols(rotational_indices_str)

A = Function('A')
A20 = Function('A20')
A30 = Function('A30')
A40 = Function('A40')
B40 = Function('B40')
A21 = Function('A21')
A31 = Function('A31')
A02 = Function('A02')
A12 = Function('A12')
A22 = Function('A22')
A32 = Function('A32')
jop = Function('jop')
pop = Function('pop')
qop = Function('qop')


class GenericTerm:
    def __init__(self, coeff, m: int, n: int, type: str):
        if type=='H' or type=='S':
            self.type = type
        else:
            raise AttributeError('Unrecognized type {}.'.format(type))
        self.coeff = coeff
        if m < 0 or n < 0 or (m+n-2) < 0:
            raise AttributeError('Invalid vibration and/or rotation operators.')
        else:
            self.vib_order = m
            self.rot_order = n
            if m == 0 and n == 2:
                self.order = 1
            else:
                self.order = m + n - 2

    def __repr__(self):
        return "{}*{}({},{})".format(self.coeff, self.type, self.vib_order, self.rot_order)

    def __add__(self, other):
        if isinstance(other, GenericTerm):
            if self.combinesWith(other):
                self.coeff = self.coeff + other.coeff
                return self
            else:
                return GenericExpression([self, other])
        elif isinstance(other, GenericExpression):
            items_list = other.items
            items_list.append(self)
            return GenericExpression(items_list)
        elif isinstance(other, GenericCommutator):
            return GenericExpression([other, GenericCommutator])
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
                new.coeff = new.coeff*other
                return new
            except:
                raise TypeError

    def __rmul__(self, other):
        if isinstance(other, (GenericTerm, GenericExpression, GenericCommutator)):
            raise NotImplementedError
        else:
            try:
                new = self
                new.coeff = new.coeff*other
                return new
            except:
                raise TypeError

    def __eq__(self, other):
        if self.combinesWith(other) and self.coeff == other.coeff:
            return True
        else:
            return False

    def combinesWith(self,other):
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

    def rotCommutator(self,other):
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

    def select(self, m, n):
        if self.vib_order == m and self.rot_order == n:
            return self
        else:
            return 0


class GenericExpression:
    def __init__(self, items_list):
        self.items = []
        for item in items_list:
            if isinstance(item, GenericExpression):
                for sub_item in item.items:
                    self.items.append(sub_item)
            elif isinstance(item, (GenericTerm, GenericCommutator)):
                if item.order >= 0 and not ([item.vib_order, item.rot_order] in [[0,0],[1,0],[0,1],[1,1]]):
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
            if item.coeff != 0:
                final_items.append(item)

    def __repr__(self):
        return "{}".format([item for item in self.items])

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
            raise NotImplementedError('Need to define expression_sort function first.')
        else:
            return False

    def select(self, m, n):
        selected = []
        for item in self.items:
            if item.vib_order == m and item.rot_order == n:
                selected.append(item)
        return GenericExpression(selected)


class GenericCommutator:
    def __init__(self, coeff, term1, term2, type: str):
        if type == 'V' or type == 'R':
            self.type = type
        else:
            raise TypeError
        if isinstance(term1, (GenericTerm, GenericCommutator)) and isinstance(term2, (GenericTerm, GenericCommutator)):
            self.coeff = coeff*term1.coeff*term2.coeff
            term1.coeff = 1
            term2.coeff = 1
            self.term1 = term1
            self.term2 = term2
            if self.type == 'V':
                self.vib_order = term1.vib_order + term2.vib_order - 2
                self.rot_order = term1.rot_order + term2.rot_order
            else:
                self.vib_order = term1.vib_order + term2.vib_order
                self.rot_order = term1.rot_order + term2.rot_order - 1
            self.order = self.vib_order + self.rot_order - 2
        elif isinstance(term1, GenericExpression) or isinstance(term2, GenericExpression):
            raise TypeError('Use object commutator function to evaluate the commutator with an expression.')
        else:
            raise NotImplementedError

    def __repr__(self):
        return "{}*[{},{}]{}".format(self.coeff, self.term1, self.term2, self.type)

    def __add__(self, other):
        if isinstance(other, GenericCommutator):
            if self.combinesWith(other):
                new = self
                new.coeff = new.coeff + other.coeff
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
        if isinstance(other, (GenericTerm, GenericExpression, GenericCommutator)):
            raise NotImplementedError
        else:
            try:
                new = self
                new.coeff = new.coeff * other
                return new
            except:
                raise TypeError

    def __rmul__(self, other):
        if isinstance(other, (GenericTerm, GenericExpression, GenericCommutator)):
            raise NotImplementedError
        else:
            try:
                new = self
                new.coeff = new.coeff * other
                return new
            except:
                raise TypeError

    def __eq__(self, other):
        if self.combinesWith(other) and self.coeff == other.coeff:
            return True
        else:
            return False

    def combinesWith(self, other):
        if isinstance(other, GenericCommutator):
            if self.term1 == other.term1 and self.term2 == other.term2 and self.type == other.type:
                return True
            else:
                return False
        else:
            return False

    def select(self, m, n):
        if self.vib_order == m and self.rot_order == n:
            return self
        else:
            return 0


class Term:
    def __init__(self, vib_op_list: list, rot_op_list: list, coefficient, vib_indices: list, rot_indices: list):
        self._n_vib_op = len(vib_op_list)
        self._n_rot_op = len(rot_op_list)
        self._vib_indices = vib_indices
        self._rot_indices = rot_indices
        self._vib_op = vib_op_list
        self._rot_op = rot_op_list
        self._coeff = coefficient

    @property
    def vib_op(self):
        return self._vib_op

    @property
    def rot_op(self):
        return self._rot_op

    @property
    def coeff(self):
        return self._coeff

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
        return "Sum(({})*{})".format(self.coeff, "*".join(operators))

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

            left_vib_indices = v[0:n_left_vib_indices]
            right_vib_indices = v[n_left_vib_indices:(n_left_vib_indices+n_right_vib_indices)]
            left_rot_indices = r[0:n_left_rot_indices]
            right_rot_indices = r[n_left_rot_indices:(n_left_rot_indices+n_right_rot_indices)]

            new_left = self.changeIndices(left_vib_indices, left_rot_indices)
            new_right = other.changeIndices(right_vib_indices, right_rot_indices)

            final_vib_indices = left_vib_indices + right_vib_indices
            final_rot_indices = left_rot_indices + right_rot_indices
            final_coeff = new_left.coeff * new_right.coeff
            final_vib_op = new_left.vib_op + new_right.vib_op
            final_rot_op = new_left.rot_op + new_right.rot_op

            return Term(final_vib_op, final_rot_op, final_coeff, final_vib_indices, final_rot_indices)
        else:
            try:
                new_coefficient = self.coeff*other
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
                    self.coeff == other.coeff]):
                return True
            else:
                return False
        else:
            return False

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
        new_coeff = self.coeff.subs(substitution_rules, simultaneous=True)
        return Term(new_vib_op, self.rot_op, new_coeff, new_vib_indices_list, self.rot_indices)

    def changeRotIndices(self, new_rot_indices_list):
        if len(new_rot_indices_list) != len(self.rot_indices):
            raise ValueError('Unequal length of rot indices lists.')
        substitution_rules = {}
        for i in range(0, len(self.rot_indices)):
            substitution_rules[self.rot_indices[i]] = new_rot_indices_list[i]
        new_rot_op = [x.subs(substitution_rules, simultaneous=True) for x in self.rot_op]
        new_coeff = self.coeff.subs(substitution_rules, simultaneous=True)
        return Term(self.vib_op, new_rot_op, new_coeff, self.vib_indices, new_rot_indices_list)

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
            # Changing smaller term to match indices of larger term.  That way the additional indices of the larger
            # term are unaffected.  If same, then doesn't matter.
            combined_vib_indices = larger_vib_term.vib_indices
            combined_vib_op = larger_vib_term.vib_op
            vib_substitutions = {}
            for i in range(0, smaller_vib_term.n_vib_op):
                smaller_vib_index = [x for x in preorder_traversal(smaller_vib_term.vib_op[i])][1]
                larger_vib_index = [x for x in preorder_traversal(larger_vib_term.vib_op[i])][1]
                vib_substitutions[smaller_vib_index] = larger_vib_index
            new_vib_indices = [x.subs(vib_substitutions, simultaneous=True) for x in smaller_vib_term.vib_indices]
            new_smaller_vib_term = smaller_vib_term.changeVibIndices(new_vib_indices)
            if len(new_smaller_vib_term.rot_indices) <= len(larger_vib_term.rot_indices):
                smaller_rot_term = new_smaller_vib_term
                larger_rot_term = larger_vib_term
            else:
                smaller_rot_term = larger_vib_term
                larger_rot_term = new_smaller_vib_term
            combined_rot_indices = larger_rot_term.rot_indices
            combined_rot_op = larger_rot_term.rot_op
            rot_substitutions = {}
            for i in range(0, smaller_rot_term.n_rot_op):
                smaller_rot_index = [x for x in preorder_traversal(smaller_rot_term.rot_op[i])][1]
                larger_rot_index = [x for x in preorder_traversal(larger_rot_term.rot_op[i])][1]
                rot_substitutions[smaller_rot_index] = larger_rot_index
            new_rot_indices = [x.subs(rot_substitutions, simultaneous=True) for x in smaller_rot_term.rot_indices]
            new_smaller_rot_term = smaller_rot_term.changeRotIndices(new_rot_indices)
            combined_coefficient = larger_rot_term.coeff + new_smaller_rot_term.coeff
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

                left_vib_indices = v[0:n_left_vib_indices]
                right_vib_indices = v[n_left_vib_indices:(n_left_vib_indices + n_right_vib_indices)]
                left_rot_indices = r[0:n_left_rot_indices]
                right_rot_indices = r[n_left_rot_indices:(n_left_rot_indices + n_right_rot_indices)]

                new_left = self.changeIndices(left_vib_indices, left_rot_indices)
                new_right = other.changeIndices(right_vib_indices, right_rot_indices)

                final_vib_indices = left_vib_indices + right_vib_indices
                final_rot_indices = left_rot_indices + right_rot_indices
                final_coeff = (1/2) * new_left.coeff * new_right.coeff

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
                                                multiplier*final_coeff,
                                                final_vib_indices,
                                                final_rot_indices)
                                new_term = new_term.changeVibIndices(kronecker_delta_rules)
                                final_terms.append(new_term)

                return Expression(final_terms)
        else:
            return 0


def pure_vibrational_commutator(left: list, right: list):
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
                return [[], kronecker_delta_rules, multiplier]
            elif a.func == pop and b.func == qop:
                kronecker_delta_rules = {a_index: b_index}
                multiplier = (-1)*I
                return [[], kronecker_delta_rules, multiplier]
            else:
                raise ValueError
        elif len(right) > 1:
            pass
        else:
            raise ValueError
    elif len(left) > 1:
        pass
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
            if term.coeff != 0:
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
        if isinstance(other, (Term, Expression)):
            raise NotImplementedError
        else:
            try:
                value = Expression([item * other for item in self.items])
                return value
            except:
                raise TypeError('Cannot multiply GenericExpression object by {}'.format(other))

    def __rmul__(self, other):
        return self*other

    def __getitem__(self, item):
        return self.items[item]

    def __eq__(self, other):
        if isinstance(other, Expression):
            raise NotImplementedError('Need to define expression_sort function first.')
        else:
            return False


h20 = (Term([qop(v[0]), qop(v[1])], [], A20(v[0], v[1]), [v[0], v[1]], []) +
       Term([pop(v[0]), pop(v[1])], [], A20(v[0], v[1]), [v[0], v[1]], []))
h30 = Term([qop(v[0]), qop(v[1]), qop(v[2])], [], A30(v[0], v[1], v[2]), [v[0], v[1], v[2]], [])
h40 = (Term([qop(v[0]), qop(v[1]), qop(v[2]), qop(v[3])], [], A40(v[0], v[1], v[2], v[3]), [v[0], v[1], v[2], v[3]], []) +
       Term([qop(v[0]), pop(v[1]), qop(v[2]), pop(v[3])], [], B40(v[0], v[1], v[2], v[3]), [v[0], v[1], v[2], v[3]], []))
h21 = Term([qop(v[0]), pop(v[1])], [jop(r[0])], A21(v[0], v[1], r[0]), [v[0], v[1]], [r[0]])
h31 = Term([qop(v[0]), pop(v[1]), qop(v[2])], [jop(r[0])], A31(v[0], v[1], v[2], r[0]), [v[0], v[1], v[2]], [r[0]])
h02 = Term([], [jop(r[0]), jop(r[1])], A02(r[0], r[1]), [], [r[0], r[1]])
h12 = Term([qop(v[0])], [jop(r[0]), jop(r[1])], A12(v[0], r[0], r[1]), [v[0]], [r[0], r[1]])
h22 = Term([qop(v[0]), qop(v[1])], [jop(r[0]), jop(r[1])], A22(v[0], v[1], r[0], r[1]), [v[0], v[1]], [r[0], r[1]])
h32 = Term([qop(v[0]), qop(v[1]), qop(v[2])], [jop(r[0]), jop(r[1])], A32(v[0], v[1], v[2], r[0], r[1]), [v[0], v[1], v[2]], [r[0], r[1]])

evil_term = Term([qop(v[4]),qop(v[3]),qop(v[2])],[],A32(v[0],v[1],v[2],v[3],v[4]),[v[0],v[1],v[2],v[3],v[4]],[])


def H(i):
    if i == 0:
        value = GenericTerm(1, 2, 0, 'H')
    elif i == 1:
        value = GenericTerm(1, 0, 2, 'H') + GenericTerm(1, 3, 0, 'H') + GenericTerm(1, 2, 1, 'H') + GenericTerm(1, 1, 2, 'H')
    elif i > 1:
        value = GenericTerm(1, i + 2, 0, 'H') + GenericTerm(1, i + 1, 1, 'H') + GenericTerm(1, i, 2, 'H')
    else:
        raise AttributeError('Order of {} is not allowed'.format(i))
    return value


def targetExpression(M, N):
    max_order = M + N - 2
    transforms = []
    for m in range(1, M+1):
        n_start = max([3-m, 0])
        for n in range(n_start, N+1):
            transforms.append(GenericTerm(1, m, n, 'S'))
    n_transforms = len(transforms)

    def transformedH(i, level):
        if level == 0:
            return H(i)
        elif level > 0:
            if i == 0:
                return H(0)
            elif i > 0:
                value = transformedH(i, level-1)
                s_term = transforms[level-1]
                for k in range(0, i):
                    coeff = ((1j)**(i - k) / (factorial(i - k)))
                    commutator = GenericExpression([s_term.vibCommutator(transformedH(k, level - 1)),
                                                    s_term.rotCommutator(transformedH(k, level - 1))])
                    counter = 1
                    while counter < i - k:
                        commutator = GenericExpression([s_term.vibCommutator(commutator), s_term.rotCommutator(commutator)])
                        counter += 1
                    value = GenericExpression([commutator * coeff, value])
                return value
            else:
                raise AttributeError
        else:
            raise AttributeError

    up_to_max_order = GenericExpression([transformedH(i, n_transforms) for i in range(0, max_order + 1)])
    select_orders = up_to_max_order.select(M, N)
    defining_equations = []
    for level in range(0, n_transforms):
        transform = transforms[level]
        a = transform.vib_order
        b = transform.rot_order
        term_to_block_diagonalize = GenericExpression([transformedH(i, level) for i in range(0, max_order + 1)]).select(a, b)
        defining_equations.append([transform, term_to_block_diagonalize])

    return select_orders, defining_equations

