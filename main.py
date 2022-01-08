from math import factorial
import numpy as np
from sympy import I, symbols, Function, preorder_traversal, Basic, Expr, signsimp, Add, Mul, Pow
from sympy.core.function import UndefinedFunction
from itertools import permutations

# vibrational_indices_str = ' '.join(['v{:02d}'.format(i) for i in range(0, 30)])
v = [symbols('v{:02d}'.format(i)) for i in range(0, 30)]

# rotational_indices_str = ' '.join(['r{:02d}'.format(i) for i in range(0, 30)])
r = [symbols('r{:02d}'.format(i)) for i in range(0, 30)]

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
ee = Function('ee')
lop = Function('lop')
sigma = Function('sigma')
omega = Function('omega')


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

    @coeff.setter
    def coeff(self, value):
        self._coeff = value

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
        new_coeff = self.coeff.subs(substitution_rules, simultaneous=True)
        new_vib_indices = []
        for x in new_vib_indices_list:
            if x not in new_vib_indices:
                new_vib_indices.append(x)
        return Term(new_vib_op, self.rot_op, new_coeff, new_vib_indices, self.rot_indices)

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
            if len(new_smaller_vib_term.rot_indices) <= len(larger_vib_term.rot_indices):
                smaller_rot_term = new_smaller_vib_term
                larger_rot_term = larger_vib_term
            else:
                smaller_rot_term = larger_vib_term
                larger_rot_term = new_smaller_vib_term
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

                left_vib_indices = v[0:n_left_vib_indices]
                right_vib_indices = v[n_left_vib_indices:(n_left_vib_indices + n_right_vib_indices)]
                left_rot_indices = r[0:n_left_rot_indices]
                right_rot_indices = r[n_left_rot_indices:(n_left_rot_indices + n_right_rot_indices)]

                new_left = self.changeIndices(left_vib_indices, left_rot_indices)
                new_right = other.changeIndices(right_vib_indices, right_rot_indices)

                final_vib_indices = left_vib_indices + right_vib_indices
                final_rot_indices = left_rot_indices + right_rot_indices
                final_coeff = (1/2) * new_left.coeff * new_right.coeff

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
                                                multiplier*final_coeff,
                                                final_vib_indices,
                                                final_rot_indices + list(new_indices))
                                final_terms.append(new_term)
                final_expression = Expression(final_terms)
                for term in final_expression:
                    ee_atoms = []
                    for atom in list(preorder_traversal(term.coeff)):
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
                    new_coefficient = term.coeff.subs(ee_replace_rules, simultaneous=True)
                    term.coeff = new_coefficient
                return final_expression
        else:
            return 0

    def toLadder(self):
        new_coefficient = self.coeff
        ladder_operators = []
        for vib_op in self.vib_op:
            index = list(preorder_traversal(vib_op))[1]
            ladder_operators.append(lop(index, sigma(index)))
            if vib_op.func == qop:
                new_coefficient = new_coefficient*(1/2)
            elif vib_op.func == pop:
                new_coefficient = new_coefficient*(1/2)*I*sigma(index)
        return LadderTerm(ladder_operators, self.rot_op, new_coefficient, self.vib_indices, self.rot_indices)

    def printProperties(self):
        print('vib_op (n={}): {}'.format(self.n_vib_op, self.vib_op))
        print('rot_op (n={}): {}'.format(self.n_rot_op, self.rot_op))
        print('vib_indices: {}'.format(self.vib_indices))
        print('rot_indices: {}'.format(self.rot_indices))
        print('coefficient: {}'.format(self.coeff))


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
                unused_indices = [i for i in r if i not in used_indices]
                new_index = unused_indices[0]
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

    @coeff.setter
    def coeff(self, value):
        self._coeff = value

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

            return LadderTerm(final_vib_op, final_rot_op, final_coeff, final_vib_indices, final_rot_indices)
        else:
            try:
                new_coefficient = self.coeff*other
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
                    self.coeff == other.coeff]):
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
        new_coeff = self.coeff.subs(substitution_rules, simultaneous=True)
        new_vib_indices = []
        for x in new_vib_indices_list:
            if x not in new_vib_indices:
                new_vib_indices.append(x)
        return LadderTerm(new_vib_op, self.rot_op, new_coeff, new_vib_indices, self.rot_indices)

    def changeRotIndices(self, new_rot_indices_list):
        if len(new_rot_indices_list) != len(self.rot_indices):
            raise ValueError('Unequal length of rot indices lists.')
        substitution_rules = {}
        for i in range(0, len(self.rot_indices)):
            substitution_rules[self.rot_indices[i]] = new_rot_indices_list[i]
        new_rot_op = [x.subs(substitution_rules, simultaneous=True) for x in self.rot_op]
        new_coeff = self.coeff.subs(substitution_rules, simultaneous=True)
        return LadderTerm(self.vib_op, new_rot_op, new_coeff, self.vib_indices, new_rot_indices_list)

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
            combined_term = LadderTerm(combined_vib_op,
                                       combined_rot_op,
                                       combined_coefficient,
                                       combined_vib_indices,
                                       combined_rot_indices)
            return combined_term
        else:
            return LadderExpression([self, other])

    def toOperator(self):
        transform_list = [[self.coeff, self.vib_op]]
        for index in self.vib_indices:
            new_transform_list = []
            for pair in transform_list:
                coeff = pair[0]
                ops = pair[1]
                new_transform_list.append([coeff.subs({sigma(index): 1}, simultaneous=True),
                                           [op.subs({sigma(index): 1}, simultaneous=True) for op in ops]])
                new_transform_list.append([coeff.subs({sigma(index): -1}, simultaneous=True),
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
        print('coefficient: {}'.format(self.coeff))


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
            if term.coeff != 0:
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


def transform_solution(operator_expression):
    if isinstance(operator_expression, Term):
        ladder_expression = Expression([operator_expression]).toLadder()
    elif isinstance(operator_expression, Expression):
        ladder_expression = operator_expression.toLadder()
    else:
        raise TypeError
    final_ladder_terms = []
    for ladder_term in ladder_expression:
        reversed_ladder_term = ladder_term.changeVibIndices(list(reversed(ladder_term.vib_indices)))
        forward_coefficient = ladder_term.coeff
        reversed_coefficient = reversed_ladder_term.coeff
        operator_indices = [list(preorder_traversal(x))[1] for x in ladder_term.vib_op]
        denominator = sum([sigma(index)*omega(index) for index in operator_indices])
        new_coefficient = -I*(1/2)*(forward_coefficient+reversed_coefficient)/denominator
        new_term = LadderTerm(ladder_term.vib_op,
                              ladder_term.rot_op,
                              new_coefficient,
                              ladder_term.vib_indices,
                              ladder_term.rot_indices)
        final_ladder_terms.append(new_term)
    final_ladder_expression = LadderExpression(final_ladder_terms)
    new_operator_expression = final_ladder_expression.toOperator()
    for term in new_operator_expression:
        term.coeff = signsimp(term.coeff)
    final_expression = Expression([term for term in new_operator_expression])
    return final_expression


def transform_solution_simplified(operator_expression, base_name: UndefinedFunction):
    if isinstance(operator_expression, Term):
        ladder_expression = Expression([operator_expression]).toLadder()
    elif isinstance(operator_expression, Expression):
        ladder_expression = operator_expression.toLadder()
    else:
        raise TypeError
    final_ladder_terms = []
    for i in range(0, len(ladder_expression)):
        ladder_term = ladder_expression[i]
        reversed_ladder_term = ladder_term.changeVibIndices(list(reversed(ladder_term.vib_indices)))
        forward_coefficient = ladder_term.coeff
        reversed_coefficient = reversed_ladder_term.coeff
        operator_indices = [list(preorder_traversal(x))[1] for x in ladder_term.vib_op]
        denominator = sum([sigma(index)*omega(index) for index in operator_indices])
        new_coefficient = -I*(1/2)*(forward_coefficient+reversed_coefficient)/denominator
        new_term = LadderTerm(ladder_term.vib_op,
                              ladder_term.rot_op,
                              new_coefficient,
                              ladder_term.vib_indices,
                              ladder_term.rot_indices)
        final_ladder_terms.append(new_term)
    final_ladder_expression = LadderExpression(final_ladder_terms)
    new_operator_expression = final_ladder_expression.toOperator()
    subbed_terms = []
    definitions = {}
    for term in new_operator_expression:
        term.coeff = signsimp(term.coeff)
    new_operator_expression = Expression([term for term in new_operator_expression])
    for i in range(0, len(new_operator_expression)):
        term = new_operator_expression[i]
        sub_name = Function(str(base_name)+'_'+str(i))
        full_coefficient = term.coeff
        vib_indices = list(set([x for x in preorder_traversal(full_coefficient) if x in v]))
        rot_indices = list(set([x for x in preorder_traversal(full_coefficient) if x in r]))
        indices = vib_indices_sorter(vib_indices) + rot_indices_sorter(rot_indices)
        new_coefficient = sub_name(*indices)
        definitions[new_coefficient] = full_coefficient
        new_term = Term(term.vib_op,
                        term.rot_op,
                        new_coefficient,
                        term.vib_indices,
                        term.rot_indices)
        subbed_terms.append(new_term)
    final_expression = Expression([term for term in subbed_terms])
    return final_expression, definitions


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


h20 = (Term([qop(v[0]), qop(v[1])], [], A20(v[0], v[1]), [v[0], v[1]], []) +
       Term([pop(v[0]), pop(v[1])], [], A20(v[0], v[1]), [v[0], v[1]], []))


def A20_simplify(term):
    def do_simplify(coefficient):
        A20_terms = []
        for item in preorder_traversal(coefficient):
            if item.func == A20:
                A20_terms.append(item)

        if len(A20_terms) == 0:
            return coefficient

        sub_rules = {}
        for item in A20_terms:
            va, vb = item.args
            str_list = [str(va), str(vb)]
            if str_list != sorted(str_list):
                sub_rules[item] = (omega(va) / omega(vb)) * A20(vb, va)

        new_coefficient = coefficient.subs(sub_rules, simultaneous=True)
        return new_coefficient

    if isinstance(term, Term):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return Term(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    if isinstance(term, LadderTerm):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return LadderTerm(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, Expression):
        return Expression([A20_simplify(item) for item in term])
    elif isinstance(term, LadderExpression):
        return LadderExpression([A20_simplify(item) for item in term])
    else:
        try:
            coeff = term
            new_coeff = do_simplify(coeff)
            return new_coeff
        except:
            raise TypeError('Unrecognized object type.')


h30 = Term([qop(v[0]), qop(v[1]), qop(v[2])], [], A30(v[0], v[1], v[2]), [v[0], v[1], v[2]], [])


def A30_simplify(term):
    def do_simplify(coefficient):
        A30_terms = []
        for item in preorder_traversal(coefficient):
            if item.func == A30:
                A30_terms.append(item)

        if len(A30_terms) == 0:
            return coefficient

        sub_rules = {}
        for item in A30_terms:
            va, vb, vc = item.args
            if [va, vb, vc] != vib_indices_sorter([va, vb, vc]):
                sub_rules[item] = A30(*vib_indices_sorter([va, vb, vc]))

        new_coefficient = coefficient.subs(sub_rules, simultaneous=True)
        return new_coefficient

    if isinstance(term, Term):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return Term(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, LadderTerm):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return LadderTerm(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, Expression):
        return Expression([A30_simplify(item) for item in term])
    elif isinstance(term, LadderExpression):
        return LadderExpression([A30_simplify(item) for item in term])
    else:
        try:
            coeff = term
            new_coeff = do_simplify(coeff)
            return new_coeff
        except:
            raise TypeError('Unrecognized object type.')


h40 = (Term([qop(v[0]), qop(v[1]), qop(v[2]), qop(v[3])], [], A40(v[0], v[1], v[2], v[3]), [v[0], v[1], v[2], v[3]], []) +
       Term([qop(v[0]), pop(v[1]), qop(v[2]), pop(v[3])], [], B40(v[0], v[1], v[2], v[3]), [v[0], v[1], v[2], v[3]], []))


def A40_simplify(term):
    def do_simplify(coefficient):
        A40_terms = []
        for item in preorder_traversal(coefficient):
            if item.func == A40:
                A40_terms.append(item)

        if len(A40_terms) == 0:
            return coefficient

        sub_rules = {}
        for item in A40_terms:
            indices = list(item.args)
            if indices != vib_indices_sorter(indices):
                sub_rules[item] = A40(*vib_indices_sorter(indices))

        new_coefficient = coefficient.subs(sub_rules, simultaneous=True)
        return new_coefficient

    if isinstance(term, Term):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return Term(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, LadderTerm):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return LadderTerm(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, Expression):
        return Expression([A40_simplify(item) for item in term])
    elif isinstance(term, LadderExpression):
        return LadderExpression([A40_simplify(item) for item in term])
    else:
        try:
            coeff = term
            new_coeff = do_simplify(coeff)
            return new_coeff
        except:
            raise TypeError('Unrecognized object type.')


def B40_simplify(term):
    def do_simplify(coefficient):
        B40_terms = []
        for item in preorder_traversal(coefficient):
            if item.func == B40:
                B40_terms.append(item)

        if len(B40_terms) == 0:
            return coefficient

        sub_rules = {}
        for item in B40_terms:
            indices = list(item.args)
            pair1 = indices[:2]
            pair2 = indices[2:]
            sorted_pair1 = vib_indices_sorter(pair1)
            sorted_pair2 = vib_indices_sorter(pair2)
            first_sorted_pair1 = sorted_pair1[0]
            first_sorted_pair2 = sorted_pair2[0]
            if [first_sorted_pair1, first_sorted_pair2] == vib_indices_sorter([first_sorted_pair1, first_sorted_pair2]):
                sorted_indices = sorted_pair1 + sorted_pair2
                needs_pair_swap = False
            else:
                sorted_indices = sorted_pair2 + sorted_pair1
                needs_pair_swap = True
            va, vb, vc, vd = sorted_indices
            if pair1 == sorted_pair1:
                needs_swap1 = False
            else:
                needs_swap1 = True
            if pair2 == sorted_pair2:
                needs_swap2 = False
            else:
                needs_swap2 = True
            if needs_pair_swap == False and needs_swap1 == False and needs_swap2 == False:
                pass
            elif (needs_pair_swap == False and needs_swap1 == False and needs_swap2 == True)\
                    or (needs_pair_swap == True and needs_swap1 == True and needs_swap2 == False):
                sub_rules[item] = -(omega(vc)/omega(vd))*B40(va, vb, vc, vd)
            elif needs_pair_swap == False and needs_swap1 == True and needs_swap2 == False\
                    or (needs_pair_swap == True and needs_swap1 == False and needs_swap2 == True):
                sub_rules[item] = -(omega(va)/omega(vb))*B40(va, vb, vc, vd)
            elif needs_pair_swap == True and needs_swap1 == False and needs_swap2 == False:
                sub_rules[item] = B40(va, vb, vc, vd)
            else:  # (needs_pair_swap == False and needs_swap1 == True and needs_swap2 == True)
                #    or (needs_pair_swap == True and needs_swap1 == True and needs_swap2 == True)
                sub_rules[item] = (omega(va)*omega(vc)/(omega(vb)*omega(vd)))*B40(va, vb, vc, vd)

        new_coefficient = coefficient.subs(sub_rules, simultaneous=True)
        return new_coefficient

    if isinstance(term, Term):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return Term(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, LadderTerm):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return LadderTerm(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, Expression):
        return Expression([B40_simplify(item) for item in term])
    elif isinstance(term, LadderExpression):
        return LadderExpression([B40_simplify(item) for item in term])
    else:
        try:
            coeff = term
            new_coeff = do_simplify(coeff)
            return new_coeff
        except:
            raise TypeError('Unrecognized object type.')


h21 = Term([qop(v[0]), pop(v[1])], [jop(r[0])], A21(v[0], v[1], r[0]), [v[0], v[1]], [r[0]])


def A21_simplify(term):
    def do_simplify(coefficient):
        A21_terms = []
        for item in preorder_traversal(coefficient):
            if item.func == A21:
                A21_terms.append(item)

        if len(A21_terms) == 0:
            return coefficient

        sub_rules = {}
        for item in A21_terms:
            va, vb, ra = item.args
            if [va, vb] != vib_indices_sorter([va, vb]):
                sub_rules[item] = -(omega(va)/omega(vb))*A21(vb, va, ra)

        new_coefficient = coefficient.subs(sub_rules, simultaneous=True)
        return new_coefficient

    if isinstance(term, Term):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return Term(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, LadderTerm):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return LadderTerm(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, Expression):
        return Expression([A21_simplify(item) for item in term])
    elif isinstance(term, LadderExpression):
        return LadderExpression([A21_simplify(item) for item in term])
    else:
        try:
            coeff = term
            new_coeff = do_simplify(coeff)
            return new_coeff
        except:
            raise TypeError('Unrecognized object type.')


h31 = Term([qop(v[0]), pop(v[1]), qop(v[2])], [jop(r[0])], A31(v[0], v[1], v[2], r[0]), [v[0], v[1], v[2]], [r[0]])


def A31_simplify(term):
    def do_simplify(coefficient):
        A31_terms = []
        for item in preorder_traversal(coefficient):
            if item.func == A31:
                A31_terms.append(item)

        if len(A31_terms) == 0:
            return coefficient

        sub_rules = {}
        for item in A31_terms:
            va, vb, vc, ra = item.args
            if [va, vc] != vib_indices_sorter([va, vc]):
                sub_rules[item] = A31(vc, vb, va, ra)

        new_coefficient = coefficient.subs(sub_rules, simultaneous=True)
        return new_coefficient

    if isinstance(term, Term):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return Term(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, LadderTerm):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return LadderTerm(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, Expression):
        return Expression([A31_simplify(item) for item in term])
    elif isinstance(term, LadderExpression):
        return LadderExpression([A31_simplify(item) for item in term])
    else:
        try:
            coeff = term
            new_coeff = do_simplify(coeff)
            return new_coeff
        except:
            raise TypeError('Unrecognized object type.')


h02 = Term([], [jop(r[0]), jop(r[1])], A02(r[0], r[1]), [], [r[0], r[1]])

# simplification requires introduction of the rotational constants directly into the coefficient expression.

h12 = Term([qop(v[0])], [jop(r[0]), jop(r[1])], A12(v[0], r[0], r[1]), [v[0]], [r[0], r[1]])


def A12_simplify(term):
    def do_simplify(coefficient):
        A12_items = []
        for item in preorder_traversal(coefficient):
            if item.func == A12:
                A12_items.append(item)

        if len(A12_items) == 0:
            return coefficient

        sub_rules = {}
        for item in A12_items:
            va, ra, rb = item.args
            if [ra, rb] != rot_indices_sorter([ra, rb]):
                sub_rules[item] = A12(va, rb, ra)

        new_coefficient = coefficient.subs(sub_rules, simultaneous=True)
        return new_coefficient

    if isinstance(term, Term):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return Term(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, LadderTerm):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return LadderTerm(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, Expression):
        return Expression([A12_simplify(item) for item in term])
    elif isinstance(term, LadderExpression):
        return LadderExpression([A12_simplify(item) for item in term])
    else:
        try:
            coeff = term
            new_coeff = do_simplify(coeff)
            return new_coeff
        except:
            raise TypeError('Unrecognized object type.')


h22 = Term([qop(v[0]), qop(v[1])], [jop(r[0]), jop(r[1])], A22(v[0], v[1], r[0], r[1]), [v[0], v[1]], [r[0], r[1]])


def A22_simplify(term):
    def do_simplify(coefficient):
        A22_items = []
        for item in preorder_traversal(coefficient):
            if item.func == A22:
                A22_items.append(item)

        if len(A22_items) == 0:
            return coefficient

        sub_rules = {}
        for item in A22_items:
            va, vb, ra, rb = item.args
            if [va, vb] == vib_indices_sorter([va, vb]):
                swap_vib = False
            else:
                swap_vib = True
            if [ra, rb] == rot_indices_sorter([ra, rb]):
                swap_rot = False
            else:
                swap_rot = True
            sorted_indices = vib_indices_sorter([va, vb]) + rot_indices_sorter([ra, rb])
            if swap_vib == False and swap_rot == False:
                pass
            else:
                sub_rules[item] = A22(*sorted_indices)

        new_coefficient = coefficient.subs(sub_rules, simultaneous=True)
        return new_coefficient

    if isinstance(term, Term):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return Term(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, LadderTerm):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return LadderTerm(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, Expression):
        return Expression([A22_simplify(item) for item in term])
    elif isinstance(term, LadderExpression):
        return LadderExpression([A22_simplify(item) for item in term])
    else:
        try:
            coeff = term
            new_coeff = do_simplify(coeff)
            return new_coeff
        except:
            raise TypeError('Unrecognized object type.')


h32 = Term([qop(v[0]), qop(v[1]), qop(v[2])], [jop(r[0]), jop(r[1])], A32(v[0], v[1], v[2], r[0], r[1]), [v[0], v[1], v[2]], [r[0], r[1]])


def A32_simplify(term):
    def do_simplify(coefficient):
        A32_terms = []
        for item in preorder_traversal(coefficient):
            if item.func == A32:
                A32_terms.append(item)

        if len(A32_terms) == 0:
            return coefficient

        sub_rules = {}
        for item in A32_terms:
            indices = list(item.args)
            va, vb, vc, ra, rb = indices
            if [va, vc] == vib_indices_sorter([va, vc]):
                sorted_indices = [va, vb, vc] + rot_indices_sorter([ra, rb])
            else:
                sorted_indices = [vc, vb, va] + rot_indices_sorter([ra, rb])

            if indices != sorted_indices:
                sub_rules[item] = A32(*sorted_indices)

        new_coefficient = coefficient.subs(sub_rules, simultaneous=True)
        return new_coefficient

    if isinstance(term, Term):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return Term(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, LadderTerm):
        coeff = term.coeff
        new_coeff = do_simplify(coeff)
        return LadderTerm(term.vib_op, term.rot_op, new_coeff, term.vib_indices, term.rot_indices)
    elif isinstance(term, Expression):
        return Expression([A32_simplify(item) for item in term])
    elif isinstance(term, LadderExpression):
        return LadderExpression([A32_simplify(item) for item in term])
    else:
        try:
            coeff = term
            new_coeff = do_simplify(coeff)
            return new_coeff
        except:
            raise TypeError('Unrecognized object type.')


def AMN_simplify(term):
    for simple_function in [A30_simplify, A40_simplify, A31_simplify, A12_simplify, A22_simplify, A32_simplify]:
        term = simple_function(term)

    # for complex_function in [A20_simplify, B40_simplify, A21_simplify]:
    #     term = complex_function(term)

    return term


def custom_expand(coefficient):
    if isinstance(coefficient, (Add, Mul)):
        power_args = [x for x in coefficient.args if isinstance(x, Pow)]
        nonpower_args = [x.expand(deep=False) for x in coefficient.args if x not in power_args]
        power_part = coefficient.func(*power_args)
        nonpower_part = coefficient.func(*nonpower_args)
        new_coefficient = coefficient.func(nonpower_part, power_part)
        return new_coefficient
    else:
        return coefficient


def AMN_collection(coefficient):
    functions = []
    for item in preorder_traversal(coefficient):
        if isinstance(item, Function):
            if item not in functions:
                functions.append(item)

    AMN_functions = []
    omega_functions = []
    for item in functions:
        if item.func not in [omega, ee]:
            AMN_functions.append(item)
        elif item.func == omega:
            omega_functions.append(item)

    if len(AMN_functions) == 0:
        return coefficient.collect(omega_functions)

    collection_dict = coefficient.collect(AMN_functions, evaluate=False)
    summation = 0
    for key, value in collection_dict.items():
        summation += key*AMN_collection(value)

    return summation


def definition_permutation_finder(variable, definition):
    if not isinstance(variable, Function):
        raise TypeError
    indices = variable.args
    vib_indices = [x for x in indices if x in v]
    rot_indices = [x for x in indices if x in r]
    if (len(vib_indices) + len(rot_indices)) != len(indices):
        raise ValueError('Unrecognized indices')
    sorted_vib_indices = vib_indices_sorter(vib_indices)
    sorted_rot_indices = rot_indices_sorter(rot_indices)
    sub_rules = {}
    for i in range(len(vib_indices)):
        sub_rules[vib_indices[i]] = sorted_vib_indices[i]
    for i in range(len(rot_indices)):
        sub_rules[rot_indices[i]] = sorted_rot_indices[i]
    canonical_variable = variable.subs(sub_rules, simultaneous=True)
    canonical_definition = definition.subs(sub_rules, simultaneous=True)
    permutation_rules = {}
    vib_permutations = [list(x) for x in permutations(sorted_vib_indices)]
    rot_permutations = [list(x) for x in permutations(sorted_rot_indices)]
    for vib_perm in vib_permutations:
        for rot_perm in rot_permutations:
            sub_rules = {}
            if len(vib_perm) > 0:
                for i in range(len(sorted_vib_indices)):
                    sub_rules[sorted_vib_indices[i]] = vib_perm[i]
            if len(rot_perm) > 0:
                for i in range(len(sorted_rot_indices)):
                    sub_rules[sorted_rot_indices[i]] = rot_perm[i]
            perm_variable = canonical_variable.subs(sub_rules, simultaneous=True)
            if perm_variable == canonical_variable:
                continue
            else:
                perm_definition = canonical_definition.subs(sub_rules, simultaneous=True)
                difference = signsimp(canonical_definition) - signsimp(perm_definition)
                combination = signsimp(canonical_definition) + signsimp(perm_definition)
                if difference == 0:
                    permutation_rules[perm_variable] = canonical_variable
                elif combination == 0:
                    permutation_rules[perm_variable] = -1*canonical_variable
    return permutation_rules


def definition_substitution(coefficient, definitions):
    definition_functions = [key.func for key, value in definitions.items()]
    coefficient_functions = []
    for item in preorder_traversal(coefficient):
        try:
            if item.func in definition_functions:
                if item not in coefficient_functions:
                    coefficient_functions.append(item)
        except AttributeError:
            continue
    all_sub_rules = {}
    for key, value in definitions.items():
        for item in coefficient_functions:
            if key.func == item.func:
                vib_indices = [x for x in item.args if x in v]
                rot_indices = [x for x in item.args if x in r]
                sorted_indices = vib_indices_sorter(vib_indices) + rot_indices_sorter(rot_indices)
                if list(item.args) != sorted_indices:
                    sub_rules = {}
                    for i in range(len(key.args)):
                        sub_rules[key.args[i]] = item.args[i]
                    new_value = signsimp(value.subs(sub_rules, simultaneous=True))
                    all_sub_rules[item] = new_value
    new_coefficient = coefficient.subs(all_sub_rules, simultaneous=True)
    return new_coefficient


definitions = {}

s12, s12_definitions = transform_solution_simplified(h12, Function('s12'))
definitions = {**definitions, **s12_definitions}

s21, s21_definitions = transform_solution_simplified(h21, Function('s21'))
definitions = {**definitions, **s21_definitions}


s22_def_a = s21.vibCommutator(h21)*I
s22_def_b = s21.rotCommutator(h02)*I
s22_def_c = s21.vibCommutator(s21.vibCommutator(h20))*(-1/2)
s22_def_d = s12.vibCommutator(h30)*I
s22_def_e = h22

s22_def = (
    s21.vibCommutator(h21)*I +
    s21.rotCommutator(h02)*I +
    s21.vibCommutator(s21.vibCommutator(h20))*(-1/2) +
    s12.vibCommutator(h30)*I +
    h22
)

permutation_rules = {}
for key, value in definitions.items():
    permutation_rules = {**permutation_rules, **definition_permutation_finder(key, value)}

simple_s22_terms = []
for term in s22_def:
    new_term = term
    new_term.coeff = signsimp(AMN_simplify(new_term.coeff).subs(permutation_rules, simultaneous=True)).expand().simplify()
    simple_s22_terms.append(new_term)

simplified_s22_def = Expression(simple_s22_terms)

s22, s22_definitions = transform_solution_simplified(simplified_s22_def, Function('s22'))
definitions = {**definitions, **s22_definitions}

# s22 = transform_solution(s22_def)

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

