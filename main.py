from math import factorial
import numpy as np

class Term:
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
        if isinstance(other, Term):
            if self.combinesWith(other):
                self.coeff = self.coeff + other.coeff
                return self
            else:
                return Expression([self, other])
        elif isinstance(other, Expression):
            items_list = other.items
            items_list.append(self)
            return Expression(items_list)
        elif isinstance(other, Commutator):
            return Expression([other, Commutator])
        else:
            raise NotImplementedError

    def __sub__(self, other):
        return self + (other*(-1))

    def __mul__(self, other):
        if isinstance(other, (Term, Expression, Commutator)):
            raise NotImplementedError
        else:
            try:
                new = self
                new.coeff = new.coeff*other
                return new
            except:
                raise TypeError

    def __rmul__(self, other):
        if isinstance(other, (Term, Expression, Commutator)):
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
        if isinstance(other, Term):
            if self.type == other.type and self.vib_order == other.vib_order and self.rot_order == other.rot_order:
                return True
            else:
                return False
        else:
            return False

    def vibCommutator(self, other):
        if isinstance(other, (Term, Commutator)):
            if self.vib_order == 0 or other.vib_order == 0:
                return 0
            else:
                return Commutator(1, self, other, 'V')
        elif isinstance(other, Expression):
            item_list = [self.vibCommutator(item) for item in other.items]
            return Expression(item_list)
        else:
            raise NotImplementedError

    def rotCommutator(self,other):
        if isinstance(other, (Term, Commutator)):
            if self.rot_order == 0 or other.rot_order == 0:
                return 0
            else:
                return Commutator(1, self, other, 'R')
        elif isinstance(other, Expression):
            item_list = [self.rotCommutator(item) for item in other.items]
            return Expression(item_list)
        else:
            raise NotImplementedError


class Expression:
    def __init__(self, items_list):
        self.items = []
        for item in items_list:
            if isinstance(item, Expression):
                for sub_item in item.items:
                    self.items.append(sub_item)
            elif isinstance(item, (Term, Commutator)):
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
        if isinstance(other, (Term, Commutator)):
            return Expression([*self.items, other])
        elif isinstance(other, Expression):
            return Expression([*self.items, *other.items])
        elif other == 0:
            return self
        else:
            raise TypeError('Cannot add {} to Expression object'.format(other))

    def __sub__(self, other):
        return self + (other*(-1))

    def __mul__(self, other):
        if isinstance(other, (Term, Expression, Commutator)):
            raise NotImplementedError
        else:
            try:
                value = Expression([item*other for item in self.items])
                return value
            except:
                raise TypeError('Cannot multiply Expression object by {}'.format(other))

    def __rmul__(self, other):
        return self*other

    def __getitem__(self, item):
        return self.items[item]

    def __eq__(self, other):
        if isinstance(other, Expression):
            raise NotImplementedError('Need to define expression_sort function first.')
        else:
            return False


class Commutator:
    def __init__(self, coeff, term1, term2, type: str):
        if type == 'V' or type == 'R':
            self.type = type
        else:
            raise TypeError
        if isinstance(term1, (Term, Commutator)) and isinstance(term2, (Term, Commutator)):
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
        elif isinstance(term1, Expression) or isinstance(term2, Expression):
            raise TypeError('Use object commutator function to evaluate the commutator with an expression.')
        else:
            raise NotImplementedError

    def __repr__(self):
        return "{}*[{},{}]{}".format(self.coeff, self.term1, self.term2, self.type)

    def __add__(self, other):
        if isinstance(other, Commutator):
            if self.combinesWith(other):
                new = self
                new.coeff = new.coeff + other.coeff
                return new
            else:
                return Expression([self, other])
        elif isinstance(other, (Term, Expression)):
            return Expression([self, other])
        else:
            raise TypeError

    def __sub__(self, other):
        return self + (other * (-1))

    def __mul__(self, other):
        if isinstance(other, (Term, Expression, Commutator)):
            raise NotImplementedError
        else:
            try:
                new = self
                new.coeff = new.coeff * other
                return new
            except:
                raise TypeError

    def __rmul__(self, other):
        if isinstance(other, (Term, Expression, Commutator)):
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
        if isinstance(other, Commutator):
            if self.term1 == other.term1 and self.term2 == other.term2 and self.type == other.type:
                return True
            else:
                return False
        else:
            return False



M = 2
N = 2

max_order = M + N - 2

transforms = []
for m in range(1,M+1):
    n_start = max([3-m,0])
    for n in range(n_start,N+1):
        transforms.append(Term(1,m,n,'S'))

n_transforms = len(transforms)


def H(i):
    if i == 0:
        value = Term(1,2,0,'H')
    elif i == 1:
        value = Term(1,0,2,'H') + Term(1,3,0,'H') + Term(1,2,1,'H') + Term(1,1,2,'H')
    elif i > 1:
        value = Term(1,i+2,0,'H') + Term(1,i+1,1,'H') + Term(1,i,2,'H')
    else:
        raise AttributeError('Order of {} is not allowed'.format(i))
    return value


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
                commutator = Expression([s_term.vibCommutator(transformedH(k, level - 1)),
                                         s_term.rotCommutator(transformedH(k, level - 1))])
                counter = 1
                while counter < i - k:
                    commutator = Expression([s_term.vibCommutator(commutator), s_term.rotCommutator(commutator)])
                    counter += 1
                value = Expression([commutator * coeff, value])
            return value
        else:
            raise AttributeError
    else:
        raise AttributeError
