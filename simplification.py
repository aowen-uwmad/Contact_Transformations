from main import *

AMN_permutations = {
    A(3, 0): [A(3, 0)(v(0), v(1), v(2)),  # canonical form
              {A(3, 0)(*list(x)): A(3, 0)(v(0), v(1), v(2)) for x in permutations([v(0), v(1), v(2)])}],  # permutations
    A(4, 0): [A(4, 0)(v(0), v(1), v(2), v(3)),
              {A(4, 0)(*list(x)): A(4, 0)(v(0), v(1), v(2), v(3)) for x in permutations([v(0), v(1), v(2), v(3)])}],
    A(3, 1): [A(3, 1)(v(0), v(1), v(2), r(0)),
              {A(3, 1)(v(2), v(1), v(0), r(0)): A(3, 1)(v(0), v(1), v(2), r(0))}],
    A(1, 2): [A(1, 2)(v(0), r(0), r(1)),
              {A(1, 2)(v(0), r(1), r(0)): A(1, 2)(v(0), r(0), r(1))}],
    A(2, 2): [A(2, 2)(v(0), v(1), r(0), r(1)),
              {A(2, 2)(v(1), v(0), r(0), r(1)): A(2, 2)(v(0), v(1), r(0), r(1)),
               A(2, 2)(v(0), v(1), r(1), r(0)): A(2, 2)(v(0), v(1), r(0), r(1)),
               A(2, 2)(v(1), v(0), r(1), r(0)): A(2, 2)(v(0), v(1), r(0), r(1))}],
    A(3, 2): [A(3, 2)(v(0), v(1), v(2), r(0), r(1)),
              {A(3, 2)(v(2), v(1), v(0), r(0), r(1)): A(3, 2)(v(0), v(1), v(2), r(0), r(1)),
               A(3, 2)(v(0), v(1), v(2), r(1), r(0)): A(3, 2)(v(0), v(1), v(2), r(0), r(1)),
               A(3, 2)(v(2), v(1), v(0), r(1), r(0)): A(3, 2)(v(0), v(1), v(2), r(0), r(1))}]
}


def AMN_simplify(some_object, permutations_dictionary: dict):
    def do_simplify(object_coefficient):
        sub_rules = {}
        for key in permutations_dictionary.keys():
            key_terms = []
            for item in preorder_traversal(object_coefficient):
                if item.func == key:
                    key_terms.append(item)
            if len(key_terms) <= 1:
                continue
            else:
                canonical_terms = key_terms[0]
                for term in key_terms[1:]:
                    if term not in canonical_terms:
                        substituted_terms = [term.subs(j, k) for j, k in permutations_dictionary[key][1].items()]
                        new_canonical_form = sorted([term, *substituted_terms], key=lambda i: str(i))[0]
                        canonical_terms.append(new_canonical_form)
                for term in canonical_terms:
                    permutation_args = list(permutations_dictionary[key][0].args)
                    new_args = list(term.args)
                    transform_rules = {permutation_args[i]: new_args[i] for i in range(len(new_args))}
                    for j, k in permutations_dictionary[key][1].items():
                        new_j = j.subs(transform_rules, simultaneous=True)
                        new_k = k.subs(transform_rules, simultaneous=True)
                        sub_rules[new_j] = new_k
        simplified_coefficient = object_coefficient.subs(sub_rules, simultaneous=True)
        return simplified_coefficient

    if isinstance(some_object, Term):
        new_coefficient = do_simplify(some_object.coefficient)
        return Term(some_object.vib_op,
                    some_object.rot_op,
                    new_coefficient,
                    some_object.vib_indices,
                    some_object.rot_indices)
    elif isinstance(some_object, LadderTerm):
        new_coefficient = do_simplify(some_object.coefficient)
        return LadderTerm(some_object.vib_op,
                          some_object.rot_op,
                          new_coefficient,
                          some_object.vib_indices,
                          some_object.rot_indices)
    elif isinstance(some_object, Expression):
        return Expression([AMN_simplify(item, permutations_dictionary) for item in some_object.items])
    elif isinstance(some_object, LadderExpression):
        return LadderExpression([AMN_simplify(item, permutations_dictionary) for item in some_object.items])
    elif isinstance(some_object, (GenericTerm, GenericCommutator, GenericExpression)):
        raise NotImplementedError
    else:
        try:
            return do_simplify(some_object)
        except:
            raise TypeError('Unrecognized object type.')


