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


