from enum import Enum
from .elements.basic import (
    OnePortElement,
    TwoPortElement,
    Transformer,
    Element_I,
    Element_C,
)
from .junctions import Junction, JunctionEqualEffort, JunctionEqualFlow
from .common import Causality
from typing import List, Set
from sympy import Symbol, Function, Equality
import logging


class Bond:
    def __init__(
        self, node_from: OnePortElement | Junction, node_to: OnePortElement | Junction
    ):
        self.node_from = node_from
        self.node_to = node_to
        self.num = None
        self.effort_in_at_to = None
        self.flow_symbol = None
        self.effort_symbol = None


class Graph:
    def __init__(self):
        self.bonds: List[Bond] = []
        self.elements: List[OnePortElement] = []
        self.junctions: List[Junction] = []
        self.two_port_elements: List[TwoPortElement] = []
        self.parameters: Set[Symbol] = set()

    def all_causalities_set(self):
        for bond in self.bonds:
            if bond.effort_in_at_to == None:
                return False
        return True

    def add(self, bond: Bond):
        if isinstance(bond.node_from, OnePortElement):
            if bond.node_from in self.elements:
                raise Exception(
                    f"OnePortElement {bond.node_from} can only be bonded once!"
                )
            if not bond.node_from in self.elements:
                self.elements.append(bond.node_from)
            bond.node_from.bond = bond
            self.parameters.update(bond.node_from.parameter_symbols())
        elif isinstance(bond.node_from, Junction):
            if not bond.node_from in self.junctions:
                self.junctions.append(bond.node_from)
            if not bond in bond.node_from.bonds:
                bond.node_from.bonds.append(bond)
        elif isinstance(bond.node_from, TwoPortElement):
            if not bond.node_from in self.two_port_elements:
                self.two_port_elements.append(bond.node_from)
            bond.node_from.bond_2 = bond
            self.parameters.update(bond.node_from.parameter_symbols())

        if isinstance(bond.node_to, OnePortElement):
            if bond.node_to in self.elements:
                raise Exception(
                    f"OnePortElement {bond.node_to} can only be bonded once!"
                )
            if not bond.node_to in self.elements:
                self.elements.append(bond.node_to)
            bond.node_to.bond = bond
            self.parameters.update(bond.node_to.parameter_symbols())
        elif isinstance(bond.node_to, Junction):
            if not bond.node_to in self.junctions:
                self.junctions.append(bond.node_to)
            if not bond in bond.node_to.bonds:
                bond.node_to.bonds.append(bond)
        elif isinstance(bond.node_to, TwoPortElement):
            if not bond.node_to in self.two_port_elements:
                self.two_port_elements.append(bond.node_to)
            bond.node_to.bond_1 = bond
            self.parameters.update(bond.node_to.parameter_symbols())

        bond.num = len(self.bonds) + 1
        bond.flow_symbol = Function(f"f_{bond.num}")
        bond.effort_symbol = Function(f"e_{bond.num}")

        self.bonds.append(bond)

    def assign_fixed_causalities(self):
        for bond in self.bonds:
            if (
                bond.node_to.causality_policy() == Causality.FixedEffortIn
                or bond.node_from.causality_policy() == Causality.FixedEffortOut
            ):
                bond.effort_in_at_to = True
                logging.debug(
                    f"Set fixed effort-in causality at {bond.node_to} (vs {bond.node_from})"
                )
            elif (
                bond.node_to.causality_policy() == Causality.FixedEffortOut
                or bond.node_from.causality_policy() == Causality.FixedEffortIn
            ):
                bond.effort_in_at_to = False
                logging.debug(
                    f"Set fixed effort-in causality at {bond.node_from} (vs {bond.node_to})"
                )

    def try_assign_constraint_causalities(self):
        something_happened = False
        for junction in self.junctions:
            if isinstance(junction, JunctionEqualEffort):
                for bond in junction.bonds:
                    if (bond.node_to == junction and bond.effort_in_at_to == True) or (
                        bond.node_from == junction and bond.effort_in_at_to == False
                    ):
                        if (
                            junction.effort_in_bond is not None
                            and junction.effort_in_bond is not bond
                        ):
                            raise Exception(
                                f"Causality conflict! Multiple bonds have declared effort-in causality at 0-junction {junction}"
                            )
                        junction.effort_in_bond = bond
                if junction.effort_in_bond is not None:
                    for bond in junction.bonds:
                        if bond is junction.effort_in_bond:
                            continue
                        if bond.node_to == junction and bond.effort_in_at_to == None:
                            bond.effort_in_at_to = False
                            logging.debug(
                                f"Set constraint effort-in causality at {bond.node_from} (vs {bond.node_to}) due to equal-effort junction {junction} "
                            )
                            something_happened = True
                        elif (
                            bond.node_from == junction and bond.effort_in_at_to == None
                        ):
                            bond.effort_in_at_to = True
                            logging.debug(
                                f"Set constraint effort-in causality at {bond.node_to} (vs {bond.node_from}) due to equal-effort junction {junction}"
                            )
                            something_happened = True
            elif isinstance(junction, JunctionEqualFlow):
                for bond in junction.bonds:
                    if (bond.node_to == junction and bond.effort_in_at_to == False) or (
                        bond.node_from == junction and bond.effort_in_at_to == True
                    ):
                        if (
                            junction.effort_out_bond is not None
                            and junction.effort_out_bond is not bond
                        ):
                            raise Exception(
                                f"Causality conflict! Multiple bonds have declared effort-out causality at 1-junction {junction}"
                            )
                        junction.effort_out_bond = bond
                if junction.effort_out_bond is not None:
                    for bond in junction.bonds:
                        if bond is junction.effort_out_bond:
                            continue
                        if bond.node_to == junction and bond.effort_in_at_to == None:
                            bond.effort_in_at_to = True
                            logging.debug(
                                f"Set constraint effort-in causality at {bond.node_to} (vs {bond.node_from}) due to equal-flow junction {junction}"
                            )
                            something_happened = True
                        elif (
                            bond.node_from == junction and bond.effort_in_at_to == None
                        ):
                            bond.effort_in_at_to = False
                            logging.debug(
                                f"Set constraint effort-in causality at {bond.node_from} (vs {bond.node_to}) due to equal-flow junction {junction}"
                            )
                            something_happened = True
        for element in self.two_port_elements:
            if isinstance(element, Transformer):
                trans: Transformer = element
                if (
                    trans.bond_1.effort_in_at_to == True
                    or trans.bond_2.effort_in_at_to == True
                ):
                    if (
                        trans.effort_in_bond is not None
                        and trans.effort_in_bond is not trans.bond_1
                    ):
                        raise Exception(
                            f"Causality conflict! Multiple bonds have declared effort-in causality at transformer {trans}"
                        )
                    trans.effort_in_bond = trans.bond_1
                    if trans.bond_1.effort_in_at_to == None:
                        trans.bond_1.effort_in_at_to = True
                        something_happened = True
                        logging.debug(
                            f"Set constraint effort-in causality at {trans.bond_1.node_to} (vs {trans.bond_1.node_from}) due to transformer {trans}"
                        )
                    if trans.bond_2.effort_in_at_to == None:
                        trans.bond_2.effort_in_at_to = True
                        something_happened = True
                        logging.debug(
                            f"Set constraint effort-in causality at {trans.bond_2.node_to} (vs {trans.bond_2.node_from}) due to transformer {trans}"
                        )
                elif (
                    trans.bond_1.effort_in_at_to == False
                    or trans.bond_2.effort_in_at_to == False
                ):
                    if (
                        trans.effort_in_bond is not None
                        and trans.effort_in_bond is not trans.bond_2
                    ):
                        raise Exception(
                            f"Causality conflict! Multiple bonds have declared effort-in causality at transformer {trans}"
                        )
                    trans.effort_in_bond = trans.bond_2
                    if trans.bond_1.effort_in_at_to == None:
                        trans.bond_1.effort_in_at_to = False
                        something_happened = True
                        logging.debug(
                            f"Set constraint effort-in causality at {trans.bond_1.node_from} (vs {trans.bond_1.node_to}) due to transformer {trans}"
                        )
                    if trans.bond_2.effort_in_at_to == None:
                        trans.bond_2.effort_in_at_to = False
                        something_happened = True
                        logging.debug(
                            f"Set constraint effort-in causality at {trans.bond_2.node_from} (vs {trans.bond_2.node_to}) due to transformer {trans}"
                        )
        return something_happened

    def try_assign_preferred_causality(self):
        something_happened = False
        for element in self.elements:
            if element.bond.effort_in_at_to != None:
                continue
            if element.causality_policy() == Causality.PreferEffortIn:
                if element.bond.node_to == element:
                    element.bond.effort_in_at_to = True
                    logging.debug(
                        f"Set preferred effort-in causality at {element.bond.node_to} (vs {element.bond.node_from})"
                    )
                    something_happened = True
                else:
                    element.bond.effort_in_at_to = False
                    logging.debug(
                        f"Set preferred effort-in causality at {element.bond.node_from} (vs {element.bond.node_to})"
                    )
                    something_happened = True
            elif element.causality_policy() == Causality.PreferEffortOut:
                if element.bond.node_to == element:
                    element.bond.effort_in_at_to = False
                    logging.debug(
                        f"Set preferred effort-out causality at {element.bond.node_to} (vs {element.bond.node_from})"
                    )
                    something_happened = True
                else:
                    element.bond.effort_in_at_to = True
                    logging.debug(
                        f"Set preferred effort-out causality at {element.bond.node_from} (vs {element.bond.node_to})"
                    )
                    something_happened = True
        return something_happened

    def assign_causalities(self):
        self.assign_fixed_causalities()

        while True:
            if self.try_assign_constraint_causalities():
                continue
            elif self.try_assign_preferred_causality():
                continue
            # elif self.try_assign_arbitrary_causality():
            #     continue
            else:
                break

        if self.all_causalities_set():
            logging.debug("Graph is causal!")
        else:
            logging.warning("Graph is not causal...")

    def get_state_equations(self):
        self.assign_causalities()

        time_symbol = Symbol("t")

        state_equations = dict()
        state_variables = dict()
        state_counter = 1
        other_equations = []
        for element in self.elements:
            eq = element.equations(
                element.bond.effort_symbol, element.bond.flow_symbol, time_symbol
            )
            if isinstance(element, Element_I):
                state_equations[element.bond.flow_symbol(time_symbol)] = eq
                state_variables[state_counter] = (
                    element.name,
                    element.bond.flow_symbol(time_symbol),
                )
                state_counter += 1
            elif isinstance(element, Element_C):
                state_equations[element.bond.effort_symbol(time_symbol)] = eq
                state_variables[state_counter] = (
                    element.name,
                    element.bond.effort_symbol(time_symbol),
                )
                state_counter += 1
            else:
                other_equations.append(eq)

        for element in self.two_port_elements:
            other_equations += element.equations(
                element.bond_1.effort_symbol,
                element.bond_2.effort_symbol,
                element.bond_1.flow_symbol,
                element.bond_2.flow_symbol,
                time_symbol,
            )

        for junction in self.junctions:
            if isinstance(junction, JunctionEqualEffort):
                # Add new equation for setting effort-in bond's flow symbol equal to the rest of the flows.
                new_eq = Equality(junction.effort_in_bond.flow_symbol(time_symbol), 0)
                for bond in junction.bonds:
                    if bond is junction.effort_in_bond:
                        continue
                    if junction == bond.node_to:
                        new_eq = Equality(
                            new_eq.lhs, new_eq.rhs + bond.flow_symbol(time_symbol)
                        )
                    else:
                        new_eq = Equality(
                            new_eq.lhs, new_eq.rhs - bond.flow_symbol(time_symbol)
                        )
                other_equations.append(new_eq)
            elif isinstance(junction, JunctionEqualFlow):
                # Add new equation for setting effort-out bond's effort symbol equal to the rest of the efforts
                new_eq = Equality(
                    junction.effort_out_bond.effort_symbol(time_symbol), 0
                )
                for bond in junction.bonds:
                    if bond is junction.effort_out_bond:
                        continue
                    if junction == bond.node_to:
                        new_eq = Equality(
                            new_eq.lhs, new_eq.rhs + bond.effort_symbol(time_symbol)
                        )
                    else:
                        new_eq = Equality(
                            new_eq.lhs, new_eq.rhs - bond.effort_symbol(time_symbol)
                        )
                other_equations.append(new_eq)

        substitutions_made = True
        while substitutions_made:
            substitutions_made = False
            for junction in self.junctions:
                substitutions = []
                if isinstance(junction, JunctionEqualEffort):
                    # Substitute all effort symbols with the effort-in bond's effort symbol.
                    for bond in junction.bonds:
                        if bond is junction.effort_in_bond:
                            continue
                        substitutions.append(
                            (
                                bond.effort_symbol(time_symbol),
                                junction.effort_in_bond.effort_symbol(time_symbol),
                            )
                        )
                elif isinstance(junction, JunctionEqualFlow):
                    # Substitute all flow symbols with the effort-out bond's flow symbol
                    for bond in junction.bonds:
                        if bond is junction.effort_out_bond:
                            continue
                        substitutions.append(
                            (
                                bond.flow_symbol(time_symbol),
                                junction.effort_out_bond.flow_symbol(time_symbol),
                            )
                        )
                for index, eq in enumerate(other_equations):
                    before = other_equations[index]
                    other_equations[index] = Equality(
                        eq.lhs, eq.rhs.subs(substitutions)
                    )
                    if other_equations[index] != before:
                        substitutions_made = True
                for key, val in state_equations.items():
                    before = state_equations[key]
                    state_equations[key] = val.subs(substitutions)
                    if state_equations[key] != before:
                        substitutions_made = True

        ordered_equations = dict()
        for eq in other_equations:
            ordered_equations[eq.lhs] = eq.rhs

        state_symbols = dict()
        state_num = dict()
        state_names = dict()

        for num, (name, var) in state_variables.items():
            state_symbols[var] = Symbol(f"x_{num}")
            state_num[var] = num
            state_names[num] = name

        diff_eq_sys = dict()
        for var, eq in state_equations.items():
            diff_eq = eq.rhs.diff(time_symbol)

            before = None
            while before != diff_eq:
                before = diff_eq
                diff_eq = diff_eq.subs(ordered_equations)

            diff_eq_sys[state_num[var]] = diff_eq.subs(state_symbols)

        return (
            diff_eq_sys,
            [var.subs(state_symbols) for (name, var) in state_variables.values()],
            state_names,
        )
