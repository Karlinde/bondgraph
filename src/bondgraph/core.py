from bondgraph.elements.basic import (
    OnePortElement,
    TwoPortElement,
)
from bondgraph.junctions import Junction, JunctionEqualEffort, JunctionEqualFlow
from bondgraph.common import Causality, Bond
from typing import List, Set
from sympy import Symbol, Function, Equality
import logging

_BG_STATE_INIT = 0
_BG_STATE_CAUSALITIES_DONE = 1


class BondGraph:
    def __init__(self, time_symbol: Symbol):
        self._bonds: List[Bond] = []
        self._elements: List[OnePortElement] = []
        self._junctions: List[Junction] = []
        self._two_port_elements: List[TwoPortElement] = []
        self._parameters: Set[Symbol] = set()
        self._time_symbol = time_symbol
        self._state = _BG_STATE_INIT

    def all_causalities_set(self):
        for bond in self._bonds:
            if not bond.has_causality_set():
                return False
        return True

    def preferred_causalities_valid(self):
        success = True
        for bond in self._bonds:
            if (
                bond.effort_in_at_to == True
                and (
                    bond.node_to.causality_policy() == Causality.PreferEffortOut
                    or bond.node_from.causality_policy() == Causality.PreferEffortIn
                )
            ) or (
                bond.effort_in_at_to == False
                and (
                    bond.node_to.causality_policy() == Causality.PreferEffortIn
                    or bond.node_from.causality_policy() == Causality.PreferEffortOut
                )
            ):
                logging.warning(
                    f"Bond from {bond.node_from.name} to {bond.node_to.name} has non-preferred causality"
                )
                success = False
            if (
                bond.effort_in_at_to == True
                and (
                    bond.node_to.causality_policy() == Causality.FixedEffortOut
                    or bond.node_from.causality_policy() == Causality.FixedEffortIn
                )
            ) or (
                bond.effort_in_at_to == False
                and (
                    bond.node_to.causality_policy() == Causality.FixedEffortIn
                    or bond.node_from.causality_policy() == Causality.FixedEffortOut
                )
            ):
                logging.error(
                    f"Bond from {bond.node_from.name} to {bond.node_to.name} has invalid causality due to fixed causality requirements"
                )
                success = False
        return success

    def add(self, bond: Bond):
        if isinstance(bond.node_from, OnePortElement):
            if bond.node_from in self._elements:
                raise Exception(
                    f"OnePortElement {bond.node_from} can only be bonded once!"
                )
            if not bond.node_from in self._elements:
                self._elements.append(bond.node_from)
            bond.node_from.bond = bond
            self._parameters.update(bond.node_from.parameter_symbols())
        elif isinstance(bond.node_from, Junction):
            if not bond.node_from in self._junctions:
                self._junctions.append(bond.node_from)
            if not bond in bond.node_from.bonds:
                bond.node_from.bonds.append(bond)
        elif isinstance(bond.node_from, TwoPortElement):
            if not bond.node_from in self._two_port_elements:
                self._two_port_elements.append(bond.node_from)
            bond.node_from.bond_2 = bond
            self._parameters.update(bond.node_from.parameter_symbols())

        if isinstance(bond.node_to, OnePortElement):
            if bond.node_to in self._elements:
                raise Exception(
                    f"OnePortElement {bond.node_to} can only be bonded once!"
                )
            if not bond.node_to in self._elements:
                self._elements.append(bond.node_to)
            bond.node_to.bond = bond
            self._parameters.update(bond.node_to.parameter_symbols())
        elif isinstance(bond.node_to, Junction):
            if not bond.node_to in self._junctions:
                self._junctions.append(bond.node_to)
            if not bond in bond.node_to.bonds:
                bond.node_to.bonds.append(bond)
        elif isinstance(bond.node_to, TwoPortElement):
            if not bond.node_to in self._two_port_elements:
                self._two_port_elements.append(bond.node_to)
            bond.node_to.bond_1 = bond
            self._parameters.update(bond.node_to.parameter_symbols())

        bond.num = len(self._bonds) + 1
        bond.flow_symbol = Function(f"f_{bond.num}")
        bond.effort_symbol = Function(f"e_{bond.num}")

        self._bonds.append(bond)

    def assign_fixed_causalities(self):
        for bond in self._bonds:
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
        for junction in self._junctions:
            if junction.assign_constraint_causality():
                something_happened = True
        for element in self._two_port_elements:
            if hasattr(element, "assign_constraint_causality"):
                if element.assign_constraint_causality():
                    something_happened = True
        return something_happened

    def try_assign_preferred_causality(self):
        something_happened = False
        for element in self._elements:
            if hasattr(element, "assign_preferred_causality"):
                if element.assign_preferred_causality():
                    something_happened = True
                    # Only assign causality for one element at a time
                    break
        return something_happened

    def try_assign_arbitrary_causality(self):
        something_happened = False
        for element in self._elements:
            if hasattr(element, "assign_arbitrary_causality"):
                if element.assign_arbitrary_causality():
                    something_happened = True
                    # Only assign causality for one element at a time
                    break
        return something_happened

    def assign_causalities(self):
        self.assign_fixed_causalities()

        while True:
            if self.try_assign_constraint_causalities():
                continue
            elif self.try_assign_preferred_causality():
                continue
            elif self.try_assign_arbitrary_causality():
                continue
            else:
                break

        if self.all_causalities_set():
            logging.debug("All causalities set, graph is causal")
        else:
            logging.error("Graph is not causal")
            raise Exception("Non-causal graph detected")

        if not self.preferred_causalities_valid():
            raise Exception("Unsupported causalities detected")
        self._state = _BG_STATE_CAUSALITIES_DONE

    def get_state_equations(self):
        if self._state < _BG_STATE_CAUSALITIES_DONE:
            self.assign_causalities()

        state_equations = dict()
        state_variables = dict()
        state_counter = 1
        other_equations = []
        logging.debug("Formulating equations for one-port elements...")
        for element in self._elements:
            if hasattr(element, "equations"):
                other_equations += element.equations(
                    element.bond.effort_symbol,
                    element.bond.flow_symbol,
                    self._time_symbol,
                )
            if hasattr(element, "state_equations"):
                state_eqs = element.state_equations(
                    element.bond.effort_symbol,
                    element.bond.flow_symbol,
                    self._time_symbol,
                )
                for eq in state_eqs:
                    if eq[0] in state_equations:
                        raise Exception(f"Duplicate state symbol encountered: {eq[0]}")
                    state_variables[state_counter] = eq[0]
                    state_equations[eq[0]] = Equality(eq[0], eq[1])
                    state_counter += 1

        logging.debug("Formulating equations for two-port elements...")
        for element in self._two_port_elements:
            other_equations += element.equations(
                element.bond_1.effort_symbol,
                element.bond_2.effort_symbol,
                element.bond_1.flow_symbol,
                element.bond_2.flow_symbol,
                self._time_symbol,
            )

        logging.debug("Formulating equations for junctions...")
        for junction in self._junctions:
            if isinstance(junction, JunctionEqualEffort):
                # Add new equation for setting effort-in bond's flow symbol equal to the rest of the flows.
                new_eq = Equality(
                    junction.effort_in_bond.flow_symbol(self._time_symbol), 0
                )
                for bond in junction.bonds:
                    if bond is junction.effort_in_bond:
                        continue
                    if junction == bond.node_to:
                        new_eq = Equality(
                            new_eq.lhs, new_eq.rhs + bond.flow_symbol(self._time_symbol)
                        )
                    else:
                        new_eq = Equality(
                            new_eq.lhs, new_eq.rhs - bond.flow_symbol(self._time_symbol)
                        )
                if junction.effort_in_bond.node_to == junction:
                    new_eq = Equality(new_eq.lhs, -new_eq.rhs)
                other_equations.append(new_eq)
            elif isinstance(junction, JunctionEqualFlow):
                # Add new equation for setting effort-out bond's effort symbol equal to the rest of the efforts
                new_eq = Equality(
                    junction.effort_out_bond.effort_symbol(self._time_symbol), 0
                )
                for bond in junction.bonds:
                    if bond is junction.effort_out_bond:
                        continue
                    if junction == bond.node_to:
                        new_eq = Equality(
                            new_eq.lhs,
                            new_eq.rhs + bond.effort_symbol(self._time_symbol),
                        )
                    else:
                        new_eq = Equality(
                            new_eq.lhs,
                            new_eq.rhs - bond.effort_symbol(self._time_symbol),
                        )
                if junction.effort_out_bond.node_to == junction:
                    new_eq = Equality(new_eq.lhs, -new_eq.rhs)
                other_equations.append(new_eq)

        logging.debug("Substituting in junction equations...")
        substitutions_made = True
        while substitutions_made:
            substitutions_made = False
            for junction in self._junctions:
                substitutions = dict()
                if isinstance(junction, JunctionEqualEffort):
                    # Substitute all effort symbols with the effort-in bond's effort symbol.
                    for bond in junction.bonds:
                        if bond is junction.effort_in_bond:
                            continue
                        substitutions[
                            bond.effort_symbol(self._time_symbol)
                        ] = junction.effort_in_bond.effort_symbol(self._time_symbol)
                elif isinstance(junction, JunctionEqualFlow):
                    # Substitute all flow symbols with the effort-out bond's flow symbol
                    for bond in junction.bonds:
                        if bond is junction.effort_out_bond:
                            continue
                        substitutions[
                            bond.flow_symbol(self._time_symbol)
                        ] = junction.effort_out_bond.flow_symbol(self._time_symbol)
                for index, eq in enumerate(other_equations):
                    before = other_equations[index]
                    other_equations[index] = Equality(
                        eq.lhs, eq.rhs.xreplace(substitutions)
                    )
                    if other_equations[index] != before:
                        substitutions_made = True
                for key, val in state_equations.items():
                    before = state_equations[key]
                    state_equations[key] = val.xreplace(substitutions)
                    if state_equations[key] != before:
                        substitutions_made = True

        ordered_equations = dict()
        for eq in other_equations:
            ordered_equations[eq.lhs] = eq.rhs

        logging.debug("Substituting in other equations...")
        substitutions_made = True
        while substitutions_made:
            substitutions_made = False
            substituted_equations = dict()
            for lhs, rhs in ordered_equations.items():
                before = rhs
                substituted_equations[lhs] = rhs.xreplace(ordered_equations)
                if substituted_equations[lhs] != before:
                    substitutions_made = True
            if substitutions_made:
                ordered_equations = substituted_equations

        state_num = dict()
        for num, var in state_variables.items():
            state_num[var] = num

        logging.debug("Generating differential equations...")
        diff_eq_sys = dict()
        for var, eq in state_equations.items():
            diff_eq = eq.rhs.diff(self._time_symbol)

            before = None
            while before != diff_eq:
                before = diff_eq
                diff_eq = diff_eq.xreplace(ordered_equations)

            diff_eq_sys[state_num[var]] = diff_eq

        return (
            diff_eq_sys,
            list(state_variables.values()),
        )

    def get_nodes(self):
        import itertools

        node_list = []
        for element in itertools.chain(
            self._elements, self._junctions, self._two_port_elements
        ):
            node_list.append(element)
        return node_list
