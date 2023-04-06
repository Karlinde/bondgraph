from bondgraph.elements.basic import (
    OnePortElement,
    TwoPortElement,
)
from bondgraph.junctions import Junction, JunctionEqualEffort, JunctionEqualFlow
from bondgraph.common import Causality, Bond, HasStateEquations, Node
from typing import Dict, List, Set
from sympy import Expr, Symbol, Function, Equality
import logging

_BG_STATE_INIT = 0
_BG_STATE_CAUSALITIES_DONE = 1


def _populate_one_port_equations(
    state_equations: Dict[Symbol, Expr],
    state_variables: Dict[int, Symbol],
    state_counter: int,
    other_equations: List[Equality],
    elements: List[OnePortElement],
    time_symbol: Symbol,
):
    for one_port_element in elements:
        if (
            one_port_element.bond is None
            or one_port_element.bond.effort_symbol is None
            or one_port_element.bond.flow_symbol is None
        ):
            continue
        other_equations += one_port_element.equations(
            one_port_element.bond.effort_symbol,
            one_port_element.bond.flow_symbol,
            time_symbol,
        )
        if isinstance(one_port_element, HasStateEquations):
            state_eqs = one_port_element.integrated_state_equations(
                one_port_element.bond.effort_symbol,
                one_port_element.bond.flow_symbol,
                time_symbol,
            )
            for state_eq in state_eqs:
                if state_eq[0] in state_equations:
                    raise Exception(
                        f"Duplicate state symbol encountered: {state_eq[0]}"
                    )
                state_variables[state_counter] = state_eq[0]
                state_equations[state_eq[0]] = state_eq[1]
                state_counter += 1


def _populate_two_port_equations(
    other_equations: List[Equality], elements: List[TwoPortElement], time_symbol: Symbol
):
    for two_port_element in elements:
        if (
            two_port_element.bond_1 is None
            or two_port_element.bond_2 is None
            or two_port_element.bond_1.effort_symbol is None
            or two_port_element.bond_1.flow_symbol is None
            or two_port_element.bond_2.effort_symbol is None
            or two_port_element.bond_2.flow_symbol is None
        ):
            continue
        other_equations += two_port_element.equations(
            two_port_element.bond_1.effort_symbol,
            two_port_element.bond_2.effort_symbol,
            two_port_element.bond_1.flow_symbol,
            two_port_element.bond_2.flow_symbol,
            time_symbol,
        )


def _populate_junction_equations(
    other_equations: List[Equality], junctions: List[Junction], time_symbol: Symbol
):
    for junction in junctions:
        if isinstance(junction, JunctionEqualEffort):
            if (
                junction.effort_in_bond is None
                or junction.effort_in_bond.flow_symbol is None
            ):
                continue
            # Add new equation for setting effort-in bond's flow symbol equal to the rest of the flows.
            new_eq: Equality = Equality(
                junction.effort_in_bond.flow_symbol(time_symbol), 0  # type: ignore
            )
            if not isinstance(new_eq, Equality):
                raise Exception("Invalid equation at junction")
            for bond in junction.bonds:
                if bond.flow_symbol is None:
                    continue
                if bond is junction.effort_in_bond:
                    continue
                if junction == bond.node_to:
                    new_eq = Equality(
                        new_eq.lhs, new_eq.rhs + bond.flow_symbol(time_symbol)  # type: ignore
                    )
                else:
                    new_eq = Equality(
                        new_eq.lhs, new_eq.rhs - bond.flow_symbol(time_symbol)  # type: ignore
                    )
            if junction.effort_in_bond.node_to == junction:
                new_eq = Equality(new_eq.lhs, -new_eq.rhs)  # type: ignore
            other_equations.append(new_eq)
        elif isinstance(junction, JunctionEqualFlow):
            if (
                junction.effort_out_bond is None
                or junction.effort_out_bond.flow_symbol is None
            ):
                continue
            # Add new equation for setting effort-out bond's effort symbol equal to the rest of the efforts
            new_eq = Equality(
                junction.effort_out_bond.effort_symbol(time_symbol), 0  # type: ignore
            )
            for bond in junction.bonds:
                if bond is junction.effort_out_bond:
                    continue
                if junction == bond.node_to:
                    new_eq = Equality(
                        new_eq.lhs,
                        new_eq.rhs + bond.effort_symbol(time_symbol),  # type: ignore
                    )
                else:
                    new_eq = Equality(
                        new_eq.lhs,
                        new_eq.rhs - bond.effort_symbol(time_symbol),  # type: ignore
                    )
            if junction.effort_out_bond.node_to == junction:
                new_eq = Equality(new_eq.lhs, -new_eq.rhs)  # type: ignore
            other_equations.append(new_eq)


def _substitute_junction_equations(
    junctions: List[Junction],
    state_equations: Dict[Symbol, Expr],
    other_equations: List[Equality],
    time_symbol: Symbol,
):
    substitutions_made = True
    while substitutions_made:
        substitutions_made = False
        for junction in junctions:
            substitutions = dict()
            if isinstance(junction, JunctionEqualEffort):
                # Substitute all effort symbols with the effort-in bond's effort symbol.
                for bond in junction.bonds:
                    if (
                        bond is junction.effort_in_bond
                        or bond.effort_symbol is None
                        or junction.effort_in_bond is None
                        or junction.effort_in_bond.effort_symbol is None
                    ):
                        continue
                    substitutions[
                        bond.effort_symbol(time_symbol)  # type: ignore
                    ] = junction.effort_in_bond.effort_symbol(
                        time_symbol
                    )  # type: ignore
            elif isinstance(junction, JunctionEqualFlow):
                # Substitute all flow symbols with the effort-out bond's flow symbol
                for bond in junction.bonds:
                    if (
                        bond is junction.effort_out_bond
                        or bond.flow_symbol is None
                        or junction.effort_out_bond is None
                        or junction.effort_out_bond.flow_symbol is None
                    ):
                        continue
                    if bond is junction.effort_out_bond:
                        continue
                    substitutions[
                        bond.flow_symbol(time_symbol)  # type: ignore
                    ] = junction.effort_out_bond.flow_symbol(
                        time_symbol
                    )  # type: ignore
            for index, eq in enumerate(other_equations):
                before_eq = other_equations[index]
                replaced_eq = Equality(eq.lhs, eq.rhs.xreplace(substitutions))
                if isinstance(replaced_eq, Equality):
                    other_equations[index] = replaced_eq
                if other_equations[index] != before_eq:
                    substitutions_made = True
            for key, val in state_equations.items():
                before_state_eq = state_equations[key]
                replacement = val.xreplace(substitutions)
                if isinstance(replacement, Expr):
                    state_equations[key] = replacement
                if state_equations[key] != before_state_eq:
                    substitutions_made = True


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
            if bond.node_to is None or bond.node_from is None:
                continue
            failed = False
            if isinstance(bond.node_to, OnePortElement):
                if (
                    bond.effort_in_at_to is True
                    and bond.node_to == Causality.PreferEffortOut
                ) or (
                    bond.effort_in_at_to is False
                    and bond.node_to == Causality.PreferEffortIn
                ):
                    failed = True
            if isinstance(bond.node_from, OnePortElement):
                if (
                    bond.effort_in_at_to is True
                    and bond.node_from == Causality.PreferEffortIn
                ) or (
                    bond.effort_in_at_to is False
                    and bond.node_from == Causality.PreferEffortOut
                ):
                    failed = True
            if failed:
                logging.warning(
                    f"Bond from {bond.node_from.name} to {bond.node_to.name} has non-preferred causality"
                )
                success = False
        return success

    def add(self, bond: Bond):
        if isinstance(bond.node_from, OnePortElement):
            if bond.node_from in self._elements:
                raise Exception(
                    f"OnePortElement {bond.node_from} can only be bonded once!"
                )
            if bond.node_from not in self._elements:
                self._elements.append(bond.node_from)
            bond.node_from.bond = bond
            self._parameters.update(bond.node_from.parameter_symbols())
        elif isinstance(bond.node_from, Junction):
            if bond.node_from not in self._junctions:
                self._junctions.append(bond.node_from)
            if bond not in bond.node_from.bonds:
                bond.node_from.bonds.append(bond)
        elif isinstance(bond.node_from, TwoPortElement):
            if bond.node_from not in self._two_port_elements:
                self._two_port_elements.append(bond.node_from)
            bond.node_from.bond_2 = bond
            self._parameters.update(bond.node_from.parameter_symbols())

        if isinstance(bond.node_to, OnePortElement):
            if bond.node_to in self._elements:
                raise Exception(
                    f"OnePortElement {bond.node_to} can only be bonded once!"
                )
            if bond.node_to not in self._elements:
                self._elements.append(bond.node_to)
            bond.node_to.bond = bond
            self._parameters.update(bond.node_to.parameter_symbols())
        elif isinstance(bond.node_to, Junction):
            if bond.node_to not in self._junctions:
                self._junctions.append(bond.node_to)
            if bond not in bond.node_to.bonds:
                bond.node_to.bonds.append(bond)
        elif isinstance(bond.node_to, TwoPortElement):
            if bond.node_to not in self._two_port_elements:
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
                isinstance(bond.node_to, OnePortElement)
                and bond.node_to.causality_policy() == Causality.FixedEffortIn
                or isinstance(bond.node_from, OnePortElement)
                and bond.node_from.causality_policy() == Causality.FixedEffortOut
            ):
                bond.effort_in_at_to = True
                logging.debug(
                    f"Set fixed effort-in causality at {bond.node_to} (vs {bond.node_from})"
                )
            elif (
                isinstance(bond.node_to, OnePortElement)
                and bond.node_to.causality_policy() == Causality.FixedEffortOut
                or isinstance(bond.node_from, OnePortElement)
                and bond.node_from.causality_policy() == Causality.FixedEffortIn
            ):
                bond.effort_in_at_to = False
                logging.debug(
                    f"Set fixed effort-in causality at {bond.node_from} (vs {bond.node_to})"
                )

    def try_assign_constraint_causalities(self) -> bool:
        something_happened = False
        for junction in self._junctions:
            if junction.assign_constraint_causality():
                something_happened = True
        for element in self._two_port_elements:
            if hasattr(element, "assign_constraint_causality"):
                if element.assign_constraint_causality():
                    something_happened = True
        return something_happened

    def try_assign_preferred_causality(self) -> bool:
        something_happened = False
        for element in self._elements:
            if element.assign_preferred_causality():
                something_happened = True
                # Only assign causality for one element at a time
                break
        return something_happened

    def try_assign_arbitrary_causality(self) -> bool:
        something_happened = False
        for element in self._elements:
            if element.assign_arbitrary_causality():
                something_happened = True
                # Only assign causality for one element at a time
                break
        return something_happened

    def assign_causalities(self) -> None:
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

    def get_state_equations(self) -> Dict[Symbol, Expr]:
        if self._state < _BG_STATE_CAUSALITIES_DONE:
            self.assign_causalities()

        state_equations: Dict[Symbol, Expr] = dict()
        state_variables: Dict[int, Symbol] = dict()
        state_counter = 1
        other_equations: List[Equality] = []
        logging.debug("Formulating equations for one-port elements...")
        _populate_one_port_equations(
            state_equations,
            state_variables,
            state_counter,
            other_equations,
            self._elements,
            self._time_symbol,
        )

        logging.debug("Formulating equations for two-port elements...")
        _populate_two_port_equations(
            other_equations, self._two_port_elements, self._time_symbol
        )

        logging.debug("Formulating equations for junctions...")
        _populate_junction_equations(
            other_equations, self._junctions, self._time_symbol
        )

        logging.debug("Substituting in junction equations...")
        _substitute_junction_equations(
            self._junctions, state_equations, other_equations, self._time_symbol
        )

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

        logging.debug("Generating differential equations...")
        diff_eq_sys: Dict[Symbol, Expr] = dict()
        for var, rhs in state_equations.items():
            diff_eq = rhs.diff(self._time_symbol)

            before = None
            while before != diff_eq:
                before = diff_eq
                diff_eq = diff_eq.xreplace(ordered_equations)

            if isinstance(diff_eq, Expr):
                diff_eq_sys[var] = diff_eq

        return diff_eq_sys

    def get_nodes(self) -> List[Node]:
        import itertools

        node_list = []
        for element in itertools.chain(
            self._elements, self._junctions, self._two_port_elements
        ):
            node_list.append(element)
        return node_list
