from sympy import Symbol, Equality, Expr, Function
from typing import List, Set, Tuple
from bondgraph.common import Causality, Node, Bond, HasStateEquations

import logging


class OnePortElement(Node):
    def __init__(self, name: str):
        super().__init__(name)
        self.bond: Bond | None = None

    def equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Equality]:
        raise NotImplementedError()

    @staticmethod
    def causality_policy() -> Causality:
        raise NotImplementedError()

    def assign_preferred_causality(self) -> bool:
        return False

    def assign_arbitrary_causality(self) -> bool:
        return False

    def __str__(self) -> str:
        return self.name

    def parameter_symbols(self) -> Set[Symbol]:
        raise NotImplementedError()


class TwoPortElement(Node):
    def __init__(self, name: str):
        super().__init__(name)
        self.bond_1: Bond | None = None
        self.bond_2: Bond | None = None

    def equations(
        self,
        effort_1: Function,
        effort_2: Function,
        flow_1: Function,
        flow_2: Function,
        time: Symbol,
    ) -> List[Equality]:
        raise NotImplementedError()

    def assign_constraint_causality(self) -> bool:
        return False

    def __str__(self) -> str:
        return self.name

    def parameter_symbols(self) -> Set[Symbol]:
        raise NotImplementedError()


class Element_R(OnePortElement):
    def __init__(self, name: str, symbol: Symbol):
        super().__init__(name)
        self.symbol = symbol

    def equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Equality]:
        if self.bond is None:
            return []
        if self.bond.effort_in_at_to is True and self.bond.node_to == self:
            return [Equality(flow(time), effort(time) / self.symbol)]  # type: ignore
        else:
            return [Equality(effort(time), self.symbol * flow(time))]  # type: ignore

    @staticmethod
    def causality_policy() -> Causality:
        return Causality.Indifferent

    def assign_arbitrary_causality(self) -> bool:
        if self.bond is None:
            return False
        if not self.bond.has_causality_set():
            self.bond.effort_in_at_to = True
            logging.debug(
                f"Set indifferent effort-in causality at {self.bond.node_to} (vs {self.bond.node_from})"
            )
            return True
        return False

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.symbol}

    def visualization_label(self) -> str:
        return f"R: {self.symbol}"


class Element_C(OnePortElement, HasStateEquations):
    def __init__(self, name: str, compliance: Symbol, displacement: Symbol):
        super().__init__(name)
        self._compliance = compliance
        self._displacement = displacement

    def equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Equality]:
        return [Equality(effort(time), self._displacement / self._compliance)]  # type: ignore

    def integrated_state_equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Tuple[Symbol, Expr]]:
        return [(self._displacement, flow(time).integrate(time))]  # type: ignore

    @staticmethod
    def causality_policy():
        return Causality.PreferEffortOut

    def assign_preferred_causality(self) -> bool:
        if self.bond is None:
            return False
        if not self.bond.has_causality_set():
            if self.bond.node_from == self:
                self.bond.effort_in_at_to = True
                logging.debug(
                    f"Set preferred effort-out causality at {self.bond.node_from} (vs {self.bond.node_to})"
                )
                return True
            elif self.bond.node_to == self:
                self.bond.effort_in_at_to = False
                logging.debug(
                    f"Set preferred effort-out causality at {self.bond.node_to} (vs {self.bond.node_from})"
                )
                return True
        return False

    def parameter_symbols(self) -> Set[Symbol]:
        return {self._compliance}

    def visualization_label(self) -> str:
        return f"C: {self._compliance}"


class Element_I(OnePortElement, HasStateEquations):
    def __init__(self, name: str, inertia: Symbol, momentum: Symbol):
        super().__init__(name)
        self._inertia = inertia
        self._momentum = momentum

    def equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Equality]:
        return [Equality(flow(time), self._momentum / self._inertia)]  # type: ignore

    def integrated_state_equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Tuple[Symbol, Expr]]:
        return [(self._momentum, effort(time).integrate(time))]  # type: ignore

    @staticmethod
    def causality_policy():
        return Causality.PreferEffortIn

    def assign_preferred_causality(self) -> bool:
        if self.bond is None:
            return False
        if not self.bond.has_causality_set():
            if self.bond.node_from == self:
                self.bond.effort_in_at_to = False
                logging.debug(
                    f"Set preferred effort-in causality at {self.bond.node_from} (vs {self.bond.node_to})"
                )
                return True
            elif self.bond.node_to == self:
                self.bond.effort_in_at_to = True
                logging.debug(
                    f"Set preferred effort-in causality at {self.bond.node_to} (vs {self.bond.node_from})"
                )
                return True
        return False

    def parameter_symbols(self) -> Set[Symbol]:
        return {self._inertia}

    def visualization_label(self) -> str:
        return f"I: {self._inertia}"


class Source_effort(OnePortElement):
    def __init__(self, name: str, symbol: Symbol):
        super().__init__(name)
        self.symbol = symbol

    def equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Equality]:
        return [Equality(effort(time), self.symbol)]  # type: ignore

    @staticmethod
    def causality_policy():
        return Causality.FixedEffortOut

    def assign_fixed_causality(self):
        if self.bond is None:
            return False
        if not self.bond.has_causality_set():
            if self.bond.node_from == self:
                self.bond.effort_in_at_to = True
                logging.debug(
                    f"Set fixed effort-out causality at {self.bond.node_from} (vs {self.bond.node_to})"
                )
                return True
            elif self.bond.node_to == self:
                self.bond.effort_in_at_to = False
                logging.debug(
                    f"Set fixed effort-out causality at {self.bond.node_to} (vs {self.bond.node_from})"
                )
                return True
        return False

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.symbol}

    def visualization_label(self) -> str:
        return f"Se: {self.symbol}"


class Source_flow(OnePortElement):
    def __init__(self, name: str, symbol: Symbol):
        super().__init__(name)
        self.symbol = symbol

    def equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Equality]:
        return [Equality(flow(time), self.symbol)]  # type: ignore

    @staticmethod
    def causality_policy():
        return Causality.FixedEffortIn

    def assign_fixed_causality(self):
        if self.bond is None:
            return False
        if not self.bond.has_causality_set():
            if self.bond.node_from == self:
                self.bond.effort_in_at_to = False
                logging.debug(
                    f"Set fixed effort-in causality at {self.bond.node_from} (vs {self.bond.node_to})"
                )
                return True
            elif self.bond.node_to == self:
                self.bond.effort_in_at_to = True
                logging.debug(
                    f"Set fixed effort-in causality at {self.bond.node_to} (vs {self.bond.node_from})"
                )
                return True
        return False

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.symbol}

    def visualization_label(self) -> str:
        return f"Sf: {self.symbol}"


class Transformer(TwoPortElement):
    def __init__(self, name: str, ratio: Symbol):
        super().__init__(name)
        self.ratio = ratio

    def equations(
        self,
        effort_1: Function,
        effort_2: Function,
        flow_1: Function,
        flow_2: Function,
        time: Symbol,
    ) -> List[Equality]:
        if self.bond_1 is None or self.bond_2 is None:
            raise Exception("Transformer is not fully connected")
        if self.bond_1.effort_in_at_to:
            return [
                Equality(flow_1(time), flow_2(time) / self.ratio),  # type: ignore
                Equality(effort_2(time), effort_1(time) / self.ratio),  # type: ignore
            ]
        elif not self.bond_1.effort_in_at_to:
            return [
                Equality(flow_2(time), flow_1(time) * self.ratio),  # type: ignore
                Equality(effort_1(time), effort_2(time) * self.ratio),  # type: ignore
            ]
        else:
            raise Exception(f"Invalid causality at transformer {self.name}")

    def assign_constraint_causality(self):
        if self.bond_1 is None or self.bond_2 is None:
            return False
        if self.bond_1.has_causality_set() and not self.bond_2.has_causality_set():
            self.bond_2.effort_in_at_to = self.bond_1.effort_in_at_to
            msg_dir = "effort-in" if self.bond_2.effort_in_at_to else "effort-out"
            logging.debug(
                f"Set constrained {msg_dir} causality at {self.bond_2.node_to} "
                + f"(vs {self.bond_2.node_from}) due to transformer {self}"
            )
            return True
        elif self.bond_2.has_causality_set() and not self.bond_1.has_causality_set():
            self.bond_1.effort_in_at_to = self.bond_2.effort_in_at_to
            msg_dir = "effort-in" if self.bond_1.effort_in_at_to else "effort-out"
            logging.debug(
                f"Set constrained {msg_dir} causality at {self.bond_1.node_to} "
                + f"(vs {self.bond_1.node_from}) due to transformer {self}"
            )
            return True
        return False

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.ratio}

    def visualization_label(self) -> str:
        return f"TF\n{self.ratio}"


class Gyrator(TwoPortElement):
    def __init__(self, name: str, ratio: Symbol):
        super().__init__(name)
        self.ratio = ratio

    def equations(
        self,
        effort_1: Function,
        effort_2: Function,
        flow_1: Function,
        flow_2: Function,
        time: Symbol,
    ):
        if self.bond_1 is None or self.bond_2 is None:
            raise Exception("Gyrator is not fully connected")

        if self.bond_1.effort_in_at_to:
            return [
                Equality(flow_1(time), effort_2(time) / self.ratio),  # type: ignore
                Equality(flow_2(time), effort_1(time) / self.ratio),  # type: ignore
            ]
        elif not self.bond_1.effort_in_at_to:
            return [
                Equality(effort_2(time), flow_1(time) * self.ratio),  # type: ignore
                Equality(effort_1(time), flow_2(time) * self.ratio),  # type: ignore
            ]

    def assign_constraint_causality(self):
        if self.bond_1 is None or self.bond_2 is None:
            return False
        if self.bond_1.has_causality_set() and not self.bond_2.has_causality_set():
            self.bond_2.effort_in_at_to = not self.bond_1.effort_in_at_to
            msg_dir = "effort-in" if self.bond_2.effort_in_at_to else "effort-out"
            logging.debug(
                f"Set constrained {msg_dir} causality at {self.bond_2.node_to} "
                + f"(vs {self.bond_2.node_from}) due to gyrator {self}"
            )
            return True
        elif self.bond_2.has_causality_set() and not self.bond_1.has_causality_set():
            self.bond_1.effort_in_at_to = not self.bond_2.effort_in_at_to
            msg_dir = "effort-in" if self.bond_1.effort_in_at_to else "effort-out"
            logging.debug(
                f"Set constrained {msg_dir} causality at {self.bond_1.node_to} "
                + f"(vs {self.bond_1.node_from}) due to gyrator {self}"
            )
            return True
        return False

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.ratio}

    def visualization_label(self) -> str:
        return f"GY\n{self.ratio}"
