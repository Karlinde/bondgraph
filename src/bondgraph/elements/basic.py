from sympy import Function, Symbol, Equality, Expr
from typing import List, Set, Tuple
from bondgraph.common import Causality, Node, Bond

import logging


class OnePortElement(Node):
    def __init__(self, name, visualization_symbol):
        super().__init__(name, visualization_symbol)
        self.bond: Bond = None

    def equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Equality]:
        pass

    @staticmethod
    def causality_policy():
        pass

    def __str__(self):
        return self.name

    def parameter_symbols(self) -> Set[Symbol]:
        pass


class TwoPortElement(Node):
    def __init__(self, name, visualization_symbol):
        super().__init__(name, visualization_symbol)
        self.bond_1: Bond = None
        self.bond_2: Bond = None

    def equations(
        self,
        effort_1: Function,
        effort_2: Function,
        flow_1: Function,
        flow_2: Function,
        time: Symbol,
    ) -> Tuple[Expr, Expr]:
        pass

    @staticmethod
    def causality_policy():
        pass

    def __str__(self):
        return self.name

    def parameter_symbols(self) -> Set[Symbol]:
        pass


class Element_R(OnePortElement):
    def __init__(self, name: str, symbol: Symbol):
        super().__init__(name, "R")
        self.symbol = symbol

    def equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Equality]:
        if self.bond.effort_in_at_to == True and self.bond.node_to == self:
            return [Equality(flow(time), effort(time) / self.symbol)]
        else:
            return [Equality(effort(time), self.symbol * flow(time))]

    @staticmethod
    def causality_policy():
        return Causality.Indifferent

    def assign_arbitrary_causality(self) -> bool:
        if not self.bond.has_causality_set():
            self.bond.effort_in_at_to = True
            logging.debug(
                f"Set indifferent effort-in causality at {self.bond.node_to} (vs {self.bond.node_from})"
            )
            return True
        return False

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.symbol}


class Element_C(OnePortElement):
    def __init__(self, name: str, compliance: Symbol, displacement: Symbol):
        super().__init__(name, "C")
        self._compliance = compliance
        self._displacement = displacement

    def equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Equality]:
        return [Equality(effort(time), self._displacement / self._compliance)]

    def state_equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Tuple[Symbol, Expr]]:
        return [(self._displacement, flow(time).integrate(time))]

    @staticmethod
    def causality_policy():
        return Causality.PreferEffortOut

    def assign_preferred_causality(self) -> bool:
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


class Element_I(OnePortElement):
    def __init__(self, name: str, inertia: Symbol, momentum: Symbol):
        super().__init__(name, "I")
        self._inertia = inertia
        self._momentum = momentum

    def equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Equality]:
        return [Equality(flow(time), self._momentum / self._inertia)]

    def state_equations(
        self, effort: Function, flow: Function, time: Symbol
    ) -> List[Tuple[Symbol, Expr]]:
        return [(self._momentum, effort(time).integrate(time))]

    @staticmethod
    def causality_policy():
        return Causality.PreferEffortIn

    def assign_preferred_causality(self) -> bool:
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


class Source_effort(OnePortElement):
    def __init__(self, name: str, symbol: Symbol):
        super().__init__(name, "Se")
        self.symbol = symbol

    def equations(self, effort: Function, flow: Function, time: Symbol):
        return [Equality(effort(time), self.symbol)]

    @staticmethod
    def causality_policy():
        return Causality.FixedEffortOut

    def assign_fixed_causality(self):
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


class Source_flow(OnePortElement):
    def __init__(self, name: str, symbol: Symbol):
        super().__init__(name, "Sf")
        self.symbol = symbol

    def equations(self, effort: Function, flow: Function, time: Symbol):
        return [Equality(flow(time), self.symbol)]

    @staticmethod
    def causality_policy():
        return Causality.FixedEffortIn

    def assign_fixed_causality(self):
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


class Transformer(TwoPortElement):
    def __init__(self, name: str, ratio: Symbol):
        super().__init__(name, "TF")
        self.ratio = ratio

    def equations(
        self,
        effort_1: Function,
        effort_2: Function,
        flow_1: Function,
        flow_2: Function,
        time: Symbol,
    ):
        if self.bond_1.effort_in_at_to:
            return [
                Equality(flow_1(time), flow_2(time) / self.ratio),
                Equality(effort_2(time), effort_1(time) / self.ratio),
            ]
        elif not self.bond_1.effort_in_at_to:
            return [
                Equality(flow_2(time), flow_1(time) * self.ratio),
                Equality(effort_1(time), effort_2(time) * self.ratio),
            ]
        else:
            raise Exception(f"Invalid causality at transformer {self.name}")

    def assign_constraint_causality(self):
        if self.bond_1.has_causality_set() and not self.bond_2.has_causality_set():
            self.bond_2.effort_in_at_to = self.bond_1.effort_in_at_to
            msg_dir = "effort-in" if self.bond_2.effort_in_at_to else "effort-out"
            logging.debug(
                f"Set constrained {msg_dir} causality at {self.bond_2.node_to} (vs {self.bond_2.node_from}) due to transformer {self}"
            )
            return True
        elif self.bond_2.has_causality_set() and not self.bond_1.has_causality_set():
            self.bond_1.effort_in_at_to = self.bond_2.effort_in_at_to
            msg_dir = "effort-in" if self.bond_1.effort_in_at_to else "effort-out"
            logging.debug(
                f"Set constrained {msg_dir} causality at {self.bond_1.node_to} (vs {self.bond_1.node_from}) due to transformer {self}"
            )
            return True
        return False

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.ratio}


class Gyrator(TwoPortElement):
    def __init__(self, name: str, ratio: Symbol):
        super().__init__(name, "GY")
        self.ratio = ratio

    def equations(
        self,
        effort_1: Function,
        effort_2: Function,
        flow_1: Function,
        flow_2: Function,
        time: Symbol,
    ):
        if self.bond_1.effort_in_at_to:
            return [
                Equality(flow_1(time), effort_2(time) / self.ratio),
                Equality(flow_2(time), effort_1(time) / self.ratio),
            ]
        elif not self.bond_1.effort_in_at_to:
            return [
                Equality(effort_2(time), flow_1(time) * self.ratio),
                Equality(effort_1(time), flow_2(time) * self.ratio),
            ]

    def assign_constraint_causality(self):
        if self.bond_1.has_causality_set() and not self.bond_2.has_causality_set():
            self.bond_2.effort_in_at_to = not self.bond_1.effort_in_at_to
            msg_dir = "effort-in" if self.bond_2.effort_in_at_to else "effort-out"
            logging.debug(
                f"Set constrained {msg_dir} causality at {self.bond_2.node_to} (vs {self.bond_2.node_from}) due to gyrator {self}"
            )
            return True
        elif self.bond_2.has_causality_set() and not self.bond_1.has_causality_set():
            self.bond_1.effort_in_at_to = not self.bond_2.effort_in_at_to
            msg_dir = "effort-in" if self.bond_1.effort_in_at_to else "effort-out"
            logging.debug(
                f"Set constrained {msg_dir} causality at {self.bond_1.node_to} (vs {self.bond_1.node_from}) due to gyrator {self}"
            )
            return True
        return False

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.ratio}
