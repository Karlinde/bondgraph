from sympy import Function, Symbol, Equality, Expr
from typing import Set, Tuple
from bondgraph.common import Causality


class OnePortElement:
    def __init__(self, name):
        self.bond = None
        self.name = name

    def equations(self, effort: Function, flow: Function, time: Symbol):
        pass

    @staticmethod
    def causality_policy():
        pass

    def __str__(self):
        return self.name

    def parameter_symbols(self) -> Set[Symbol]:
        pass


class TwoPortElement:
    def __init__(self, name):
        self.bond_1 = None
        self.bond_2 = None
        self.name = name

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
        super().__init__(name)
        self.symbol = symbol

    def equations(self, effort: Function, flow: Function, time: Symbol):
        if self.bond.effort_in_at_to == True and self.bond.node_to == self:
            return Equality(flow(time), effort(time) / self.symbol)
        else:
            return Equality(effort(time), self.symbol * flow(time))

    @staticmethod
    def causality_policy():
        return Causality.Indifferent

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.symbol}


class Element_C(OnePortElement):
    def __init__(self, name: str, symbol: Symbol):
        super().__init__(name)
        self.symbol = symbol

    def equations(self, effort: Function, flow: Function, time: Symbol):
        return Equality(effort(time), 1 / self.symbol * flow(time).integrate(time))

    @staticmethod
    def causality_policy():
        return Causality.PreferEffortOut

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.symbol}


class Element_I(OnePortElement):
    def __init__(self, name: str, symbol: Symbol):
        super().__init__(name)
        self.symbol = symbol

    def equations(self, effort: Function, flow: Function, time: Symbol):
        return Equality(flow(time), 1 / self.symbol * effort(time).integrate(time))

    @staticmethod
    def causality_policy():
        return Causality.PreferEffortIn

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.symbol}


class Source_effort(OnePortElement):
    def __init__(self, name: str, symbol: Symbol):
        super().__init__(name)
        self.symbol = symbol

    def equations(self, effort: Function, flow: Function, time: Symbol):
        return Equality(effort(time), self.symbol)

    @staticmethod
    def causality_policy():
        return Causality.FixedEffortOut

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.symbol}


class Source_flow(OnePortElement):
    def __init__(self, name: str, symbol: Symbol):
        super().__init__(name)
        self.symbol = symbol

    def equations(self, effort: Function, flow: Function, time: Symbol):
        return Equality(flow(time), self.symbol)

    @staticmethod
    def causality_policy():
        return Causality.FixedEffortIn

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.symbol}


class Transformer(TwoPortElement):
    def __init__(self, name: str, ratio: Symbol):
        super().__init__(name)
        self.ratio = ratio
        self.effort_in_bond = None

    def equations(
        self,
        effort_1: Function,
        effort_2: Function,
        flow_1: Function,
        flow_2: Function,
        time: Symbol,
    ):
        if (self.bond_1.effort_in_at_to == True and self.bond_1.node_to == self) or (
            self.bond_1.effort_in_at_to == False and self.bond_1.node_to != self
        ):
            return [
                Equality(flow_1(time), flow_2(time) / self.ratio),
                Equality(effort_2(time), effort_1(time) / self.ratio),
            ]
        elif (self.bond_2.effort_in_at_to == True and self.bond_2.node_to == self) or (
            self.bond_2.effort_in_at_to == False and self.bond_2.node_to != self
        ):
            return [
                Equality(flow_2(time), flow_1(time) * self.ratio),
                Equality(effort_1(time), effort_2(time) * self.ratio),
            ]

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.ratio}
