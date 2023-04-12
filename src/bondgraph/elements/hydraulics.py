from bondgraph.common import Causality
from bondgraph.elements.basic import OnePortElement
from sympy import Symbol, Equality, tanh
from typing import Set


class Element_ReliefValve(OnePortElement):
    def __init__(self, name: str, ro: Symbol, rc: Symbol, k: Symbol, d: Symbol):
        super().__init__(name)
        self.ro = ro
        self.rc = rc
        self.d = d
        self.k = k

    def equations(self, effort: Symbol, flow: Symbol):
        return [
            Equality(
                flow,
                effort
                / (
                    self.rc * (0.5 - 0.5 * tanh(self.k * (effort - self.d)))  # type: ignore
                    + self.ro * (0.5 + 0.5 * tanh(self.k * (effort - self.d)))  # type: ignore
                ),
            )
        ]

    @staticmethod
    def causality_policy():
        return Causality.FixedEffortIn

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.ro, self.rc, self.d, self.k}

    def visualization_label(self) -> str:
        return "Rv"
