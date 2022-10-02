from bondgraph.core import Causality
from bondgraph.elements.basic import OnePortElement
from sympy import Symbol, Equality, Function, tanh
from typing import Set


class Element_ReliefValve(OnePortElement):
    def __init__(self, name: str, ro: Symbol, rc: Symbol, k: Symbol, d: Symbol):
        super().__init__(name)
        self.ro = ro
        self.rc = rc
        self.d = d
        self.k = k

    def equations(self, effort: Function, flow: Function, time: Symbol):
        return Equality(
            flow(time),
            effort(time)
            / (
                self.rc * (0.5 - 0.5 * tanh(self.k * (effort(time) - self.d)))
                + self.ro * (0.5 + 0.5 * tanh(self.k * (effort(time) - self.d)))
            ),
        )

    @staticmethod
    def causality_policy():
        return Causality.FixedEffortIn

    def parameter_symbols(self) -> Set[Symbol]:
        return {self.ro, self.rc, self.d, self.k}
