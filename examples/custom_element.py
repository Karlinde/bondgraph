"""
This is an example of how to implement a custom one-port element.
"""
from bondgraph.common import Causality
from bondgraph.elements import OnePortElement
from sympy import Symbol, Equality, tanh
from typing import List


class Element_ReliefValve(OnePortElement):
    """
    A hydraulic relief valve, which opens if pressure difference across it is too high.

    Implemented as an R-element with variable resistance based on the logistic
    function, to keep the expressions differentiable while approximating a step
    function.

    Since the element is a variation of an R-element, it's intended to be
    connected to an equal-flow junction representing the path through the relief
    valve.
    """
    def __init__(self, name: str, ro: Symbol, rc: Symbol, k: Symbol, d: Symbol):
        super().__init__(name)
        self.ro = ro  # Equivalent R-element resistance when valve is open
        self.rc = rc  # Equivalent R-element resistance when valve is closed
        self.d = d    # Effort (pressure) at which to open
        self.k = k    # Smoothness of step at the opening point (k-value of logistic function)

    def equations(self, effort: Symbol, flow: Symbol) -> List[Equality]:
        """
        Generate a list of the equations describing the element
        """
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
        """
        Describe which causality this element expects. This will be used when
        validating the bond graph.
        """
        return Causality.FixedEffortIn

    def visualization_label(self) -> str:
        """
        Specify how this element should be visualized in a graphical
        representation of the graph.
        """
        return "Rv"
