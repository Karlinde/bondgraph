from abc import ABC, abstractmethod
from enum import Enum
from typing import List, Tuple

from sympy import Expr, Symbol


class Node:
    def __init__(self, name: str):
        self.name = name

    def visualization_label(self) -> str:
        raise NotImplementedError()


class Bond:
    def __init__(self, node_from: Node, node_to: Node):
        self.node_from: Node | None = node_from
        self.node_to: Node | None = node_to
        self.num: int | None = None
        self.effort_in_at_to: bool | None = None
        self.flow_symbol: Symbol | None = None
        self.effort_symbol: Symbol | None = None

    def has_causality_set(self) -> bool:
        return self.effort_in_at_to is not None


class Causality(Enum):
    Indifferent = 0
    PreferEffortIn = 1
    PreferEffortOut = 2
    FixedEffortIn = 3
    FixedEffortOut = 4


class HasStateEquations(ABC):
    @abstractmethod
    def state_equations(
        self, effort: Symbol, flow: Symbol
    ) -> List[Tuple[Symbol, Expr]]:
        """
        Return the governing equations for this element as a list of tuples
        where the first item of each tuple is the symbol of a state variable and
        the second item is the right-hand side of the corresponding state-space
        equation.
        """


def count_bonds_with_causalities_set(bonds: List[Bond]) -> int:
    count = 0
    for b in bonds:
        if b.has_causality_set():
            count += 1
    return count


class AlgebraicLoopError(Exception):
    pass
