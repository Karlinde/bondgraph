from enum import Enum
from typing import List


class Bond:
    def __init__(self, node_from, node_to):
        self.node_from = node_from
        self.node_to = node_to
        self.num = None
        self.effort_in_at_to = None
        self.flow_symbol = None
        self.effort_symbol = None

    def has_causality_set(self) -> bool:
        return self.effort_in_at_to is not None


class Node:
    def __init__(self, name: str, visualization_symbol: str):
        self.name = name
        self.visualization_symbol = visualization_symbol


class Causality(Enum):
    Indifferent = 0
    PreferEffortIn = 1
    PreferEffortOut = 2
    FixedEffortIn = 3
    FixedEffortOut = 4


def count_bonds_with_causalities_set(bonds: List[Bond]) -> int:
    count = 0
    for b in bonds:
        if b.has_causality_set():
            count += 1
    return count
