from enum import Enum


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
