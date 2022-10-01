from enum import Enum


class Causality(Enum):
    Indifferent = 0
    PreferEffortIn = 1
    PreferEffortOut = 2
    FixedEffortIn = 3
    FixedEffortOut = 4
