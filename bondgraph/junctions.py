from bondgraph.common import Node

class Junction(Node):
    def __init__(self, name, visualization_symbol):
        super().__init__(name, visualization_symbol)
        self.bonds = []

    @staticmethod
    def causality_policy():
        return None

    def __str__(self):
        return self.name


class JunctionEqualEffort(Junction):
    def __init__(self, name: str):
        super().__init__(name, '0')
        self.effort_in_bond = None


class JunctionEqualFlow(Junction):
    def __init__(self, name: str):
        super().__init__(name, '1')
        self.effort_out_bond = None
