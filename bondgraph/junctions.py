class Junction:
    def __init__(self, name):
        self.bonds = []
        self.name = name

    @staticmethod
    def causality_policy():
        return None

    def __str__(self):
        return self.name


class JunctionEqualEffort(Junction):
    def __init__(self, name: str):
        super().__init__(name)
        self.effort_in_bond = None


class JunctionEqualFlow(Junction):
    def __init__(self, name: str):
        super().__init__(name)
        self.effort_out_bond = None
