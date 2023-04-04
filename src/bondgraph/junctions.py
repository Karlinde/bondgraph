from bondgraph.common import Node, Bond, Causality, count_bonds_with_causalities_set
from typing import List
import logging


class Junction(Node):
    def __init__(self, name: str):
        super().__init__(name)
        self.bonds: List[Bond] = []

    @staticmethod
    def causality_policy() -> Causality:
        raise NotImplementedError()

    def assign_constraint_causality(self) -> bool:
        raise NotImplementedError()

    def __str__(self) -> str:
        return self.name


class JunctionEqualEffort(Junction):
    def __init__(self, name: str):
        super().__init__(name)
        self.effort_in_bond: Bond | None = None

    def assign_constraint_causality(self) -> bool:
        num_bonds_with_causality = count_bonds_with_causalities_set(self.bonds)
        if num_bonds_with_causality == len(self.bonds):
            return False

        # First, check if there is already a dominant effort-in causality
        for bond in self.bonds:
            if (bond.node_to == self and bond.effort_in_at_to == True) or (
                bond.node_from == self and bond.effort_in_at_to == False
            ):
                self.effort_in_bond = bond
                break

        something_happened = False

        if self.effort_in_bond is not None:
            # There is a dominant causality, set all other bonds as effort-out
            for bond in self.bonds:
                if bond == self.effort_in_bond:
                    continue
                if bond.node_to == self:
                    bond.effort_in_at_to = False
                    logging.debug(
                        f"Set constraint effort-out causality at {bond.node_to} (vs {bond.node_from}) due to equal-effort junction {self}"
                    )
                    something_happened = True
                elif bond.node_from == self:
                    bond.effort_in_at_to = True
                    logging.debug(
                        f"Set constraint effort-out causality at {bond.node_from} (vs {bond.node_to}) due to equal-effort junction {self}"
                    )
                    something_happened = True

        elif num_bonds_with_causality == (len(self.bonds) - 1):
            # No dominant causality but all except for one bond already have causalities,
            # the remaining one must be effort-in
            for bond in self.bonds:
                if bond.effort_in_at_to is None:
                    if bond.node_to == self:
                        bond.effort_in_at_to = True
                        self.effort_in_bond = bond
                        logging.debug(
                            f"Set constraint effort-in causality at {bond.node_to} (vs {bond.node_from}) due to equal-effort junction {self}"
                        )
                        something_happened = True
                    elif bond.node_from == self:
                        bond.effort_in_at_to = False
                        self.effort_in_bond = bond
                        logging.debug(
                            f"Set constraint effort-in causality at {bond.node_from} (vs {bond.node_to}) due to equal-effort junction {self}"
                        )
                        something_happened = True
                    break

        return something_happened

    def visualization_label(self) -> str:
        return "0"


class JunctionEqualFlow(Junction):
    def __init__(self, name: str):
        super().__init__(name)
        self.effort_out_bond: Bond | None = None

    def assign_constraint_causality(self) -> bool:
        num_bonds_with_causality = count_bonds_with_causalities_set(self.bonds)
        if num_bonds_with_causality == len(self.bonds):
            return False

        # First, check if there is already a dominant effort-out causality
        for bond in self.bonds:
            if (bond.node_to == self and bond.effort_in_at_to == False) or (
                bond.node_from == self and bond.effort_in_at_to == True
            ):
                self.effort_out_bond = bond
                break

        something_happened = False

        if self.effort_out_bond is not None:
            # There is a dominant causality, set all other bonds as effort-in
            for bond in self.bonds:
                if bond == self.effort_out_bond:
                    continue
                if bond.node_to == self:
                    bond.effort_in_at_to = True
                    logging.debug(
                        f"Set constraint effort-in causality at {bond.node_to} (vs {bond.node_from}) due to equal-effort junction {self}"
                    )
                    something_happened = True
                elif bond.node_from == self:
                    bond.effort_in_at_to = False
                    logging.debug(
                        f"Set constraint effort-in causality at {bond.node_from} (vs {bond.node_to}) due to equal-effort junction {self}"
                    )
                    something_happened = True

        elif num_bonds_with_causality == (len(self.bonds) - 1):
            # No dominant causality but all except for one bond already have causalities,
            # the remaining one must be effort-out
            for bond in self.bonds:
                if bond.effort_in_at_to is None:
                    if bond.node_to == self:
                        bond.effort_in_at_to = False
                        self.effort_out_bond = bond
                        logging.debug(
                            f"Set constraint effort-out causality at {bond.node_to} (vs {bond.node_from}) due to equal-effort junction {self}"
                        )
                        something_happened = True
                    elif bond.node_from == self:
                        bond.effort_in_at_to = True
                        self.effort_out_bond = bond
                        logging.debug(
                            f"Set constraint effort-out causality at {bond.node_from} (vs {bond.node_to}) due to equal-effort junction {self}"
                        )
                        something_happened = True
                    break

        return something_happened

    def visualization_label(self) -> str:
        return "1"
