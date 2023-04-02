"""
This example generates the following simple bond graph:

Se ---> 1 ---> I
        |
        v
        R
"""
from bondgraph.core import Bond, BondGraph
from bondgraph.junctions import JunctionEqualFlow
from bondgraph.elements.basic import Element_R, Element_I, Source_effort
from sympy import Symbol

force = Source_effort("force", Symbol("F"))
friction = Element_R("friction", Symbol("k_f"))
inertia = Element_I("inertia", Symbol("m"), Symbol("p"))
mass_object = JunctionEqualFlow("mass_object")

graph = BondGraph(Symbol("t"))
graph.add(Bond(force, mass_object))
graph.add(Bond(mass_object, friction))
graph.add(Bond(mass_object, inertia))

state_equations, state_vars = graph.get_state_equations()

# Print the dictionary of state numbers and the right hand side of their state equations:
print(state_equations)
