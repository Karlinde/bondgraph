# bondgraph
A python library for creating [bond graphs](https://en.wikipedia.org/wiki/Bond_graph) of 
physical systems and generating differential equations.

It heavily relies on `sympy` for manipulation and simplification of equations symbolically.

## Installation
Install via pip:
```
pip install bondgraph
```

To install with extra dependencies for visualization using graphviz:
```
pip install bondgraph[visualization]
``` 

## Usage
Create a `BondGraph` object and add `Bond` objects connecting various elements:
```python
# This example generates the following simple bond graph:
#
# Se ---> 1 ---> I
#         |
#         v
#         R

from bondgraph.core import Bond, BondGraph
from bondgraph.junctions import Junction, JunctionEqualFlow
from bondgraph.elements.basic import Element_R, Element_I, Source_effort
from sympy import Symbol

force = Source_effort('force', Symbol('F'))
friction = Element_R('friction', Symbol('k_f'))
inertia = Element_I('inertia', Symbol('m'))
mass_object = JunctionEqualFlow('mass_object')

graph = BondGraph(Symbol('t'))
graph.add(Bond(force, mass_object))
graph.add(Bond(mass_object, friction))
graph.add(Bond(mass_object, inertia))

state_equations, state_vars, state_names = graph.get_state_equations()

# Print the dictionary of state numbers and the right hand side of their state equations: 
print(state_equations)
# {1: (F - inertia_state*k_f)/m}

# If using the optional visualization features, generate and display a graphviz graph:
from bondgraph.visualization import gen_graphviz

output_graph = gen_graphviz(graph)
output.view()
```

## Limitations
- Currently only a few standard elements are implemented.
- Algebraic loops are not handled at all and will result in failure to generate equations.
- Non-integrating (differential) causality for C or I elements is not currently possible.

