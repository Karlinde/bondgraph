from sympy import Symbol
from bondgraph.core import BondGraph, Bond
from bondgraph.common import Node

import graphviz

from typing import Dict

def gen_graphviz(bond_graph: BondGraph) -> graphviz.Digraph:
    g = graphviz.Digraph(
        node_attr={
            'shape': 'none'
        },
        edge_attr={
            'dir': 'both'
        })

    bond_graph.assign_causalities()
    
    # Map from bondgraph node to graphviz node name
    node_map: Dict[object, str] = dict()
    
    n: Node
    for n in bond_graph.get_nodes():
        label = n.visualization_symbol
        if hasattr(n, 'symbol'):
            label += f" : {n.symbol}"
        g.node(n.name, label)    
        node_map[n] =  n.name
    
    b: Bond
    for b in bond_graph._bonds:
        if b.effort_in_at_to == True:
            arrowhead_attr = 'teelvee'
            arrowtail_attr = 'none'
        elif b.effort_in_at_to == False:
            arrowhead_attr = 'lvee'
            arrowtail_attr = 'tee'
        else:
            arrowhead_attr = 'lvee'
            arrowtail_attr = 'none'
        g.edge(node_map[b.node_from], node_map[b.node_to], arrowhead=arrowhead_attr, arrowtail=arrowtail_attr, len='2')
    
    return g
