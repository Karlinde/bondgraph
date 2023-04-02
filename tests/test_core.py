from bondgraph.core import Bond, BondGraph
from bondgraph.junctions import Junction, JunctionEqualEffort, JunctionEqualFlow
from bondgraph.elements.basic import Element_R, Element_I, Element_C, Transformer, Source_effort, Source_flow

from sympy import Symbol as _

def test_basic_i_element():
    sym_t = _('t')
    sym_F = _('F')
    sym_r = _('r')
    sym_i = _('i')
    sym_p = _('p')

    se = Source_effort('f', sym_F)
    r = Element_R('r', sym_r)
    i = Element_I('i', sym_i, sym_p)
    j = JunctionEqualFlow('j')

    g = BondGraph(sym_t)
    g.add(Bond(se, j))
    g.add(Bond(j, r))
    g.add(Bond(j, i))

    eqs, vars = g.get_state_equations()
    assert eqs[1] == (sym_F - sym_r * sym_p / sym_i)
    assert vars[0] == sym_p

def test_basic_c_element():
    sym_t = _('t')
    sym_F = _('F')
    sym_r = _('r')
    sym_c = _('c')
    sym_q = _('q')

    se = Source_effort('F', sym_F)
    r = Element_R('r', sym_r)
    c = Element_C('c', sym_c, sym_q)
    j = JunctionEqualFlow('j')

    g = BondGraph(sym_t)
    g.add(Bond(se, j))
    g.add(Bond(j, r))
    g.add(Bond(j, c))

    eqs, vars = g.get_state_equations()
    assert eqs[1] == (sym_F - sym_q / sym_c) / sym_r
    assert vars[0] == sym_q

def test_more_complex():
    t = _('t')
    F = _('F')
    v = _('v')
    r = _('r')
    i = _('i')
    c = _('c')
    q = _('q')
    p = _('p')

    e_se = Source_effort('F', F)
    j1 = JunctionEqualFlow('j1')
    e_i = Element_I('i', i, p)
    j2 = JunctionEqualEffort('j2')
    j3 = JunctionEqualFlow('j3')
    e_c = Element_C('c', c, q)
    e_r = Element_R('r', r)
    e_sf = Source_flow('v', v)

    g = BondGraph(t)
    g.add(Bond(e_se, j1))
    g.add(Bond(j1, e_i))
    g.add(Bond(j1, j2))
    g.add(Bond(j2, j3))
    g.add(Bond(e_sf, j2))
    g.add(Bond(j3, e_c))
    g.add(Bond(j3, e_r))

    eqs, vars = g.get_state_equations()
    assert vars == [p, q]
    assert eqs[1].equals(F - r * v - r * p / i - q / c)
    assert eqs[2].equals(v + p / i)
