/*
 *  rootedcladisticancestrygraph.cpp
 *
 *   Created on: 29-sep-2015
 *       Author: M. El-Kebir
 */

#include "rootedcladisticancestrygraph.h"

#include <lemon/connectivity.h>

namespace gm {
  
RootedCladisticAncestryGraph::RootedCladisticAncestryGraph(const RealTensor& F,
                                                           const StateTreeVector& S)
  : _F(F)
  , _S(S)
  , _G()
  , _root(lemon::INVALID)
  , _charStateToNode(_F.n(), NodeVector(_F.k(), lemon::INVALID))
  , _nodeToCharState(_G)
  , _nodeToCharStateList(_G)
  , _label(_G)
{
}
  
void RootedCladisticAncestryGraph::makeMonoclonal(Node v_ci)
{
  std::list<Arc> toRemove;
  for (OutArcIt a_00dj(_G, _root); a_00dj != lemon::INVALID; ++a_00dj)
  {
    Node v_dj = _G.target(a_00dj);
    if (v_dj != v_ci)
    {
      toRemove.push_back(a_00dj);
    }
  }
  
  for (InArcIt a_djci(_G, v_ci); a_djci != lemon::INVALID; ++a_djci)
  {
    Node v_dj = _G.source(a_djci);
    if (v_dj != _root)
    {
      toRemove.push_back(a_djci);
    }
  }
  
  for (Arc a_00dj : toRemove)
  {
    _G.erase(a_00dj);
  }
  // TODO: what happens if G becomes disconnected?
  // it shouldn't
}
  
Digraph::Node RootedCladisticAncestryGraph::collapse(const IntPairSet& X)
{
  if (X.empty())
    return lemon::INVALID;
  
  std::list<IntPair> sortedX;
  sortedX.insert(sortedX.begin(), X.begin(), X.end());
  sortedX.sort(Compare(_S));
  
  const IntPair& firstCharState = sortedX.front();
  const IntPair& lastCharState = sortedX.back();
  
  // determine vertex set to collapse
  NodeSet V_X;

  // and determine outgoing and incoming vertex set
  NodeSet delta_plus, delta_minus;
  for (const IntPair& ci : X)
  {
    Node v_ci = charStateToNode(ci.first, ci.second);
    assert(v_ci != lemon::INVALID);
    V_X.insert(v_ci);
//    
//    for (OutArcIt a_cidj(_G, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
//    {
//      Node v_dj = _G.target(a_cidj);
//      const IntPairSet& XX = nodeToCharState(v_dj);
//
//      IntPairSet XXX;
//      std::set_intersection(X.begin(), X.end(), XX.begin(), XX.end(), std::inserter(XXX, XXX.begin()));
//      if (XXX.empty())
//      {
//        // v_dj not in V_X
//        delta_plus.insert(v_dj);
//      }
//    }
//    
//    for (InArcIt a_djci(_G, v_ci); a_djci != lemon::INVALID; ++a_djci)
//    {
//      Node v_dj = _G.source(a_djci);
//      const IntPairSet& XX = nodeToCharState(v_dj);
//      
//      IntPairSet XXX;
//      std::set_intersection(X.begin(), X.end(), XX.begin(), XX.end(), std::inserter(XXX, XXX.begin()));
//      if (XXX.empty())
//      {
//        delta_minus.insert(v_dj);
//      }
//    }
  }
  
  for (InArcIt a_cidj(_G, charStateToNode(firstCharState.first, firstCharState.second)); a_cidj != lemon::INVALID; ++a_cidj)
  {
    Node v_ci = _G.source(a_cidj);
    
    if (V_X.count(v_ci))
      continue;
    
    delta_minus.insert(v_ci);
  }
  for (OutArcIt a_cidj(_G, charStateToNode(lastCharState.first, lastCharState.second)); a_cidj != lemon::INVALID; ++a_cidj)
  {
    Node v_dj = _G.target(a_cidj);
    
    if (V_X.count(v_dj))
      continue;
    
    delta_plus.insert(v_dj);
  }
  
  // erase all vertices in V_X
  for (Node v_ci : V_X)
  {
    _G.erase(v_ci);
  }
  
  // insert new vertex
  Node v_ci = _G.addNode();
  _nodeToCharState[v_ci] = X;
  for (const IntPair& ci : X)
  {
    _charStateToNode[ci.first][ci.second] = v_ci;
  }
  
  // insert arcs
  for (Node v_dj : delta_minus)
  {
    // check if ok
    const IntPairSet& X_dj = nodeToCharState(v_dj);
    bool ok = true;
    for (const IntPair& dj : X_dj)
    {
      for (const IntPair& ci : X)
      {
        int pi_i = _S[ci.first].parent(ci.second);
        if (dj.first == ci.first && dj.second != pi_i && !X.count(IntPair(ci.first, pi_i)))
        {
          ok = false;
        }
      }
    }
    
    if (ok)
    {
      _G.addArc(v_dj, v_ci);
    }
  }
  for (Node v_dj : delta_plus)
  {
    // check if ok
    const IntPairSet& X_dj = nodeToCharState(v_dj);
    bool ok = true;
    for (const IntPair& dj : X_dj)
    {
      for (const IntPair& ci : X)
      {
        int pi_j = _S[dj.first].parent(dj.second);
        if (dj.first == ci.first && ci.second != pi_j && !X.count(IntPair(dj.first, pi_j)))
        {
          ok = false;
        }
      }
    }
    
    if (ok)
    {
      _G.addArc(v_ci, v_dj);
    }
  }
  
  _nodeToCharStateList[v_ci] = sortedX;
  
  assert(lemon::simpleGraph(_G));
  return v_ci;
}
  
void RootedCladisticAncestryGraph::writeDOT(std::ostream& out) const
{
  const int m = _F.m();
  out << "digraph G {" << std::endl;
  out.precision(3);

  for (NodeIt v_ci(_G); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPairSet& X_ci = _nodeToCharState[v_ci];
    out << "\t" << _G.id(v_ci) << " [label=\"" << _label[v_ci] << "\\n";

    for (const IntPair& ci : X_ci)
    {
      const IntSet& D_ci = _S[ci.first].D(ci.second);
      out << "{";
      for (IntSetIt it = D_ci.begin(); it != D_ci.end(); ++it)
      {
        out << " " << *it;
      }
      out << " } ";
    }
    out << "\\n";

    for (int p = 0; p < m; ++p)
    {
      bool first = true;
      for (const IntPair& ci : X_ci)
      {
        if (first)
          first = false;
        else
          out << " | ";
        out << _F.getCumFreq(p, ci.first, _S[ci.first].D(ci.second));
      }
      out << "\\n";
    }
    out << "\"]" << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    out << "\t" << _G.id(_G.source(a)) << " -> " << _G.id(_G.target(a)) << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void RootedCladisticAncestryGraph::init()
{
  const int k = _F.k();
  const int m = _F.m();
  const int n = _F.n();
  
  // add root node
  _root = _G.addNode();
  for (int c = 0; c < n; ++c)
  {
    _nodeToCharState[_root].insert(IntPair(c, 0));
    _charStateToNode[c][0] = _root;
  }
  
  // add the other n(k-1) nodes
  for (int c = 0; c < n; ++c)
  {
    for (int i = 1; i < k; ++i)
    {
      // we should not be removing states present in the state tree!
      if (!_S[c].isPresent(i))
        continue;

      Node v_ci = _G.addNode();
      _charStateToNode[c][i] = v_ci;
      _nodeToCharState[v_ci].insert(std::make_pair(c, i));
      _nodeToCharStateList[v_ci].push_back(std::make_pair(c, i));
    }
  }

  // let's add the edges according to the state trees
  for (int c = 0; c < n; ++c)
  {
    for (int i = 0; i < k; ++i)
    {
      Node v_ci = _charStateToNode[c][i];
      if (v_ci == lemon::INVALID) continue;
      for (int j = 1; j < k; ++j)
      {
        if (i == j) continue;
        
        Node v_cj = _charStateToNode[c][j];
        if (v_cj == lemon::INVALID) continue;
        if (_S[c].isParent(i, j))
        {
          _G.addArc(v_ci, v_cj);
        }
      }
    }
  }
  
  // now let's add edges according to frequency tensor (for distinct characters)
  for (int c = 0; c < n; ++c)
  {
    for (int d = 0; d < n; ++d)
    {
      if (c == d) continue;
      
      for (int i = 1; i < k; ++i)
      {
        // there's only one root vertex
        if (i == 0 && c != 0) continue;
        if (_charStateToNode[c][i] == lemon::INVALID) continue;
        
        for (int j = 1; j < k; ++j)
        {
          if (j == 0 && d != 0) continue;
          if (_charStateToNode[d][j] == lemon::INVALID) continue;
          
          // respect the state tree, also for the root vertex
          if (c == 0 && i == 0 && _S[d].parent(j) != 0)
            continue;
          
          bool ok = true;
          for (int p = 0; p < m; ++p)
          {
            double F_p_ci = _F.getCumFreq(p, c, _S[c].D(i));
            double F_p_dj = _F.getCumFreq(p, d, _S[d].D(j));

            ok &= !g_tol.less(F_p_ci, F_p_dj);
            if (!ok) break;
          }
          
          if (ok)
          {
            _G.addArc(_charStateToNode[c][i], _charStateToNode[d][j]);
          }
        }
      }
    }
  }
}
  
} // namespace gm
