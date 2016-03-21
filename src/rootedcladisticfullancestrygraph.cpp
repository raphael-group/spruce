/*
 *  rootedcladisticfullancestrygraph.cpp
 *
 *   Created on: 18-ocy-2015
 *       Author: M. El-Kebir
 */

#include "rootedcladisticfullancestrygraph.h"

namespace gm {
  
RootedCladisticFullAncestryGraph::RootedCladisticFullAncestryGraph(const StateTreeVector& S)
  : _S(S)
  , _G()
  , _root(lemon::INVALID)
  , _charStateToNode(_S.size())
  , _nodeToCharState(_G)
{
  int n = _S.size();
  for (int c = 0; c < n; ++c)
  {
    _charStateToNode[c] = NodeVector(_S[c].k(), lemon::INVALID);
  }
}
  
void RootedCladisticFullAncestryGraph::writeDOT(std::ostream& out) const
{
  out << "digraph G {" << std::endl;

  for (NodeIt v_ci(_G); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPair& ci = _nodeToCharState[v_ci];
    if (ci == std::make_pair(0, 0))
      out << "\t" << _G.id(v_ci) << " [label=\"(*," << _S[0].label(0) << ")\\n";
    else
      out << "\t" << _G.id(v_ci) << " [label=\"(" << ci.first << "," << _S[ci.first].label(ci.second) << ")\\n";

    const IntSet& D_ci = _S[ci.first].D(ci.second);
    out << "{";
    for (IntSetIt it = D_ci.begin(); it != D_ci.end(); ++it)
    {
      out << " " << *it;
    }
    out << " }\"]" << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    out << "\t" << _G.id(_G.source(a)) << " -> " << _G.id(_G.target(a)) << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void RootedCladisticFullAncestryGraph::init()
{
  const int n = _S.size();
  
  // add root node
  _root = _G.addNode();
  _nodeToCharState[_root] = std::make_pair(0, 0);
  for (int c = 0; c < n; ++c)
  {
    _charStateToNode[c][0] = _root;
  }
  
  // add the other n(k-1) nodes
  for (int c = 0; c < n; ++c)
  {
    const int k = _S[c].k();
    for (int i = 1; i < k; ++i)
    {
      if (_S[c].isPresent(i))
      {
        Node v_ci = _G.addNode();
        _charStateToNode[c][i] = v_ci;
        _nodeToCharState[v_ci] = std::make_pair(c, i);
      }
    }
  }

  // let's add the edges according to the state trees
  for (int c = 0; c < n; ++c)
  {
    const int k = _S[c].k();
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
      
      for (int i = 1; i < _S[c].k(); ++i)
      {
        // there's only one root vertex
        if (i == 0 && c != 0) continue;
        if (_charStateToNode[c][i] == lemon::INVALID) continue;
        
        for (int j = 1; j < _S[d].k(); ++j)
        {
          if (j == 0 && d != 0) continue;
          if (_charStateToNode[d][j] == lemon::INVALID) continue;
          
          // respect the state tree, also for the root vertex
          if (c == 0 && i == 0 && _S[d].parent(j) != 0)
            continue;
          
          _G.addArc(_charStateToNode[c][i], _charStateToNode[d][j]);
        }
      }
    }
  }
}

  
} // namespace gm