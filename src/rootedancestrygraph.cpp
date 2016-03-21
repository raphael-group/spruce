/*
 *  rootedancestrygraph.cpp
 *
 *  Created on: 10-auh-2015
 *      Author: M. El-Kebir
 */

#include "rootedancestrygraph.h"

namespace gm {

RootedAncestryGraph::RootedAncestryGraph(const RealTensor& F)
  : _F(F)
  , _G()
  , _charStateToNode(_F.n(),
                     NodeVector(_F.k(), lemon::INVALID))
  , _nodeToCharState(_G)
  , _compatibiltyMatrix(_F.n(), IntSetPairListVector(_F.n()))
  , _theta(_G)
  , _theta2(_G)
{
  initCompatibilityMatrix();
  initG();
}
  
void RootedAncestryGraph::initCompatibilityMatrix()
{
  const int m = _F.m();
  const int n = _F.n();
  
  for (int c = 0; c < n; ++c)
  {
    for (int d = 0; d < n; ++d)
    {
      _compatibiltyMatrix[c][d].clear();
      
      if (c == d) continue;
      IntSet S;
      while (next(c, S))
      {
//        std::cout << _F.getColLabel(c) << " -> " << _F.getColLabel(d) << ": " << " {";
//        for (IntSetIt it = S.begin(); it != S.end(); ++it)
//        {
//          std::cout << " " << *it;
//        }
//        std::cout << " }" << std::endl;
        IntSet T;
        while (next(d, T))
        {
//          std::cout << j << " {";
//          for (IntSetIt it = T.begin(); it != T.end(); ++it)
//          {
//            std::cout << " " << *it;
//          }
//          std::cout << " }" << std::endl;
          
          bool ok = true;
          for (int p = 0; p < m; ++p)
          {
            double cumS = _F.getCumFreq(p, c, S);
            double cumT = _F.getCumFreq(p, d, T);
            
            ok &= !g_tol.less(cumS, cumT);
          }
          if (ok)
          {
            _compatibiltyMatrix[c][d].push_back(std::make_pair(S, T));
          }
        }
      }
    }
  }
}
  
void RootedAncestryGraph::initG()
{
  const int n = _F.n();
  const int m = _F.m();
  const int k = _F.k();
  
  // add nodes
  for (int c = 0; c < n; ++c)
  {
    for (int i = 0; i < k; ++i)
    {
      Node v_ci = _G.addNode();
      _charStateToNode[c][i] = v_ci;
      _nodeToCharState[v_ci] = std::make_pair(c, i);
    }
  }
  
  // add edges
  for (int c = 0; c < n; ++c)
  {
    for (int i = 0; i < k; ++i)
    {
      Node v_ci = _charStateToNode[c][i];
      for (int d = 0; d < n; ++d)
      {
        const IntSetPairList& L = _compatibiltyMatrix[c][d];
        for (int j = 0; j < k; ++j)
        {
          IntSetFamily LL2;
          IntSetPairList LL;
          for (IntSetPairListIt it = L.begin(); it != L.end(); ++it)
          {
            const IntSet& S = it->first;
            const IntSet& T = it->second;
            if (S.find(i) != S.end() && T.find(j) != T.end())
            {
              LL2.insert(S);
              LL.push_back(*it);
            }
          }
          
          if (!LL.empty())
          {
            Node v_dj = _charStateToNode[d][j];
            Arc a = _G.addArc(v_ci, v_dj);
            _theta[a] = LL;
            _theta2[a] = LL2;
          }
        }
      }
    }
  }
}
  
void RootedAncestryGraph::writeDOT(const SubDigraph& G,
                                   const IntSetFamilyArcMap& theta,
                                   const IntSetNodeMap& sigma,
                                   std::ostream& out) const
{
  const int m = _F.m();
  out << "digraph G {" << std::endl;
  
  for (SubNodeIt v_ci(G); v_ci != lemon::INVALID; ++v_ci)
  {
    int c = _nodeToCharState[v_ci].first;
    int i = _nodeToCharState[v_ci].second;
    
    out << "\t" << _G.id(v_ci)
        << " [label=\"(" << _F.getColLabel(c) << "," << i << ")\\n{";
    
    const IntSet& S = sigma[v_ci];
    for (IntSetIt it = S.begin(); it != S.end(); ++it)
    {
      out << " " << *it;
    }
    out << " }\\n[";
    
    for (int p = 0; p < m; ++p)
    {
      double f = 0;
      const IntSet& S = sigma[v_ci];
      for (IntSetIt it = S.begin(); it != S.end(); ++it)
      {
        f += _F(*it, p, c);
      }
      out << " " << f;
    }
    out << " ]";
    
    out << "\"]" << std::endl;
  }
  
  for (SubArcIt a(G); a != lemon::INVALID; ++a)
  {
    Node v_ci = _G.source(a);
    Node v_dj = _G.target(a);
    out << "\t" << _G.id(v_ci) << " -> " << _G.id(v_dj) << " [label=\"";
    const IntSetFamily& SS = theta[a];
    out << "{ ";
    for (IntSetFamilyIt it = SS.begin(); it != SS.end(); ++it)
    {
      const IntSet& S = *it;
      out << "{";
      for (IntSetIt it2 = S.begin(); it2 != S.end(); ++it2)
      {
        out << " " << *it2;
      }
      out << " } ";
    }
    out << "}\"";

    if (SS.empty())
    {
      out << ",color=red";
    }
    out << "]" << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void RootedAncestryGraph::writeGraph(std::ostream& out) const
{
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node v_ci = _G.source(a);
    Node v_dj = _G.target(a);
    
    int c = getCharStatePair(v_ci).first;
    int i = getCharStatePair(v_ci).second;
    int d = getCharStatePair(v_dj).first;
    int j = getCharStatePair(v_dj).second;
    
    out << "(" << _F.getColLabel(c) << "," << i << ")"
        << " -> "
        << "(" << _F.getColLabel(d) << "," << j << ")" << std::endl;

    const IntSetPairList& L = _theta[a];
    for (IntSetPairListIt it = L.begin(); it != L.end(); ++it)
    {
      const IntSet& S = it->first;
      const IntSet& T = it->second;
      
      out << "{";
      for (IntSetIt it2 = S.begin(); it2 != S.end(); ++it2)
      {
        out << " " << *it2;
      }
      out << " } ";
      out << "{";
      for (IntSetIt it2 = T.begin(); it2 != T.end(); ++it2)
      {
        out << " " << *it2;
      }
      out << " }" << std::endl;
    }
  }
}
  
void RootedAncestryGraph::writeDOT(std::ostream& out) const
{
  out << "digraph G {" << std::endl;
  
  for (NodeIt v_ci(_G); v_ci != lemon::INVALID; ++v_ci)
  {
    int c = _nodeToCharState[v_ci].first;
    int i = _nodeToCharState[v_ci].second;
    
    out << "\t" << _G.id(v_ci)
        << " [label=\"(" << _F.getColLabel(c) << "," << i << ")\"]" << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node v_ci = _G.source(a);
    Node v_dj = _G.target(a);
    out << "\t" << _G.id(v_ci) << " -> " << _G.id(v_dj) << " [label=\"";
    const IntSetPairList& L = _theta[a];
    out << "{ ";
    for (IntSetPairListIt it = L.begin(); it != L.end(); ++it)
    {
      const IntSet& S = it->first;
      const IntSet& T = it->second;
      
      out << "[ {";
      for (IntSetIt it2 = S.begin(); it2 != S.end(); ++it2)
      {
        out << " " << *it2;
      }
      out << " } ";
      out << "{";
      for (IntSetIt it2 = T.begin(); it2 != T.end(); ++it2)
      {
        out << " " << *it2;
      }
      out << " } ] ";

    }
    out << "}" << "\"]" << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void RootedAncestryGraph::writeDOT2(std::ostream& out) const
{
  out << "digraph G {" << std::endl;
  
  for (NodeIt v_ci(_G); v_ci != lemon::INVALID; ++v_ci)
  {
    int c = _nodeToCharState[v_ci].first;
    int i = _nodeToCharState[v_ci].second;
    
    out << "\t" << _G.id(v_ci)
        << " [label=\"(" << _F.getColLabel(c) << "," << i << ")\"]" << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node v_ci = _G.source(a);
    Node v_dj = _G.target(a);
    out << "\t" << _G.id(v_ci) << " -> " << _G.id(v_dj) << " [label=\"";
    const IntSetFamily& L = _theta2[a];
    out << "{";
    for (IntSetFamilyIt it = L.begin(); it != L.end(); ++it)
    {
      const IntSet& S = *it;
      
      out << " {";
      for (IntSetIt it2 = S.begin(); it2 != S.end(); ++it2)
      {
        out << " " << *it2;
      }
      out << " }";
    }
    out << " }" << "\"]" << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void RootedAncestryGraph::writeCompatibilityMatrix(std::ostream& out) const
{
  const int n = _F.n();
  for (int c = 0; c < n; ++c)
  {
    for (int d = 0; d < n; ++d)
    {
      out << _F.getColLabel(c) << " -> " << _F.getColLabel(d) << std::endl;
      const IntSetPairList& L = _compatibiltyMatrix[c][d];
      for (IntSetPairListIt it = L.begin(); it != L.end(); ++it)
      {
        out << "{";
        const IntSet& S = it->first;
        for (IntSetIt it2 = S.begin(); it2 != S.end(); ++it2)
        {
          out << " " << *it2;
        }
        out << " }";

        out << " {";
        const IntSet& T = it->second;
        for (IntSetIt it2 = T.begin(); it2 != T.end(); ++it2)
        {
          out << " " << *it2;
        }
        out << " }" << std::endl;
      }
    }
  }
}
  
bool RootedAncestryGraph::next(int c,
                               IntSet& S) const
{
  const int k = _F.k();
  
#ifdef DEBUG
  for (IntSetIt it = S.begin(); it != S.end(); ++it)
  {
    assert(0 <= *it && *it < k);
    assert(*it != 0);
  }
#endif
  
  if (S.size() == k - 1)
    return false;
  
  // convert S to int;
  int s = 0;
  for (IntSetIt it = S.begin(); it != S.end(); ++it)
  {
    int i = *it;
    
    assert(0 <= i && i < k);
    assert(i != 0);
    if (i < 0)
    {
      s |= 1 << i;
    }
    else
    {
      s |= 1 << (i-1);
    }
  }
  
  // increment s
  ++s;
  
  // convert s back to S
  S.clear();
  for (int i = 0; i < k; ++i)
  {
    if (i < 0)
    {
      if ((s & (1 << i)) != 0)
      {
        S.insert(i);
      }
    }
    else if (i > 0)
    {
      if ((s & (1 << (i-1))) != 0)
      {
        S.insert(i);
      }
    }
  }
  
  return true;
}
  
} // namespace gm
