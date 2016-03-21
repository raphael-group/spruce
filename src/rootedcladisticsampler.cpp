/*
 * rootedcladisticsampler.cpp
 *
 *  Created on: 19-oct-2015
 *      Author: M. El-Kebir
 */

#include "rootedcladisticsampler.h"
#include <lemon/bfs.h>
#include <map>
#include <stdlib.h>

namespace gm {
  
RootedCladisticSampler::RootedCladisticSampler(const RootedCladisticFullAncestryGraph& G,
                                               int seed)
  : _G(G)
  , _result()
  , _generator(seed)
{
}
  
void RootedCladisticSampler::init(SubDigraph& subG,
                                  SubDigraph& T,
                                  ArcList& F)
{
  const Digraph& G = _G.G();
  Node root = _G.root();

  T.enable(root);
  F.clear();
  for (OutArcIt a(G, root); a != lemon::INVALID; ++a)
  {
    F.push_back(a);
  }  
}
  
void RootedCladisticSampler::sample()
{
  const Digraph& G = _G.G();
  
  BoolNodeMap filterNodesT(G, false);
  BoolArcMap filterArcsT(G, false);
  SubDigraph T(G, filterNodesT, filterArcsT);
  
  BoolNodeMap filterNodesG(G, true);
  BoolArcMap filterArcsG(G, true);
  SubDigraph subG(G, filterNodesG, filterArcsG);
  
  ArcList F;
  IntSetNodeMap D(G, IntSet());
  
  init(subG, T, F);
  grow(subG, T, F);
}
  
bool RootedCladisticSampler::isAncestor(const SubDigraph& T,
                                        Node v_dj,
                                        Node v_ci) const
{
  assert(T.status(v_ci));
  
  Node root = _G.root();
  
  while (v_ci != root && v_ci != v_dj)
  {
    Arc a(SubInArcIt(T, v_ci));
    v_ci = T.source(a);
  }
  
  return v_ci == v_dj;
}
  
void RootedCladisticSampler::addArc(SubDigraph& T,
                                    Arc a_cidj) const
{
  const Node v_dj = T.target(a_cidj);
  
  // add a_cidj to T
  T.enable(v_dj);
  T.enable(a_cidj);
  
  assert(isArborescence(T));
  assert(isValid(T));
}
  
void RootedCladisticSampler::removeArc(SubDigraph& T,
                                       Arc a_cidj) const
{
  assert(T.status(a_cidj));
  
  const Node v_dj = T.target(a_cidj);
  
  // remove a_cidj from T
  T.disable(a_cidj);
  T.disable(v_dj);
  
  assert(isArborescence(T));
  assert(isValid(T));
}
  
bool RootedCladisticSampler::isArborescence(const SubDigraph& T) const
{
  assert(T.status(_G.root()));
  lemon::Bfs<SubDigraph> bfs(T);
  bfs.run(_G.root());
  
  for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    if (!bfs.reached(v_ci))
    {
      return false;
    }
  }
  
  return true;
}
  
bool RootedCladisticSampler::isValid(const SubDigraph& T) const
{
  Node root = _G.root();
  for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    if (v_ci == root)
      continue;
    
    const IntPair& ci = _G.nodeToCharState(v_ci);
    int pi_i = _G.S(ci.first).parent(ci.second);
    Node v_c_pi_i = _G.charStateToNode(ci.first, pi_i);
    if (!isFirstAncestor(T, ci.first, v_c_pi_i, v_ci))
    {
      return false;
    }
  }
  
  return true;
}

bool RootedCladisticSampler::isValid(const SubDigraph& T,
                                     Arc a_ciel) const
{
  T.enable(a_ciel);
  T.enable(_G.G().target(a_ciel));

  bool res = isValid(T);
  
  T.disable(a_ciel);
  T.disable(_G.G().target(a_ciel));
  
  return res;
}

void RootedCladisticSampler::writeDOT(std::ostream& out,
                                      const SubDigraph& T) const
{
  const Digraph& G = _G.G();
  
  out << "digraph T {" << std::endl;
  out.precision(3);
  
  for (SubNodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPair& ci = _G.nodeToCharState(v_ci);
    if (ci == std::make_pair(0, 0))
      out << "\t" << G.id(v_ci) << " [label=\"(*," << _G.S(ci.first).label(ci.second) << ")\\n";
    else
      out << "\t" << G.id(v_ci) << " [label=\"(" << ci.first << "," << _G.S(ci.first).label(ci.second) << ")\\n";
    
    const IntSet& D_ci = _G.S(ci.first).D(ci.second);
    out << "{";
    for (IntSetIt it = D_ci.begin(); it != D_ci.end(); ++it)
    {
      out << " " << *it;
    }
    out << " }\"]" << std::endl;
  }
  
  for (SubArcIt a(T); a != lemon::INVALID; ++a)
  {
    out << "\t" << G.id(G.source(a)) << " -> " << G.id(G.target(a)) << std::endl;
  }
  out << "}" << std::endl;
}
  
void RootedCladisticSampler::writeEdgeList(std::ostream& out) const
{
  const Digraph& G = _G.G();
  for (ArcListIt it = _result.begin(); it != _result.end(); ++it)
  {
    Node v_ci = G.source(*it);
    Node v_dj = G.target(*it);
    
    const IntPair& ci = _G.nodeToCharState(v_ci);
    const IntPair& dj = _G.nodeToCharState(v_dj);
    
//    if (ci == std::make_pair(0, 0))
//      out << "(*," << _G.S(ci.first).label(ci.second) << ") -> ";
//    else
    out << "(" << ci.first << "," << _G.S(ci.first).label(ci.second) << ") -> ";
    

    out << "(" << dj.first << "," << _G.S(dj.first).label(dj.second) << ")" << std::endl;
  }
}
  
void RootedCladisticSampler::grow(SubDigraph& G,
                                  SubDigraph& T,
                                  ArcList& F)
{
  typedef std::map<Node, ArcList> NodeArcMap;
  typedef NodeArcMap::const_iterator NodeArcMapIt;
  
  while (!F.empty())
  {
    // collect all source nodes of F
    NodeArcMap bySources;
    for (ArcListIt it = F.begin(); it != F.end(); ++it)
    {
      Node v_ci = T.source(*it);
      bySources[v_ci].push_back(*it);
    }
    
    // pick source node and arc at random
    int idx = rand(0, bySources.size() - 1);
    NodeArcMapIt it = bySources.begin();
    for (; idx > 0; --idx, ++it);
    Node v_ci = it->first;
    
    // pick arc
    idx = rand(0, it->second.size() - 1);
    ArcListIt it2 = it->second.begin();
    for (; idx > 0; --idx, ++it2);
    
    Arc a_cidj = *it2;
    
    // remove a_cidj from F
    for (ArcListNonConstIt it3 = F.begin(); it3 != F.end(); ++it3)
    {
      if (*it3 == a_cidj)
      {
        F.erase(it3);
        break;
      }
    }
    
    Node v_dj = G.target(a_cidj);
    assert(T.status(v_ci));
    assert(!T.status(v_dj));
    assert(!T.status(a_cidj));
    
    // add a_cidj to T
    addArc(T, a_cidj);
    
    // remove each arc wv where w in T from F
    for (ArcListNonConstIt it = F.begin(); it != F.end();)
    {
      if (G.target(*it) == v_dj)
      {
        assert(T.status(G.source(*it)));
        it = F.erase(it);
      }
      else
      {
        assert(isValid(T, *it));
        ++it;
      }
    }
    
    // push each arc a_djel where v_el not in V(T) onto F
    for (SubOutArcIt a_djel(G, v_dj); a_djel != lemon::INVALID; ++a_djel)
    {
      Node v_el = G.target(a_djel);
      const IntPair& el = _G.nodeToCharState(v_el);
      int pi_l = _G.S(el.first).parent(el.second);
      Node v_e_pi_l = _G.charStateToNode(el.first, pi_l);
      if (!T.status(v_el) && isFirstAncestor(T, el.first, v_e_pi_l, v_dj))
      {
        F.push_back(a_djel);
      }
    }
  }

  if (F.empty())
  {
    generateResult(T, _G.root());
  }
}

bool RootedCladisticSampler::isFirstAncestor(const SubDigraph& T,
                                             int c,
                                             Node v_ci,
                                             Node v_dj) const
{
  assert(T.status(v_dj));

  if (!T.status(v_ci))
  {
    return false;
  }
  if (v_ci == v_dj)
  {
    return true;
  }
  
  Arc a(SubInArcIt(T, v_dj));
  while (a != lemon::INVALID)
  {
    v_dj = T.source(a);
    a = SubInArcIt(T, v_dj);
    
    if (_G.nodeToCharState(v_dj).first == c)
      break;
  }
  
  return v_ci == v_dj;
}

  
void RootedCladisticSampler::generateResult(const SubDigraph& T,
                                            Node v_ci)
{
  for (SubOutArcIt a_cidj(T, v_ci); a_cidj != lemon::INVALID; ++a_cidj)
  {
    _result.push_back(a_cidj);
    generateResult(T, T.target(a_cidj));
  }
}
  
} // namespace gm