/*
 *  rootedcladisticfullancestrygraph.h
 *
 *   Created on: 18-ocy-2015
 *       Author: M. El-Kebir
 */

#ifndef ROOTEDCLADISTICFULLANCESTRYGRAPH_H
#define ROOTEDCLADISTICFULLANCESTRYGRAPH_H

#include "utils.h"
#include "statetree.h"
#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include <set>

namespace gm {
  
class RootedCladisticFullAncestryGraph
{
public:
  DIGRAPH_TYPEDEFS(Digraph);
  typedef std::vector<Node> NodeVector;
  typedef std::vector<NodeVector> NodeMatrix;
  
  typedef std::vector<StateTree> StateTreeVector;

  typedef Digraph::NodeMap<IntPair> IntPairNodeMap;
  typedef std::set<int> IntSet;
  typedef IntSet::const_iterator IntSetIt;
  typedef std::set<IntSet> IntSetFamily;
  typedef IntSetFamily::const_iterator IntSetFamilyIt;
  typedef IntSetFamily::iterator IntSetFamilyNonConstIt;
  typedef Digraph::ArcMap<IntSetFamily> IntSetFamilyArcMap;
  typedef std::pair<int, IntSetFamily> CharIntSetFamilyPair;
  typedef std::vector<IntSetFamily> IntSetFamilyVector;
  typedef Digraph::NodeMap<IntSetFamilyVector> IntSetFamilyVectorMap;
  typedef Digraph::NodeMap<CharIntSetFamilyPair> CharIntSetFamilyPairNodeMap;
  typedef std::set<Arc> ArcSet;
  typedef ArcSet::const_iterator ArcSetIt;
  typedef std::set<Node> NodeSet;
  typedef NodeSet::const_iterator NodeSetIt;
  typedef Digraph::NodeMap<IntSet> IntSetNodeMap;
  typedef lemon::SubDigraph<const Digraph> SubDigraph;
  typedef SubDigraph::ArcIt SubArcIt;
  typedef SubDigraph::NodeIt SubNodeIt;
  typedef SubDigraph::OutArcIt SubOutArcIt;
  typedef SubDigraph::InArcIt SubInArcIt;
  
  typedef std::pair<IntSet, IntSet> IntSetPair;
  typedef std::list<IntSetPair> IntSetPairList;
  typedef IntSetPairList::const_iterator IntSetPairListIt;
  typedef std::vector<IntSetPairList> IntSetPairListVector;
  typedef IntSetPairListVector::const_iterator IntSetPairListVectorIt;
  typedef std::vector<IntSetPairListVector> IntSetPairListMatrix;
  typedef IntSetPairListMatrix::const_iterator IntSetPairListMatrixIt;
  typedef Digraph::ArcMap<IntSetPairList> IntSetPairListArcMap;
  
  RootedCladisticFullAncestryGraph(const StateTreeVector& S);
  
  virtual void writeDOT(std::ostream& out) const;
  
  const StateTree& S(int c) const
  {
    assert(0 <= c && c < _S.size());
    return _S[c];
  }
  
  const StateTreeVector& S() const
  {
    return _S;
  }
  
  const Digraph& G() const
  {
    return _G;
  }
  
  Node root() const
  {
    return _root;
  }
  
  Node charStateToNode(int c, int i) const
  {
    assert(0 <= c && c < _S.size());
    assert(0 <= i && i < _S[c].k());
    
    return _charStateToNode[c][i];
  }
  
  const IntPair& nodeToCharState(Node v_ci) const
  {
    return _nodeToCharState[v_ci];
  }
  
  virtual void init();
  
protected:
  const StateTreeVector& _S;
  
  Digraph _G;
  Node _root;
  NodeMatrix _charStateToNode;
  IntPairNodeMap _nodeToCharState;
};
  
} // namespace gm

#endif // ROOTEDCLADISTICFULLANCESTRYGRAPH_H