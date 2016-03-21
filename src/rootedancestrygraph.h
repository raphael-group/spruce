/*
 *  rootedancestrygraph.h
 *
 *   Created on: 10-aug-2015
 *       Author: M. El-Kebir
 */

#ifndef ROOTEDANCESTRYGRAPH_H
#define ROOTEDANCESTRYGRAPH_H

#include "utils.h"
#include "realtensor.h"
#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include <set>

namespace gm {

class RootedAncestryGraph
{
public:
  typedef lemon::ListDigraph Digraph;
  DIGRAPH_TYPEDEFS(Digraph);
  typedef std::vector<Node> NodeVector;
  typedef std::vector<NodeVector> NodeMatrix;
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
  
  RootedAncestryGraph(const RealTensor& F);
  
  void writeGraph(std::ostream& out) const;
  
  void writeDOT(std::ostream& out) const;
  
  void writeDOT2(std::ostream& out) const;
  
  void writeDOT(const SubDigraph& G,
                const IntSetFamilyArcMap& theta,
                const IntSetNodeMap& sigma,
                std::ostream& out) const;
  
  void writeCompatibilityMatrix(std::ostream& out) const;
  
  
  const Digraph& getG() const { return _G; }
  
  Node getNode(int c, int i) const
  {
    return _charStateToNode[c][i];
  }
  
  const IntPair& getCharStatePair(Node v_ci) const
  {
    return _nodeToCharState[v_ci];
  }
    
  const RealTensor& getF() const
  {
    return _F;
  }
  
private:
  const RealTensor& _F;
  Digraph _G;
  
  NodeMatrix _charStateToNode;
  IntPairNodeMap _nodeToCharState;
  IntSetPairListMatrix _compatibiltyMatrix;
  
  IntSetPairListArcMap _theta;
  IntSetFamilyArcMap _theta2;
  
  void initCompatibilityMatrix();
  void initG();
  bool next(int c, IntSet& S) const;
};
  
} // namespace gm

#endif // ROOTEDANCESTRYGRAPH_H
