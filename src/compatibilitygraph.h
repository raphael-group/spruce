/*
 * compatibilitygraph.h
 *
 *  Created on: 29-jan-2016
 *      Author: M. El-Kebir
 */

#ifndef COMPATABILITYGRAPH_H
#define COMPATABILITYGRAPH_H

#include "utils.h"
#include "charactermatrix.h"
#include "bronkerbosch.h"
#include "solutionset.h"
#include <map>
#include <set>

namespace gm {
  
class CompatibilityGraph
{
public:
  typedef std::vector<StateTree> StateTreeVector;
  
  typedef std::map<int, int> IntIntMap;
  
  CompatibilityGraph(const CharacterMatrix& M);
  
  void write(std::ostream& out) const;
  
  void writeDOT(std::ostream& out) const;
  
  unsigned long combinations() const
  {
    return _combinations;
  }
  
  void init(unsigned long combination,
            RealTensor& F,
            StateTreeVector& S,
            RealTensor& F_lb,
            RealTensor& F_ub,
            StlIntVector& mapNewCharToOldChar,
            StlIntVector& mapOldCharToNewChar) const;
  
  void init(const IntPairSet& filter, int size);
  
  void init(std::ifstream& inFile);
  
protected:
  GRAPH_TYPEDEFS(Graph);
  typedef Graph::NodeMap<IntPair> IntPairNodeMap;
  
  typedef std::vector<Node> NodeVector;
  typedef std::set<Node> NodeSet;
  typedef NodeVector::const_iterator NodeVectorIt;
  typedef std::vector<NodeVector> NodeMatrix;
  typedef NodeMatrix::const_iterator NodeMatrixIt;
  
  void initVertices();
  void initEdges();
  
  bool isCompatible(const IntPair& cs, const IntPair& dt) const;
  void applyFilter(const IntPairSet& filter, const NodeMatrix& cliques);
  
private:
  const CharacterMatrix& _M;
  Graph _G;
  IntPairNodeMap _nodeToCharStateTree;
  NodeMatrix _charStateTreeToNode;
  
  NodeMatrix _cliques;
  unsigned long _combinations;
  StlIntVector _mapping;
};
  
} // namespace gm

#endif // COMPATABILITYGRAPH_H
