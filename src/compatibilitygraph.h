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

namespace gm {
  
class CompatibilityGraph
{
public:
  typedef std::vector<StateTree> StateTreeVector;
  
  CompatibilityGraph(const CharacterMatrix& M,
                     bool shuffle);
  
  CompatibilityGraph(const CharacterMatrix& M,
                     bool shuffle,
                     std::istream& inFile);
  
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
  
protected:
  GRAPH_TYPEDEFS(Graph);
  typedef Graph::NodeMap<IntPair> IntPairNodeMap;
  
  typedef std::vector<Node> NodeVector;
  typedef NodeVector::const_iterator NodeVectorIt;
  typedef std::vector<NodeVector> NodeMatrix;
  
  void init(bool initEdges);
  void solve();
  
  bool isCompatible(const IntPair& cs, const IntPair& dt) const;
  
private:
  const CharacterMatrix& _M;
  Graph _G;
  IntPairNodeMap _nodeToCharStateTree;
  NodeMatrix _charStateTreeToNode;
  
  NodeMatrix _maximumCliques;
  unsigned long _combinations;
  StlIntVector _mapping;
};
  
} // namespace gm

#endif // COMPATABILITYGRAPH_H
