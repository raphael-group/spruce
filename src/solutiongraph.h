/*
 *  rootedcladisticancestrygraph.h
 *
 *   Created on: 4-oct-2015
 *       Author: M. El-Kebir
 */

#ifndef SOLUTIONGRAPH_H
#define SOLUTIONGRAPH_H

#include <lemon/list_graph.h>
#include <lemon/core.h>
#include "utils.h"
#include "perfectphylotree.h"
#include "realtensor.h"
#include "realmatrix.h"
#include "solution.h"

namespace gm {
  
class SolutionGraph
{
public:
  typedef std::vector<StateTree> StateTreeVector;
  typedef lemon::ListBpGraph MixingGraph;
  DIGRAPH_TYPEDEFS(Digraph);
  typedef MixingGraph::Node BpNode;
  typedef MixingGraph::NodeIt BpNodeIt;
  typedef MixingGraph::BlueNode BpBlueNode;
  typedef MixingGraph::BlueNodeIt BpBlueNodeIt;
  typedef MixingGraph::RedNode BpRedNode;
  typedef MixingGraph::RedNodeIt BpRedNodeIt;
  typedef MixingGraph::Edge BpEdge;
  typedef MixingGraph::EdgeIt BpEdgeIt;
  typedef MixingGraph::RedNodeMap<int> IntBpRedNodeMap;
  typedef MixingGraph::BlueNodeMap<int> IntBpBlueNodeMap;
  typedef MixingGraph::NodeMap<StlBoolVector> BoolVectorBpNodeMap;
  
  SolutionGraph(const RealTensor& F,
                const StateTreeVector& S,
                const Solution& sol,
                double threshold,
                bool showStateVectors,
                bool showCumFreqs,
                bool showEdgeLabels);
  
  void writeUsageMatrix(std::ostream& out) const;
  
  void writeDOT(std::ostream& out) const;
  
  void writeASCII(std::ostream& out) const;
  
  const MixingGraph& getMixingGraph() const
  {
    return _G;
  }
    
private:
  typedef std::vector<Node> NodeVector;
  typedef Digraph::NodeMap<StlBoolVector> BoolVectorNodeMap;
  typedef Digraph::NodeMap<BpBlueNode> BpBlueNodeNodeMap;
  
  typedef std::vector<BpRedNode> BpRedNodeVector;
  typedef std::vector<BpBlueNode> BpBlueNodeVector;
  typedef MixingGraph::RedNodeMap<StlBoolVector> BoolVectorRedNodeMap;
  typedef MixingGraph::BlueNodeMap<Node> NodeBpBlueNodeMap;
  typedef MixingGraph::EdgeMap<double> MixingEdgeMap;
  
  const Solution& _sol;
  const RealTensor& _F;
  const StateTreeVector& _S;
  PerfectPhyloTree _T;

  BoolNodeMap _leaf;
  
  // red node: samples
  // blue node: deconvoluted samples
  lemon::ListBpGraph _G;
  IntBpRedNodeMap _bpRedNodeToSample;
  BpRedNodeVector _sampleToBpRedNode;
  BpBlueNodeNodeMap _charStateToBpBlueNode;
  NodeBpBlueNodeMap _bpBlueNodeToCharState;
  MixingEdgeMap _mixingFraction;
  double _threshold;
  bool _showStateVectors;
  bool _showCumFreqs;
  bool _showEdgeLabels;
  
  void constructMixingGraph();
  std::string toEdgeLabel(const IntPair& ci) const;
  std::string toNodeLabel(const IntPair& ci) const;
};

  
} // namespace gm

#endif // SOLUTIONGRAPH_H