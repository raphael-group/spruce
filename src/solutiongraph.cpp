/*
 *  solutiongraph.cpp
 *
 *   Created on: 4-oct-2015
 *       Author: M. El-Kebir
 */

#include "solutiongraph.h"
#include <lemon/bfs.h>
#include <iomanip>
#include <stdio.h>

namespace gm {

SolutionGraph::SolutionGraph(const RealTensor& F,
                             const StateTreeVector& S,
                             const Solution& sol,
                             double threshold,
                             bool showStateVectors,
                             bool showCumFreqs,
                             bool showEdgeLabels)
  : _sol(sol)
  , _F(F)
  , _S(S)
  , _T(sol.A(), S)
  , _leaf(_T.T())
  , _G()
  , _bpRedNodeToSample(_G)
  , _sampleToBpRedNode(F.m(), lemon::INVALID)
  , _charStateToBpBlueNode(_T.T())
  , _bpBlueNodeToCharState(_G)
  , _mixingFraction(_G)
  , _threshold(threshold)
  , _showStateVectors(showStateVectors)
  , _showCumFreqs(showCumFreqs)
  , _showEdgeLabels(showEdgeLabels)
{
  _T.setLabels(F);
  const Digraph& T = _T.T();
  for (NodeIt v(T); v != lemon::INVALID; ++v)
  {
    _leaf[v] = OutArcIt(T, v) == lemon::INVALID;
  }
  constructMixingGraph();
}

void SolutionGraph::constructMixingGraph()
{
  const Digraph& T = _T.T();
  const int m = _F.m();
  const int n = _F.n();

  lemon::mapFill(T, _charStateToBpBlueNode, lemon::INVALID);
  
  // let's first add the sample nodes
  for (int p = 0; p < m; ++p)
  {
    BpRedNode v = _G.addRedNode();
    _sampleToBpRedNode[p] = v;
    _bpRedNodeToSample[v] = p;
  }
  
  // now let's add the mixing nodes
  for (NodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPair& ci = _T.nodeToCharState(v_ci);

    for (int p = 0; p < m; ++p)
    {
      int row_idx = ci.first == 0 && ci.second == 0 ? 0 : n * (ci.second - 1) + ci.first + 1;
      double u_pci = _sol.U()(p, row_idx);
      if (u_pci != 0 && u_pci >= _threshold)
      {
        BpBlueNode v = _charStateToBpBlueNode[v_ci];
        if (v == lemon::INVALID)
        {
          v = _G.addBlueNode();
          _bpBlueNodeToCharState[v] = v_ci;
          _charStateToBpBlueNode[v_ci] = v;
        }
        
        // add edge
        BpEdge e = _G.addEdge(_sampleToBpRedNode[p], v);
        _mixingFraction[e] = u_pci;
      }
    }
  }
}
  
void SolutionGraph::writeCloneTree(std::ostream& out,
                                   bool samples) const
{
  const Digraph& T = _T.T();
  const int m = _F.m();
  const int n = _F.n();
  
  // skip adding samples to the root
  for (ArcIt a_cidj(T); a_cidj != lemon::INVALID; ++a_cidj)
  {
    // skip the root
    Node v_ci = T.source(a_cidj);
    if (v_ci == _T.root()) continue;
    
    Node v_dj = T.target(a_cidj);
    
    out << _T.label(v_ci) << " " << _T.label(v_dj) << std::endl;
    if (samples)
    {
      const IntPair& dj = _T.nodeToCharState(v_dj);
      for (int p = 0; p < m; ++p)
      {
        assert(!(dj.first == 0 && dj.second == 0));
        int row_idx = n * (dj.second - 1) + dj.first + 1;
        double u_pdj = _sol.U()(p, row_idx);
        if (g_tol.nonZero(u_pdj))
        {
          std::string sample_dj = _T.label(v_dj) + "_" + _F.getRowLabel(p);
          // add edge
          out << _T.label(v_dj) << " " << sample_dj << std::endl;
        }
      }
    }
  }
}
  
void SolutionGraph::writeLeafLabeling(std::ostream& out) const
{
  const Digraph& T = _T.T();
  const int m = _F.m();
  const int n = _F.n();
  
  // skip adding samples to the root
  for (ArcIt a_cidj(T); a_cidj != lemon::INVALID; ++a_cidj)
  {
    Node v_dj = T.target(a_cidj);
    
    const IntPair& dj = _T.nodeToCharState(v_dj);
    for (int p = 0; p < m; ++p)
    {
      assert(!(dj.first == 0 && dj.second == 0));
      int row_idx = n * (dj.second - 1) + dj.first + 1;
      double u_pdj = _sol.U()(p, row_idx);
      if (g_tol.nonZero(u_pdj))
      {
        std::string sample_dj = _T.label(v_dj) + "_" + _F.getRowLabel(p);
        // add edge
        out << sample_dj << " " << _F.getRowLabel(p) << std::endl;
      }
    }
  }
}
  
void SolutionGraph::writeASCII(std::ostream& out) const
{
  const Digraph& T = _T.T();
  
  out << _T.numOfVertices() << " #nodes" << std::endl;
  for (NodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    out << T.id(v_ci) << " " << toNodeLabel(_T.nodeToCharState(v_ci)) << std::endl;
  }

  out << _T.numOfEdges() << " #edges" << std::endl;
  for (ArcIt a_cidj(T); a_cidj != lemon::INVALID; ++a_cidj)
  {
    Node v_ci = T.source(a_cidj);
    Node v_dj = T.target(a_cidj);
    
    out << T.id(v_ci) << " -> " << T.id(v_dj) << " " << toEdgeLabel(_T.nodeToCharState(v_dj)) << std::endl;
  }
}
  
void SolutionGraph::writeUsageMatrix(std::ostream& out) const
{
  int m = _F.m();
  int n = _F.n();
  int k = _F.k();
  
  for (int p = 0; p < m; ++p)
  {
    out << "\t" << _F.getRowLabel(p);
  }
  out << std::endl;
  
  const Digraph& T = _T.T();
  for (NodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    const IntPair& ci = _T.nodeToCharState(v_ci);
    
    out << _T.label(v_ci);
    for (int p = 0; p < m; ++p)
    {
      int row_idx = ci.first == 0 && ci.second == 0 ? 0 : n * (ci.second - 1) + ci.first + 1;
      double u_pci = _sol.U()(p, row_idx);
      out << "\t" << u_pci;
    }
    out << std::endl;
  }
}

void SolutionGraph::writeDOT(std::ostream& out) const
{
  static int fontsizeBox = 60;
  static int fontsize = 35;
  static int minPenwidth = 5;
  
  const Digraph& T = _T.T();
  
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << "graph G {" << std::endl;
  
  out << "\tfontsize=" << fontsize << std::endl;
  out << "\tlabel=\"|V| = " << lemon::countNodes(T) << ", |U| = "
      << lemon::countEdges(_G) << "\"" << std::endl;
  
  out << "\tsubgraph mixed {" << std::endl;
  for (BpRedNodeIt v(_G); v != lemon::INVALID; ++v)
  {
    out << "\t\ts" << _G.id(v)
        << " [colorscheme=paired10,penwidth=" << minPenwidth
        << ",fontsize=" << fontsizeBox
        << ",color=" << (_G.id(v) % 10) + 1 << ",shape=box,label=\"";
    out << _F.getRowLabel(_bpRedNodeToSample[v]) << "\"]" << std::endl;
  }
  out << "\t}" << std::endl;
  
  lemon::Bfs<Digraph> bfs(T);
  bfs.run(_T.root());
  int maxLevel = lemon::mapMaxValue(T, bfs.distMap());
  
  out << "\tsubgraph unmixed {" << std::endl;
  for (int l = 0; l < maxLevel; ++l)
  {
    out << "\t\tsubgraph " << l << " {" << std::endl;
    out << "\t\t\trank=same" << std::endl;
    for (NodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
    {
      if (bfs.dist(v_ci) == l && !_leaf[v_ci])
      {
        const IntPair& ci = _T.nodeToCharState(v_ci);
        out << "\t\t\t" << T.id(v_ci)
            << " [penwidth=" << minPenwidth
            << ",fontsize=" << fontsize
            << ",label=\"";
        out << toNodeLabel(ci);
        
        if (_showCumFreqs)
        {
          for (int p = 0; p < _F.m(); ++p)
          {
            out << "\\n" << _F(ci.second, p, ci.first);
          }
        }
        
        out << "\"";
        out << "]" << std::endl;
      }
    }
    out << "\t\t}" << std::endl;
  }
  
  // leaves
  out << "\t\tsubgraph " << "leaves" << " {" << std::endl;
  out << "\t\t\trank=same" << std::endl;
  for (NodeIt v_ci(T); v_ci != lemon::INVALID; ++v_ci)
  {
    if (_leaf[v_ci])
    {
      const IntPair& ci = _T.nodeToCharState(v_ci);
      out << "\t\t\t" << T.id(v_ci)
          << " [penwidth=" << minPenwidth
          << ",fontsize=" << fontsize
          << ",label=\"";
      out << toNodeLabel(ci);
      
      if (_showCumFreqs)
      {
        for (int p = 0; p < _F.m(); ++p)
        {
          out << "\\n" << _F(ci.second, p, ci.first);
        }
      }
      
      out << "\"";
      out << "]" << std::endl;
    }
  }
  for (BpBlueNodeIt v(_G); v != lemon::INVALID; ++v)
  {
    Node vv = _bpBlueNodeToCharState[v];
    const IntPair& ci = _T.nodeToCharState(vv);
    if (vv != lemon::INVALID && !_leaf[vv])
    {
      out << "\t\t\tdup" << _G.id(v)
          << " [penwidth=" << minPenwidth
          << ",fontsize=" << fontsize;
      out << ",label=\"";
      out << toNodeLabel(ci);
      out << "\"";
      out << "]" << std::endl;
    }
  }
  out << "\t\t}" << std::endl;
  out << "\t}" << std::endl;
  
  for (ArcIt a(T); a != lemon::INVALID; ++a)
  {
    Node u = T.source(a);
    Node v = T.target(a);
    const IntPair& ci = _T.nodeToCharState(v);
    out << "\t" << T.id(u) << " -- " << T.id(v)
        << " [fontsize=" << fontsize;
    if (_showEdgeLabels)
    {
      out << ",label=\"" << toEdgeLabel(ci)
          << "\"";
    }
    out << ",colorscheme=accent8,penwidth=3]" << std::endl;
  }
  
  for (BpEdgeIt e(_G); e != lemon::INVALID; ++e)
  {
    assert(_G.valid(e));
    if (_mixingFraction[e] < _threshold)
      continue;
    
    BpBlueNode v = _G.blueNode(e);
    BpRedNode u = _G.redNode(e);
    Node vv = _bpBlueNodeToCharState[v];
    if (vv != lemon::INVALID && _leaf[vv])
    {
      out << "\t" << T.id(vv) << " -- "
          << "s" << _G.id(u);
    }
    else
    {
      out << "\tdup" << _G.id(v) << " -- "
          << "s" << _G.id(u);
    }
    out << " [splines=none,colorscheme=paired10,color=" << (_G.id(u) % 10) + 1 << ",minlen=4,fontsize="
        << fontsize << ",label=" << std::setprecision(2)
        << _mixingFraction[e] << ",penwidth="
        << minPenwidth + 25 * _mixingFraction[e] << "]" << std::endl;
  }
  
  for (NodeIt v(T); v != lemon::INVALID; ++v)
  {
    BpBlueNode vv = _charStateToBpBlueNode[v];
    if (!_leaf[v] && vv != lemon::INVALID)
    {
      out << "\t" << T.id(v) << " -- " << "dup" << _G.id(vv)
          << " [penwidth=" << minPenwidth << ",style=dashed]" << std::endl;
    }
  }
  
  out << "}" << std::endl;
}
  
std::string SolutionGraph::toNodeLabel(const IntPair& ci) const
{
  if (_showStateVectors)
  {
    const int n = _F.n();
    std::string label;
    for (int d = 0; d < n; ++d)
    {
      if (d != 0)
        label += ",";
      label += _S[ci.first].label(_sol.A()(ci, d));
    }
    return label;
  }
  else
  {
    if (ci.first == 0 && ci.second == 0)
    {
      return "(*,0)";
    }
    else
    {
      return _T.label(_T.charStateToNode(ci.first, ci.second));
    }
  }
}
  
std::string SolutionGraph::toEdgeLabel(const IntPair& ci) const
{
  char buf[1024];
  
  if (ci.second == 0)
  {
    return "(*,0)";
  }
  else
  {
    snprintf(buf, 1024, "(%s,%s)",
             _F.getColLabel(ci.first).c_str(),
             _S[ci.first].label(ci.second).c_str());
    return buf;
  }
}

} // namespace gm
