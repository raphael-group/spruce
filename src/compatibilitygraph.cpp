/*
 * compatibilitygraph.cpp
 *
 *  Created on: 29-jan-2016
 *      Author: M. El-Kebir
 */

#include "compatibilitygraph.h"
#include "realtensor.h"
#include "rootedcladisticnoisyancestrygraph.h"
#include "rootedcladisticnoisyenumeration.h"
#include "bronkerbosch.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace gm {
  
CompatibilityGraph::CompatibilityGraph(const CharacterMatrix& M)
  : _M(M)
  , _G()
  , _cliqueLimit(-1)
  , _nodeToCharStateTree(_G)
  , _charStateTreeToNode()
  , _combinations(0)
  , _mapping()
{
  initVertices();
}
  
void CompatibilityGraph::init(std::ifstream& inFile)
{
  g_lineNumber = 0;
  std::string line;
  
  int compatibilityCount = -1;
  gm::getline(inFile, line);
  std::stringstream ss(line.c_str());
  ss >> compatibilityCount;
  
  for (int idx = 0; idx < compatibilityCount; ++idx)
  {
    gm::getline(inFile, line);
    int c = -1, d = -1, s = -1, t = -1;
    sscanf(line.c_str(), "(%d,%d) (%d,%d)", &c, &s, &d, &t);
    
    Node v_cs = _charStateTreeToNode[c][s];
    Node v_dt = _charStateTreeToNode[d][t];
    
    _G.addEdge(v_cs, v_dt);
  }
  
  int maxCliqueSize = -1;
  gm::getline(inFile, line);
  ss.clear();
  ss.str(line);
  ss >> maxCliqueSize;
  
  int cliqueCount = -1;
  gm::getline(inFile, line);
  ss.clear();
  ss.str(line);
  ss >> cliqueCount;
  
  // read mapping
  StringVector s;
  
  // read cliques
  for (int idx = 0; idx < cliqueCount; ++idx)
  {
    s.clear();
    gm::getline(inFile, line);
    boost::split(s, line, boost::is_any_of(" "));
    
    assert(s.size() == maxCliqueSize);
    NodeVector C;
    for (const std::string& str : s)
    {
      int c = -1, s = -1;
      if (sscanf(str.c_str(), "(%d,%d)", &c, &s) != 2)
      {
        throw std::runtime_error(getLineNumber() + "Error: '" + str
                                 + "' is not a pair of a character index and a state tree index");
      }
      
      C.push_back(_charStateTreeToNode[c][s]);
    }
    _cliques.push_back(C);
  }
  
  _combinations = cliqueCount;
  
  _mapping = StlIntVector(_combinations, 0);
  for (int i = 1; i < _combinations; ++i)
  {
    _mapping[i] = _mapping[i-1] + 1;
  }
//  
//  std::mt19937 rng(0);
//  std::shuffle(_mapping.begin(), _mapping.end(), rng);
}
  
void CompatibilityGraph::write(std::ostream& out) const
{
  out << lemon::countEdges(_G) << " #compatibilities" << std::endl;
  for (EdgeIt e_csdt(_G); e_csdt != lemon::INVALID; ++e_csdt)
  {
    Node v_cs = _G.u(e_csdt);
    Node v_dt = _G.v(e_csdt);
    
    const IntPair& cs = _nodeToCharStateTree[v_cs];
    const IntPair& dt = _nodeToCharStateTree[v_dt];
    
    out << "(" << cs.first << "," << cs.second << ") ("
        << dt.first << "," << dt.second << ")" << std::endl;
  }
  
  out << _cliques.front().size() << " #clique size" << std::endl;
  out << _cliques.size() << " #cliques" << std::endl;
  
  bool first = true;
  
  for (const NodeVector& C : _cliques)
  {
    first = true;
    for (Node v_cs : C)
    {
      if (first)
        first = false;
      else
        out << " ";
      
      const IntPair& cs = _nodeToCharStateTree[v_cs];
      out << "(" << cs.first << "," << cs.second << ")";
    }
    out << std::endl;
  }
}
  
void CompatibilityGraph::init(unsigned long combination,
                              RealTensor& F,
                              StateTreeVector& S,
                              RealTensor& F_lb,
                              RealTensor& F_ub,
                              StlIntVector& mapNewCharToOldChar,
                              StlIntVector& mapOldCharToNewChar) const
{
  assert(combination < _combinations);
  
  const int n =  _cliques.front().size();
  const int m = _M.m();
  const int k = _M.k();
  
  mapNewCharToOldChar = StlIntVector(n, -1);
  mapOldCharToNewChar = StlIntVector(_M.n(), -1);
  
  // construct F and S
  F = RealTensor(k, m, n);
  F_lb = RealTensor(k, m, n);
  F_ub = RealTensor(k, m, n);
  for (int p = 0; p < m; ++p)
  {
    F.setRowLabel(p, _M(p, 0).sampleLabel());
  }
  S.clear();
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Initializing state tree combination # " << _mapping[combination] << " ..." << std::endl;
  }
  
  const NodeVector& clique = _cliques[_mapping[combination]];
  int cc = 0; // mapped c
  for (NodeVectorIt it = clique.begin(); it != clique.end(); ++it, ++cc)
  {
    Node v_cs = *it;
    const IntPair& cs = _nodeToCharStateTree[v_cs];
    mapNewCharToOldChar[cc] = cs.first;
    mapOldCharToNewChar[cs.first] = cc;
    
    S.push_back(_M.stateTree(cs.second, cs.first));
    F.setColLabel(cc, _M(0, cs.first).characterLabel());
    
    for (int i = 0; i < k; ++i)
    {
      for (int p = 0; p < m; ++p)
      {
        // I'm specifically using fLB here, point estimate VAF may not be feasible for the current state tree...
        double f_p_ci = _M.get_f_lb(cs.second, p, cs.first, i);
        double l_p_ci = f_p_ci;
        double u_p_ci =  _M.get_f_ub(cs.second, p, cs.first, i);
        assert(l_p_ci <= u_p_ci);
        
        F.set(i, p, cc, f_p_ci);
        F_lb.set(i, p, cc, l_p_ci);
        F_ub.set(i, p, cc, u_p_ci);
      }
    }
  }
}
  
void CompatibilityGraph::applyFilter(const IntPairSet& filter,
                                     const NodeMatrix& cliques)
{
  // Map filter to nodes
  NodeSet nodeFilter;
  for (const IntPair& cs : filter)
  {
    nodeFilter.insert(_charStateTreeToNode[cs.first][cs.second]);
  }
  
  for (const auto& C : cliques)
  {
    NodeSet CC(C.begin(), C.end());
    NodeSet X;
    std::set_intersection(nodeFilter.begin(), nodeFilter.end(), CC.begin(), CC.end(), std::inserter(X, X.begin()));
    
    if (X == nodeFilter)
    {
      _cliques.push_back(C);
    }
  }
}
  
void CompatibilityGraph::init(const IntPairSet& filter, int size)
{
  initEdges();

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Searching for maximum cliques ... " << std::endl;
  }
  
  BronKerbosch bk(_G, _cliqueLimit);
  bk.run(BronKerbosch::BK_PIVOT_DEGENERACY);
  bk.sortBySize();
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    IntSet cliqueSizes = bk.getCliqueSizes();
    for (int s : cliqueSizes)
    {
      std::cerr << "Size: " << s << "; #maximal cliques: " << bk.getNrCliquesBySize(s) << std::endl;
    }
    std::cerr << "Total number of maximal cliques: " << bk.getNumberOfMaximalCliques() << std::endl;
  }
  
  NodeMatrix cliques;
  if (size == -1)
  {
    bk.getMaximumCliques(cliques);
    
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cerr << "Found " << cliques.size() << " maximum cliques of size " << bk.getMaximumCliqueSize() << std::endl;
    }
  }
  else
  {
    bk.getCliquesBySize(size, cliques);
    
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cerr << "Found " << cliques.size() << " maximal cliques of size " << size << std::endl;
    }
  }
  
  if (filter.empty())
  {
    std::swap(_cliques, cliques);
  }
  else
  {
    applyFilter(filter, cliques);
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cerr << "Applying filter results in "
                << _cliques.size() << " maximal cliques of size " << size << std::endl;
    }
  }
  _combinations = _cliques.size();
  
  _mapping = StlIntVector(_combinations, 0);
  for (int i = 1; i < _combinations; ++i)
  {
    _mapping[i] = _mapping[i-1] + 1;
  }
//  
//  std::mt19937 rng(0);
//  std::shuffle(_mapping.begin(), _mapping.end(), rng);
  
  

//  bool first = true;
//  for (const NodeVector& C : bk.getMaxCliques())
//  {
//    first = true;
//    for (Node v_cs : C)
//    {
//      if (first)
//        first = false;
//      else
//        std::cerr << " ";
//      
//      const IntPair& cs = _nodeToCharStateTree[v_cs];
//      std::cerr << "(" << cs.first << "," << cs.second << ")";
//    }
//    std::cerr << std::endl;
//  }
}
  
void CompatibilityGraph::initVertices()
{
  const int n = _M.n();

  _charStateTreeToNode = NodeMatrix(n);
  for (int c = 0; c < n; ++c)
  {
    int stateTreeCount = _M.numStateTrees(c);
    _charStateTreeToNode[c] = NodeVector(stateTreeCount, lemon::INVALID);
    
    for (int s = 0; s < stateTreeCount; ++s)
    {
      Node v_cs = _G.addNode();
      _charStateTreeToNode[c][s] = v_cs;
      _nodeToCharStateTree[v_cs] = std::make_pair(c, s);
    }
  }
}
  
void CompatibilityGraph::initEdges()
{
  int conflict = 0;
  for (NodeIt v_cs(_G); v_cs != lemon::INVALID; ++v_cs)
  {
    for (NodeIt v_dt(v_cs); v_dt != lemon::INVALID; ++v_dt)
    {
      if (v_dt == v_cs) continue;
      
      const IntPair& cs = _nodeToCharStateTree[v_cs];
      const IntPair& dt = _nodeToCharStateTree[v_dt];
      
      if (cs.first == dt.first) continue;
      
      if (isCompatible(cs, dt))
      {
        _G.addEdge(v_cs, v_dt);
      }
      else
      {
        ++conflict;
        if (g_verbosity >= VERBOSE_DEBUG)
        {
          std::cerr << "Conflict between (" << cs.first << "," << cs.second
                    << ") and (" << dt.first << "," << dt.second << ")" << std::endl;
        }
      }
    }
  }
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Number of conflicts detected: " << conflict << std::endl;
    std::cerr << "Compatibility graph has " << lemon::countNodes(_G) << " nodes and " << lemon::countEdges(_G) << " edges" << std::endl;
  }
}
  
void CompatibilityGraph::writeDOT(std::ostream& out) const
{
  out << "graph G {" << std::endl;
  
  for (NodeIt v_cs(_G); v_cs != lemon::INVALID; ++v_cs)
  {
    const IntPair& cs = _nodeToCharStateTree[v_cs];
    StateTree S = _M.stateTree(cs.second, cs.first);
    
    out << "\t" << _G.id(v_cs) << " [label=\"(" << cs.first << "," << cs.second << ")";
    
    for (StateTree::ArcIt a_ij(S.S()); a_ij != lemon::INVALID; ++a_ij)
    {
      StateTree::Node v_i = S.S().source(a_ij);
      StateTree::Node v_j = S.S().target(a_ij);
      
      out << "\\n" << S.label(S.state(v_i)) << " -> " << S.label(S.state(v_j));
    }
    out << "\"]" << std::endl;
  }
  
  for (EdgeIt e_csdt(_G); e_csdt != lemon::INVALID; ++e_csdt)
  {
    Node v_cs = _G.u(e_csdt);
    Node v_dt = _G.v(e_csdt);
    
    out << "\t" << _G.id(v_cs) << " -- " << _G.id(v_dt) << std::endl;
  }
  
  out << "}" << std::endl;
}
  
bool CompatibilityGraph::isCompatible(const IntPair& cs, const IntPair& dt) const
{
  const int m = _M.m();
  const int k = _M.k();
  
  if (!_M.intervals().empty() && _M.interval(cs.first) == _M.interval(dt.first))
  {
    CharacterMatrix::IntPairPairSet C_c = _M.copyTree(cs.second, cs.first);
    CharacterMatrix::IntPairPairSet C_d = _M.copyTree(dt.second, dt.first);
    
    if (C_c != C_d)
    {
      // incompatible copy-trees
      return false;
    }
  }
  
  std::vector<StateTree> S;
  S.push_back(_M.stateTree(cs.second, cs.first));
  S.push_back(_M.stateTree(dt.second, dt.first));
  
  RealTensor F_lb(k, m, 2);
  RealTensor F_ub(k, m, 2);

  for (int p = 0; p < m; ++p)
  {
//    F_lb.setColLabel(0, _M(p, cs.first).characterLabel());
//    F_lb.setColLabel(1, _M(p, dt.first).characterLabel());
    for (int i = 0; i < k; ++i)
    {
      F_lb.set(i, p, 0, _M.get_f_lb(cs.second, p, cs.first, i));
      F_ub.set(i, p, 0, _M.get_f_ub(cs.second, p, cs.first, i));

      F_lb.set(i, p, 1, _M.get_f_lb(dt.second, p, dt.first, i));
      F_ub.set(i, p, 1, _M.get_f_ub(dt.second, p, dt.first, i));
    }
  }
  
  RootedCladisticNoisyAncestryGraph G(F_lb, S, F_lb, F_ub);
  G.init();
//  G.setLabels(F_lb);
//  G.writeDOT(std::cout);
  
  RootedCladisticNoisyEnumeration enumerate(G, -1, -1, 1, 2, false, false, IntSet());
  enumerate.run();
  
  if (enumerate.objectiveValue() == 2)
  {
    return true;
  }
  else
  {
    return false;
  }
}
  
} // namespace gm
