/*
 *  solutionset.cpp
 *
 *   Created on: 1-oct-2015
 *       Author: M. El-Kebir
 */

#include "solutionset.h"
#include "rootedcladisticancestrygraph.h"
#include "perfectphylotree.h"
#include <boost/unordered_set.hpp>

namespace gm {
  
SolutionSet::SolutionSet()
  : _sol()
{
}
  
void SolutionSet::sort()
{
  std::sort(_sol.begin(), _sol.end(), Compare());
}
  
void SolutionSet::add(const SolutionSet& sols)
{
  int n = sols.solutionCount();
  for (int i = 0; i < n; ++i)
  {
    add(sols.solution(i));
  }
}

void SolutionSet::add(const Solution& sol)
{
  _sol.push_back(sol);
}
  
void SolutionSet::writeSummaryDOT(std::ostream& out,
                                  const PerfectPhyloTree& trueT) const
{
  typedef std::map<std::string, bool> BoolMap;
  typedef std::map<std::string, int> IntMap;
  typedef IntMap::const_iterator IntMapIt;
  
  typedef std::vector<BoolMap> BoolMapVector;
  typedef std::vector<IntMap> IntMapVector;
  
  typedef std::map<std::string, BoolMapVector> BoolMapVectorMap;
  typedef BoolMapVectorMap::const_iterator BoolMapVectorMapIt;
  typedef std::map<std::string, IntMapVector> IntMapVectorMap;
  typedef IntMapVectorMap::const_iterator IntMapVectorMapIt;
  
  typedef std::vector<BoolMapVectorMap> BoolMapVectorMapVector;
  typedef std::vector<IntMapVectorMap> IntMapVectorMapVector;
  
  typedef std::vector<StlBoolMatrix> StlBool3Matrix;
  typedef std::vector<StlBool3Matrix> StlBool4Matrix;
  typedef std::vector<StlIntMatrix> StlInt3Matrix;
  typedef std::vector<StlInt3Matrix> StlInt4Matrix;
  typedef std::map<IntPair, std::string> IntPairMap;
  typedef IntPairMap::const_iterator IntPairMapIt;
  
  assert(!_sol.empty());
  
  const int n = _sol[0].inferredF().n();
  const int k = _sol[0].inferredF().k();
  
  // trueArc[c][i][d][j] = true if (c,i) -> (d,j) in true solution
  BoolMapVectorMapVector trueArc(n);
  // occArc[c][i][d][j] is #times (c,i) -> (d,j) occurs
  IntMapVectorMapVector occArc(n);
  
  static int minPenwidth = 1;
  
  IntPairMap nodes;
  IntMapVector stateIdx(n);

  // populate occArc and nodes
  for (int solIdx = 0; solIdx < solutionCount(); ++solIdx)
  {
    const Solution& sol = _sol[solIdx];
    
    PerfectPhyloTree T(sol.A(), sol.S());
    T.setLabels(sol.inferredF());
    
    for (Digraph::ArcIt a_cidj(T.T()); a_cidj != lemon::INVALID; ++a_cidj)
    {
      Digraph::Node v_ci = T.T().source(a_cidj);
      Digraph::Node v_dj = T.T().target(a_cidj);
      
      const IntPair& ci = T.nodeToCharState(v_ci);
      const IntPair& dj = T.nodeToCharState(v_dj);
      
      // assuming that states are the same!
      const std::string& label_ci = T.label(v_ci);
      const std::string& label_dj = T.label(v_dj);
      
      int real_c = -1;
      if (label_ci.length() >= 2 && label_ci[1] == '*')
      {
        real_c = 0;
      }
      else
      {
        sscanf(label_ci.c_str(), "(%d", &real_c);
      }
      
      int real_d = -1;
      sscanf(label_dj.c_str(), "(%d", &real_d);
      
      assert(0 <= real_c && real_c < n);
      assert(0 <= real_d && real_d < n);
      
      const std::string& label_i = sol.S()[ci.first].label(ci.second);
      const std::string& label_j = sol.S()[dj.first].label(dj.second);
      
      if (occArc[real_c].find(label_i) == occArc[real_c].end())
      {
        occArc[real_c][label_i] = IntMapVector(n);
      }
      ++occArc[real_c][label_i][real_d][label_j];

      nodes[std::make_pair(real_c, ci.second)] = label_ci;
      nodes[std::make_pair(real_d, dj.second)] = label_dj;
      stateIdx[real_c][label_i] = ci.second;
      stateIdx[real_d][label_j] = dj.second;
    }
  }
  
  // populate trueArc
  for (Digraph::ArcIt a_cidj(trueT.T()); a_cidj != lemon::INVALID; ++a_cidj)
  {
    Digraph::Node v_ci = trueT.T().source(a_cidj);
    Digraph::Node v_dj = trueT.T().target(a_cidj);
    
    const IntPair& ci = trueT.nodeToCharState(v_ci);
    const IntPair& dj = trueT.nodeToCharState(v_dj);
    
    const std::string& label_i = trueT.S(ci.first).label(ci.second);
    const std::string& label_j = trueT.S(dj.first).label(dj.second);
    
    if (trueArc[ci.first].find(label_i) == trueArc[ci.first].end())
    {
      trueArc[ci.first][label_i] = BoolMapVector(n);
    }
    trueArc[ci.first][label_i][dj.first][label_j] = true;
  }

//  const int m = _F.m();
  out << "digraph G {" << std::endl;
  out.precision(3);
  
  for (IntPairMapIt it = nodes.begin(); it != nodes.end(); ++it)
  {
    const IntPair& ci = it->first;

    out << "\t" << ci.first * n + ci.second << " [label=\"" << it->second << "\"]" << std::endl;
  }
  
  for (int c = 0; c < n; ++c)
  {
    const IntMapVectorMap& occArc_c = occArc[c];
    for (IntMapVectorMapIt it_i = occArc_c.begin(); it_i != occArc_c.end(); ++it_i)
    {
      const std::string& label_i = it_i->first;
      const IntMapVector& occArc_ci = it_i->second;
      for (int d = 0; d < n; ++d)
      {
        const IntMap& occArc_ci_d = occArc_ci[d];
        for (IntMapIt it_j = occArc_ci_d.begin(); it_j != occArc_ci_d.end(); ++it_j)
        {
          const std::string& label_j = it_j->first;
          int occ = it_j->second;

          if (occ > 0)
          {
            out << "\t" << c*n + stateIdx[c][label_i] << " -> " << d*n + stateIdx[d][label_j]
                << " [label=" << occ << ",penwidth="
                << minPenwidth + 5 * (occ / (double) solutionCount());

            if (trueArc[c].find(label_i) != trueArc[c].end()
                && trueArc[c][label_i][d].find(label_j) != trueArc[c][label_i][d].end())
            {
              out << ",color=red";
            }            
            out << "]" << std::endl;
          }
        }
      }
    }
  }
  out << "}" << std::endl;
}
  
void SolutionSet::initSummary()
{
  _occArc.clear();
  _nodes.clear();
  
  // populate occArc and nodes
  int nodeIdx = 0;
  for (int solIdx = 0; solIdx < solutionCount(); ++solIdx)
  {
    const Solution& sol = _sol[solIdx];
    
    PerfectPhyloTree T(sol.A(), sol.S());
    T.setLabels(sol.inferredF());
    
    for (Digraph::ArcIt a_cidj(T.T()); a_cidj != lemon::INVALID; ++a_cidj)
    {
      Digraph::Node v_ci = T.T().source(a_cidj);
      Digraph::Node v_dj = T.T().target(a_cidj);
      
      // assuming that labels are ids!
      const std::string& label_ci = T.label(v_ci);
      const std::string& label_dj = T.label(v_dj);
      
      if (_occArc[label_ci].find(label_dj) == _occArc[label_ci].end())
      {
        _occArc[label_ci][label_dj] = 0;
      }
      ++_occArc[label_ci][label_dj];
      
      if (_nodes.find(label_ci) == _nodes.end())
      {
        _nodes[label_ci] = nodeIdx++;
      }
      if (_nodes.find(label_dj) == _nodes.end())
      {
        _nodes[label_dj] = nodeIdx++;
      }
    }
  }
}
  
void SolutionSet::initSummaryGraph(Digraph& G,
                                   Digraph::Node& root,
                                   Digraph::NodeMap<std::string>& label,
                                   Digraph::NodeMap<StateGraph::CnaTriple>& state)
{
  G.clear();
  initSummary();
  
  // init nodes
  Digraph::NodeMap<int> id(G);
  std::vector<Digraph::Node> id2node(_nodes.size());
  
  for (StringIntMapIt it = _nodes.begin(); it != _nodes.end(); ++it)
  {
    Digraph::Node v_ci = G.addNode();
    id[v_ci] = it->second;
    id2node[it->second] = v_ci;
    label[v_ci] = it->first;
    
    StateGraph::CnaTriple xyz;
    sscanf(&(it->first.c_str()[it->first.length() - 8]), "(%d,%d,%d))", &xyz._x, &xyz._y, &xyz._z);
    
    state[v_ci] = xyz;
    
//    std::cout << it->first << " " << it->second << " " << xyz._x << "," << xyz._y << "," << xyz._z << std::endl;
  }
  
  // init edges
  for (StringIntMapIt it1 = _nodes.begin(); it1 != _nodes.end(); ++it1)
  {
    for (StringIntMapIt it2 = _nodes.begin(); it2 != _nodes.end(); ++it2)
    {
      if (it1 == it2) continue;
      if (_occArc.find(it1->first) != _occArc.end() &&
          _occArc[it1->first].find(it2->first) != _occArc[it1->first].end() &&
          _occArc[it1->first][it2->first] > 0)
      {
        G.addArc(id2node[it1->second], id2node[it2->second]);
      }
    }
  }
}
  
void SolutionSet::writeSummaryDOT(std::ostream& out,
                                  const std::map<std::string, std::string>& label2color)
{
  typedef std::vector<StlBoolMatrix> StlBool3Matrix;
  typedef std::vector<StlBool3Matrix> StlBool4Matrix;
  typedef std::vector<StlIntMatrix> StlInt3Matrix;
  typedef std::vector<StlInt3Matrix> StlInt4Matrix;
  typedef std::set<IntPair> IntPairSet;
  typedef IntPairSet::const_iterator IntPairSetIt;
  
  assert(!_sol.empty());
  initSummary();
  static int minPenwidth = 1;
  
  //  const int m = _F.m();
  out << "digraph G {" << std::endl;
  out.precision(3);
  
  for (StringIntMapIt it = _nodes.begin(); it != _nodes.end(); ++it)
  {
    std::string character = it->first.substr(1, it->first.find(",") - 1);
    const auto& it2 = label2color.find(character);
    out << "\t" << it->second
        << " [label=\"" << it->first << "\",shape=box,style=rounded,fixedsize=true,width=2.5";
    bool isCNA = it->first[it->first.length() - 3] == '0';
    if (it2 != label2color.end() && !isCNA)
    {
      out << ",penwidth=5,color=" << it2->second;
    }
    out << "]" << std::endl;
  }
  
  for (StringIntMapIt it1 = _nodes.begin(); it1 != _nodes.end(); ++it1)
  {
    for (StringIntMapIt it2 = _nodes.begin(); it2 != _nodes.end(); ++it2)
    {
      if (it1 == it2) continue;
      if (_occArc.find(it1->first) != _occArc.end() &&
          _occArc[it1->first].find(it2->first) != _occArc[it1->first].end() &&
          _occArc[it1->first][it2->first] > 0)
      {
        int occ = _occArc[it1->first][it2->first];
        out << "\t" << it1->second << " -> " << it2->second
            << " [label=" << occ << ",penwidth="
            << minPenwidth + 5 * (occ / (double) solutionCount());
        out << "]" << std::endl;
        
      }
    }
  }
  out << "}" << std::endl;
}
  
void SolutionSet::writeSummaryJSON(std::ostream& out,
                                   const std::map<std::string, std::string>& label2color)
{
  typedef std::vector<StlBoolMatrix> StlBool3Matrix;
  typedef std::vector<StlBool3Matrix> StlBool4Matrix;
  typedef std::vector<StlIntMatrix> StlInt3Matrix;
  typedef std::vector<StlInt3Matrix> StlInt4Matrix;
  typedef std::set<IntPair> IntPairSet;
  typedef IntPairSet::const_iterator IntPairSetIt;
  
  assert(!_sol.empty());
  initSummary();
  
  //  const int m = _F.m();
  out << "{" << std::endl;
  out << "\t\"nodes\": [" << std::endl;
  out.precision(3);

  for (StringIntMapIt it = _nodes.begin(); it != _nodes.end(); ++it)
  {
    if (it != _nodes.begin())
      out << "," << std::endl;
//    out << "\t\t\t\"" << it->second << "\": [\n\t\t\t\t{" << std::endl << "\t\t\t\t\t\"label\": \"" << it->first << "\"\n\t\t\t\t}\n\t\t\t]";
      out << "\t\t{\n" << "\t\t\t\"id\": " << it->second << ",\n\t\t\t\"label\": \"" << it->first << "\"\n\t\t}";
  }
  out << std::endl << "\t]," << std::endl;
  
  for (int solIdx = 0; solIdx < solutionCount(); ++solIdx)
  {
    if (solIdx != 0)
      out << "," << std::endl;
    out << "\t\"sol_" << solIdx << "\": [" << std::endl;
    const Solution& sol = _sol[solIdx];
    const int m = sol.observedF().m();
    const int n = sol.observedF().n();
    
    PerfectPhyloTree T(sol.A(), sol.S());
    T.setLabels(sol.inferredF());
    
//    // write usages
//    for (Digraph::NodeIt v_ci(T.T()); v_ci != lemon::INVALID; ++v_ci)
//    {
//      const std::string& label_ci = T.label(v_ci);
//      const IntPair& ci = T.nodeToCharState(v_ci);
//      
//      out << "    {\"node_id\":" << _nodes.find(label_ci)->second
//          << ",\"U\":[";
//      
//      bool first = true;
//      for (int p = 0; p < m; ++p)
//      {
//        if (first)
//          first = false;
//        else
//          out << ",";
//        
//        out << "{\"sample_id\":" << p << ",\"sample_label\":" << '"' << sol.observedF().getRowLabel(p) << '"';
//        
//        int row_idx = ci.first == 0 && ci.second == 0 ? 0 : n * (ci.second - 1) + ci.first + 1;
//        double u_pci = sol.U()(p, row_idx);
//        
//        out << ",\"usage\":" << u_pci;
//        out << "}";
//      }
//      
//      out << "]"
//          << "}" << std::endl;
//    }
    
    bool first = true;
    for (Digraph::ArcIt a_cidj(T.T()); a_cidj != lemon::INVALID; ++a_cidj)
    {
      if (first)
        first = false;
      else
        out << "," << std::endl;
      
      Digraph::Node v_ci = T.T().source(a_cidj);
      Digraph::Node v_dj = T.T().target(a_cidj);
      
      // assuming that labels are ids!
      const std::string& label_ci = T.label(v_ci);
      const std::string& label_dj = T.label(v_dj);
      
      out << "\t\t{\n\t\t\t\"source\": " << _nodes.find(label_ci)->second
          << ",\n\t\t\t\"target\": " << _nodes.find(label_dj)->second << "\n\t\t}";
    }
    out << std::endl << "\t]";
  }
  out << std::endl << "}" << std::endl;
}
  
void SolutionSet::assignDistancesByOccurenceCounts()
{
  initSummary();
  
  int solCount = _sol.size();
  for (int idx = 0; idx < solCount; ++idx)
  {
    Solution& sol = _sol[idx];
    
    int distance = 0;
    PerfectPhyloTree T(sol.A(), sol.S());
    T.setLabels(sol.inferredF());
    for (Digraph::ArcIt a_cidj(T.T()); a_cidj != lemon::INVALID; ++a_cidj)
    {
      Digraph::Node v_ci = T.T().source(a_cidj);
      Digraph::Node v_dj = T.T().target(a_cidj);
      
      // assuming that labels are ids!
      const std::string& label_ci = T.label(v_ci);
      const std::string& label_dj = T.label(v_dj);
      
      distance -= _occArc[label_ci][label_dj];
    }
    
    sol.distance() = distance;
  }
}
  
int SolutionSet::unique()
{
  typedef std::pair<IntPair, IntPair> IntPairPair;
  typedef std::set<IntPairPair> ArcSet;
  typedef std::set<ArcSet> ArcSetSet;
  
  int res = 0;
  ArcSetSet AS;
  for (SolutionVectorNonConstIt it = _sol.begin(); it != _sol.end();)
  {
    const Solution& sol = *it;
    PerfectPhyloTree T(sol.A(), sol.S());
    const Digraph& TT = T.T();
    
    ArcSet arcs;
    for (Digraph::ArcIt a_cidj(TT); a_cidj != lemon::INVALID; ++a_cidj)
    {
      Digraph::Node v_ci = TT.source(a_cidj);
      Digraph::Node v_dj = TT.target(a_cidj);
      
      const IntPair& ci = T.nodeToCharState(v_ci);
      const IntPair& dj = T.nodeToCharState(v_dj);
      
      arcs.insert(std::make_pair(ci, dj));
    }
    
    if (AS.find(arcs) == AS.end())
    {
      AS.insert(arcs);
      ++it;
    }
    else
    {
      it = _sol.erase(it);
      ++res;
    }
  }
  return res;
}
  
std::ostream& operator<<(std::ostream& out, const SolutionSet& sols)
{
  int solCount = sols.solutionCount();
  out << solCount << " # solutions" << std::endl << std::endl;
  
  for (int idx = 0; idx < solCount; ++idx)
  {
    out << sols.solution(idx);
  }
  
  return out;
}
  
std::istream& operator>>(std::istream& in, SolutionSet& sols)
{
  int solCount = -1;
  
  std::string line;
  gm::getline(in, line);
  std::stringstream ss(line);
  ss >> solCount;
  
  gm::getline(in, line);
  
  if (solCount <= 0)
  {
    throw std::runtime_error(getLineNumber() + "Error: solCount should be nonnegative");
  }
  
  sols._sol.clear();
  for (int idx = 0; idx < solCount; ++idx)
  {
    Solution sol;
    in >> sol;
    sols.add(sol);
    
    gm::getline(in, line);
  }
  
  return in;
}
  
} // namespace gm
