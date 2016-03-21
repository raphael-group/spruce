/*
 *  solutionset.h
 *
 *   Created on: 1-oct-2015
 *       Author: M. El-Kebir
 */

#ifndef SOLUTIONSET_H
#define SOLUTIONSET_H

#include "utils.h"
#include "solution.h"
#include "realtensor.h"
#include "statetree.h"
#include "perfectphylotree.h"
#include "stategraph.h"

namespace gm {
  
class SolutionSet
{
public:
  typedef std::vector<Solution> SolutionVector;
  typedef SolutionVector::const_iterator SolutionVectorIt;
  typedef SolutionVector::iterator SolutionVectorNonConstIt;
  typedef std::vector<StateTree> StateTreeVector;
  
  SolutionSet();
  
  int solutionCount() const
  {
    return _sol.size();
  }
  
  const Solution& solution(int idx) const
  {
    assert(0 <= idx && idx < _sol.size());
    
    return _sol[idx];
  }
  
  void clear()
  {
    _sol.clear();
  }
  
  void add(const Solution& sol);
  
  void add(const SolutionSet& sols);
  
  void writeSummaryJSON(std::ostream& out,
                        const std::map<std::string, std::string>& label2color);
  
  void writeSummaryDOT(std::ostream& out,
                       const std::map<std::string, std::string>& label2color);
  
  void initSummaryGraph(Digraph& G,
                        Digraph::Node& root,
                        Digraph::NodeMap<std::string>& label,
                        Digraph::NodeMap<StateGraph::CnaTriple>& state);
  
  void writeSummaryDOT(std::ostream& out, const PerfectPhyloTree& trueT) const;
  
  void sort();
  
  int unique();
  
  void assignDistancesByOccurenceCounts();
  
  friend std::ostream& operator<<(std::ostream& out, const SolutionSet& sols);
  friend std::istream& operator>>(std::istream& in, SolutionSet& sols);
  
private:
  typedef std::map<std::string, std::map<std::string, int > > Map;
  typedef std::map<std::string, int> StringIntMap;
  typedef StringIntMap::const_iterator StringIntMapIt;
  
private:
  SolutionVector _sol;
  Map _occArc;
  StringIntMap _nodes;

  void initSummary();
  
  struct Compare
  {
    bool operator()(const Solution& sol1,
                    const Solution& sol2)
    {
      return sol1.distance() < sol2.distance();
    }
  };
};
  
std::ostream& operator<<(std::ostream& out, const SolutionSet& sols);
std::istream& operator>>(std::istream& in, SolutionSet& sols);
  
} // namespace gm

#endif // SOLUTIONSET_H