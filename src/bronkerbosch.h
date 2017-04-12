/*
 * bronkerbosch.h
 *
 *  Created on: 21-jan-2014
 *      Author: M. El-Kebir
 */

#ifndef BRONKERBOSCH_H
#define BRONKERBOSCH_H

#include <iostream>
#include <vector>
#include <list>
#include <limits>
#include <set>
#include <boost/dynamic_bitset.hpp>
#include "utils.h"

namespace gm {
  
class BronKerbosch
{
public:
  GRAPH_TYPEDEFS(Graph);
  
  typedef enum
  {
    BK_CLASSIC,
    BK_PIVOT,
    BK_PIVOT_DEGENERACY
  } SolverType;

  typedef std::vector<Node> NodeVector;
  typedef std::map<int, int> IntIntMap;
  
protected:
  typedef Graph::NodeMap<size_t> BitNodeMap;
  typedef boost::dynamic_bitset<> BitSet;
  typedef Graph::NodeMap<BitSet> BitSetNodeMap;
  typedef std::list<Node> NodeList;
  typedef std::vector<NodeList> NodeListVector;
  typedef std::vector<NodeVector> NodeVectorVector;
  typedef NodeVectorVector::const_iterator NodeVectorVectorIt;
  typedef std::map<int, NodeVectorVectorIt> CliqueItMap;
  typedef std::set<int> IntSet;
  
  struct CliqueComp
  {
    bool operator()(const NodeVector& a, const NodeVector& b) const
    {
      return a.size() > b.size();
    }
  };

public:
  BronKerbosch(const Graph& g, int limit);

  void run(SolverType type);

  void print(std::ostream& out) const;

  size_t getNumberOfMaximalCliques() const { return _cliques.size(); }
  
  size_t getMaximumCliqueSize() const
  {
    if (_cliques.empty())
    {
      return 0;
    }
    else
    {
      return _cliques[0].size();
    }
  }
  
  const IntIntMap& getNrCliquesBySize() const
  {
    return _nrCliquesBySize;
  }
  
  int getNrCliquesBySize(int size) const
  {
    if (_nrCliquesBySize.count(size) == 0)
    {
      return 0;
    }
    else
    {
      return _nrCliquesBySize.find(size)->second;
    }
  }
  
  void getCliquesBySize(int size,
                        std::vector<NodeVector>& cliqueVector) const
  {
    int count = getNrCliquesBySize(size);
    if (count == 0)
      return;
    
    assert(_firstCliqueBySize.find(size) != _firstCliqueBySize.end());
    NodeVectorVectorIt firstIt = _firstCliqueBySize.find(size)->second;
    
    cliqueVector.insert(cliqueVector.end(), firstIt, firstIt + count);
  }
  
  IntSet getCliqueSizes() const
  {
    IntSet res;
    
    for (const auto& kv : _nrCliquesBySize)
    {
      res.insert(kv.first);
    }
    
    return res;
  }
  
  void getMaximumCliques(std::vector<NodeVector>& cliqueVector) const
  {
    int size = _cliques[0].size();
    getCliquesBySize(size, cliqueVector);
  }
  
  size_t getNumberOfMaximumCliques() const
  {
    size_t maximumCliqueSize = getMaximumCliqueSize();
    for (size_t i = 0; i < _cliques.size(); ++i)
    {
      if (_cliques[i].size() != maximumCliqueSize)
      {
        return i;
      }
    }
    return _cliques.size();
  }

  const NodeVectorVector& getCliques() const { return _cliques; }
  
  void sortBySize()
  {
    std::sort(_cliques.begin(), _cliques.end(), CliqueComp());
    
    if (_cliques.empty())
      return;
    
    _nrCliquesBySize.clear();
    for (size_t i = 0; i < _cliques.size(); ++i)
    {
      int s = _cliques[i].size();
      if (_nrCliquesBySize.count(s) == 0)
        _nrCliquesBySize[s] = 0;
      
      ++_nrCliquesBySize[s];
    }
    
    _firstCliqueBySize.clear();
    for (size_t i = 0; i < _cliques.size(); ++i)
    {
      int s = _cliques[i].size();
      if (i == 0)
      {
        _firstCliqueBySize[s] = _cliques.begin() + i;
      }
      else if (_cliques[i-1].size() != s)
      {
        _firstCliqueBySize[s] = _cliques.begin() + i;
      }
    }
  }

protected:
  const Graph& _g;
  const size_t _n;
  const int _limit;
  NodeVectorVector _cliques;
  NodeVector _bitToNode;
  BitNodeMap _nodeToBit;
  BitSetNodeMap _bitNeighborhood;
  IntIntMap _nrCliquesBySize;
  CliqueItMap _firstCliqueBySize;

  size_t computeDegeneracy(NodeList& order);

  void report(const BitSet& R);

  void printBitSet(const BitSet& S, std::ostream& out) const;

private:
  /// Classic Bron-Kerbosch algorithm without pivoting
  ///
  /// Reports maximal cliques in P \cup R (but not in X)
  void bkClassic(BitSet P, BitSet R, BitSet X);
  void bkPivot(BitSet P, BitSet R, BitSet X);
  void bkDegeneracy(const NodeList& order);
};

} // namespace gm

#endif // BRONKERBOSCH_H
