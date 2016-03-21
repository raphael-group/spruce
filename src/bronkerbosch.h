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
  
protected:
  typedef Graph::NodeMap<size_t> BitNodeMap;
  typedef boost::dynamic_bitset<> BitSet;
  typedef Graph::NodeMap<BitSet> BitSetNodeMap;
  typedef std::list<Node> NodeList;
  typedef std::vector<NodeList> NodeListVector;
  
  struct CliqueComp
  {
    bool operator()(const NodeVector& a, const NodeVector& b) const
    {
      return a.size() > b.size();
    }
  };

public:
  BronKerbosch(const Graph& g);

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
  
  void getCliques(const size_t count,
                  std::vector<NodeVector>& cliqueVector) const
  {
    assert(0 <= count && count <= _cliques.size());
    cliqueVector.insert(cliqueVector.end(), _cliques.begin(), _cliques.begin() + count);
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

  const std::vector<NodeVector>& getMaxCliques() const { return _cliques; }
  
  void sortBySize()
  {
    std::sort(_cliques.begin(), _cliques.end(), CliqueComp());
  }

protected:
  const Graph& _g;
  const size_t _n;
  std::vector<NodeVector> _cliques;
  std::vector<Node> _bitToNode;
  BitNodeMap _nodeToBit;
  BitSetNodeMap _bitNeighborhood;

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
