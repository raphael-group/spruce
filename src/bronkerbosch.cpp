/*
 * bronkerbosch.cpp
 *
 *  Created on: 21-jan-2014
 *      Author: M. El-Kebir
 */

#include "bronkerbosch.h"

namespace gm {

BronKerbosch::BronKerbosch(const Graph& g,
                           int limit = -1)
  : _g(g)
  , _n(static_cast<size_t>(lemon::countNodes(_g)))
  , _limit(limit)
  , _cliques()
  , _bitToNode()
  , _nodeToBit(g, std::numeric_limits<size_t>::max())
  , _bitNeighborhood(g, BitSet(_n))
{
  // initialize mappings
  _bitToNode.reserve(_n);
  size_t i = 0;
  for (NodeIt v(_g); v != lemon::INVALID; ++v, ++i)
  {
    _bitToNode.push_back(v);
    _nodeToBit[v] = i;
  }

  // initialize neighborhoods
  for (NodeIt v(_g); v != lemon::INVALID; ++v, ++i)
  {
    BitSet& neighborhood = _bitNeighborhood[v];
    for (IncEdgeIt e(_g, v); e != lemon::INVALID; ++e)
    {
      Node w = _g.oppositeNode(v, e);
      neighborhood[_nodeToBit[w]] = 1;
    }
  }
}

size_t BronKerbosch::computeDegeneracy(NodeList& order)
{
  // Requires O(|V| + |E|) time
  order.clear();

  typedef Graph::NodeMap<size_t> DegNodeMap;
  typedef Graph::NodeMap<NodeList::iterator> NodeListItMap;

  BoolNodeMap present(_g, true);
  DegNodeMap deg(_g, 0);
  size_t maxDeg = 0;
  NodeListItMap it(_g);

  // compute node degrees, O(|E|) time
  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    size_t d = 0;
    for (IncEdgeIt e(_g, v); e != lemon::INVALID; ++e, ++d);
    deg[v] = d;
    if (d > maxDeg) maxDeg = d;
  }

  // fill T, O(d) time
  NodeListVector T(maxDeg + 1, NodeList());
  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    size_t d = deg[v];
    T[d].push_front(v);
    it[v] = T[d].begin();
  }

  size_t degeneracy = 0;

  // O(|V|) time, Eppstein et al. (2010)
  const size_t n = T.size();
  size_t i = 0;
  while (i < n)
  {
    NodeList& l = T[i];
    if (T[i].size() > 0)
    {
      Node v = l.front();
      l.pop_front();
      order.push_back(v);
      present[v] = false;
      if (deg[v] > degeneracy)
      {
        degeneracy = deg[v];
      }
      //std::cout << "Removed " << _g.id(v) << std::endl;

      for (IncEdgeIt e(_g, v); e != lemon::INVALID; ++e)
      {
        Node w = _g.oppositeNode(v, e);
        if (present[w])
        {
          size_t deg_w = deg[w];
          typename NodeList::iterator it_w = it[w];

          T[deg_w - 1].splice(T[deg_w - 1].begin(), T[deg_w], it_w);
          deg[w]--;
        }
      }

      i = 0;
    }
    else
    {
      ++i;
    }
  }

  //std::cerr << "Degeneracy: " << degeneracy << std::endl;
  return degeneracy;
}

void BronKerbosch::run(SolverType type)
{

  switch (type)
  {
  case BK_CLASSIC:
  case BK_PIVOT:
    {
      BitSet P(_n), R(_n), X(_n);
      P.set();

      if (type == BK_CLASSIC)
        bkClassic(P, R, X);
      else
        bkPivot(P, R, X);
    }
    break;
  case BK_PIVOT_DEGENERACY:
    {
      NodeList order;
      computeDegeneracy(order);
      bkDegeneracy(order);
    }
    break;
  }
}

void BronKerbosch::report(const BitSet& R)
{
  NodeVector clique;
  for (size_t i = 0; i < R.size(); ++i)
  {
    if (R[i])
    {
      if (g_verbosity >= VERBOSE_DEBUG)
      {
        std::cerr << " " << _g.id(_bitToNode[i]);
      }
      clique.push_back(_bitToNode[i]);
    }
  }
  _cliques.push_back(clique);
  if (g_verbosity >= VERBOSE_DEBUG)
    std::cerr << std::endl;
}

void BronKerbosch::printBitSet(const BitSet& S, std::ostream& out) const
{
  out << "{";
  bool first = true;
  for (size_t i = 0; i < S.size(); ++i)
  {
    if (S[i])
    {
      if (!first)
      {
        out << ", ";
      }
      else
      {
        first = false;
      }
      out << _g.id(_bitToNode[i]);
    }
  }
  out << "}";
}

void BronKerbosch::bkDegeneracy(const NodeList& order)
{
  BitSet mask(_n);

  for (typename NodeList::const_iterator it = order.begin(); it != order.end(); ++it)
  {
    Node v = *it;
    // ~mask includes v but we're fine as _bitNeighborhood[v] excludes v
    BitSet P = _bitNeighborhood[v] & ~mask;
    BitSet X = _bitNeighborhood[v] & mask;
    BitSet R(_n);
    R.set(_nodeToBit[v]);

    bkPivot(P, R, X);
    mask.set(_nodeToBit[v]);
  }
}

void BronKerbosch::bkPivot(BitSet P, BitSet R, BitSet X)
{
  assert((P & X).none());
  assert((P & R).none());
  assert((R & X).none());

  // let's print P, R and X
  //std::cout << "P = ";
  //print(P, std::cout);
  //std::cout << ", R = ";
  //print(R, std::cout);
  //std::cout << ", X = ";
  //print(X, std::cout);
  //std::cout << std::endl;
  
  if (_limit != -1 && _cliques.size() > _limit)
  {
    // stop enumerating when limit is reached
    return;
  }

  // Reports maximal cliques in P \cup R (but not in X)
  BitSet P_cup_X = P | X;
  if (P_cup_X.none())
  {
    report(R);
  }
  else
  {
    // choose a pivot u from (P | X) s.t |P & N(u)| is maximum, Tomita et al. (2006)
    size_t maxBitCount = 0;
    Node max_u = lemon::INVALID;
    for (size_t i = 0; i < P.size(); ++i)
    {
      if (P_cup_X[i])
      {
        Node u = _bitToNode[i];
        BitSet P_cap_Nu = P & _bitNeighborhood[u];
        size_t s = P_cap_Nu.count();
        if (s >= maxBitCount)
        {
          max_u = u;
          maxBitCount = s;
        }
      }
    }

    assert(max_u != lemon::INVALID);
    BitSet P_diff_Nu = P - _bitNeighborhood[max_u];
    for (size_t i = 0; i < P.size(); ++i)
    {
      if (P_diff_Nu[i])
      {
        Node v = _bitToNode[i];
        BitSet R_ = R;
        R_[_nodeToBit[v]] = 1;
        // report all maximal cliques in ( (P | N[v]) & R) \ (X & N[v]) )
        bkPivot(P & _bitNeighborhood[v], R_, X & _bitNeighborhood[v]);
        P[i] = 0;
        X[i] = 1;
      }
    }
  }

  /* Invariants:
   * - R is a clique
   * - Each node v in P is adjacent to all nodes in R
   *   (nodes in P are used to extend R, upon usage it's moved to X)
   * - Each node v in X is adjacent to all nodes in R
   *   (maximal cliques containing R \cup v have already been reported)
   */
}

void BronKerbosch::bkClassic(BitSet P, BitSet R, BitSet X)
{
  assert((P & X).none());
  assert((P & R).none());
  assert((R & X).none());

  // let's print P, R and X
  //std::cout << "P = ";
  //print(P, std::cout);
  //std::cout << ", R = ";
  //print(R, std::cout);
  //std::cout << ", X = ";
  //print(X, std::cout);
  //std::cout << std::endl;

  // Reports maximal cliques in P \cup R (but not in X)
  BitSet P_cup_X = P | X;
  if (P_cup_X.none())
  {
    report(R);
  }
  else
  {
    for (size_t i = 0; i < P.size(); ++i)
    {
      if (P[i])
      {
        Node v = _bitToNode[i];
        BitSet R_ = R;
        R_[_nodeToBit[v]] = 1;
        // report all maximal cliques in ( (P | N[v]) & R) \ (X & N[v]) )
        bkClassic(P & _bitNeighborhood[v], R_, X & _bitNeighborhood[v]);
        P[i] = 0;
        X[i] = 1;
      }
    }
  }

  /* Invariants:
   * - R is a clique
   * - Each node v in P is adjacent to all nodes in R
   *   (nodes in P are used to extend R, upon usage it's moved to X)
   * - Each node v in X is adjacent to all nodes in R
   *   (maximal cliques containing R \cup v have already been reported)
   */
}

void BronKerbosch::print(std::ostream& out) const
{
  for (size_t i = 0; i < _cliques.size(); ++i)
  {
    bool first = true;
    const NodeVector& clique = _cliques[i];
    out << clique.size() << ": ";
    for (size_t j = 0; j < clique.size(); ++j)
    {
      if (!first)
      {
        out << ", ";
      }
      else
      {
        first = false;
      }
      out << _g.id(clique[j]);
    }
    out << std::endl;
  }
}
  
} // namespace gm
