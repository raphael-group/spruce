/*
 * stategraph.h
 *
 *  Created on: 11-jan-2016
 *      Author: M. El-Kebir
 */

#ifndef STATEGRAPH_H
#define STATEGRAPH_H

#include <lemon/adaptors.h>
#include <vector>
#include "utils.h"

namespace gm {

class StateGraph
{
public:
  DIGRAPH_TYPEDEFS(Digraph);
  
  StateGraph(int max_x);
  
  typedef enum { MUTATION, AMPLIFICATION, DELETION } EdgeType;
  typedef std::set<IntPair> IntPairSet;
  typedef IntPairSet::const_iterator IntPairSetIt;
  typedef IntPairSet::iterator IntPairSetNonConstIt;
  
  typedef Digraph::ArcMap<EdgeType> EdgeTypeArcMap;
  
  void writeDOT(std::ostream& out) const;
  
  int enumerate(const IntPairSet& L,
                bool includeMutationEdge);
  
private:
  typedef std::vector<Node> NodeVector;
  typedef std::vector<NodeVector> NodeMatrix;
  typedef std::vector<NodeMatrix> Node3Matrix;
  typedef std::vector<Node3Matrix> Node4Matrix;
  
  typedef std::list<Arc> ArcList;
  typedef ArcList::const_iterator ArcListIt;
  typedef ArcList::iterator ArcListNonConstIt;
  typedef ArcList::const_reverse_iterator ArcListRevIt;

  typedef lemon::SubDigraph<const Digraph> SubDigraph;
  typedef SubDigraph::ArcIt SubArcIt;
  typedef SubDigraph::NodeIt SubNodeIt;
  typedef SubDigraph::OutArcIt SubOutArcIt;
  typedef SubDigraph::InArcIt SubInArcIt;
  
  void init();
  
  void init(bool includeMutationEdge,
            SubDigraph& subG,
            SubDigraph& T,
            ArcList& F);
  
  bool finalize(const IntPairSet& L,
                const StlBoolMatrix& LL,
                const SubDigraph& T,
                const StlIntMatrix& V);
  
  void writeDOT(const SubDigraph& T,
                std::ostream& out) const;
  
  void writeDOT(const SubDigraph& T,
                const ArcList& F,
                std::ostream& out) const;
  
  void addEdge(Node v, EdgeType type, int x, int y, int xbar, int ybar)
  {
    Node w = _toNode[x][y][xbar][ybar];
    Arc vw = _G.addArc(v, w);
    _type[vw] = type;
  }
  
  bool isValid(const IntPairSet& L,
               const StlBoolMatrix& LL,
               const SubDigraph& T,
               const StlIntMatrix& V) const;
  
  bool isValid(const SubDigraph& T) const;
  
  bool isValid(const SubDigraph& T,
               Arc a);
  
  void grow(const IntPairSet& L,
            const StlBoolMatrix& LL,
            bool includeMutationEdge,
            SubDigraph& G,
            SubDigraph& T,
            ArcList& F,
            StlIntMatrix& V,
            Arc& mutationEdge);
  
public:
  struct CnaTriple
  {
    CnaTriple()
      : _x(-1)
      , _y(-1)
      , _z(-1)
    {
    }
    
    CnaTriple(int x, int y, int z)
      : _x(x)
      , _y(y)
      , _z(z)
    {
    }
    
    bool operator==(const CnaTriple& other) const
    {
      return _x == other._x && _y == other._y && _z == other._z;
    }
    
    int _x;
    int _y;
    int _z;
  };

  typedef std::set<CnaTriple> CnaTripleSet;
  typedef CnaTripleSet::const_iterator CnaTripleSetIt;
  typedef std::vector<CnaTriple> CnaTripleVector;
  
  typedef std::pair<CnaTriple, CnaTriple> StateEdge;
  typedef std::set<StateEdge> StateEdgeSet;
  typedef StateEdgeSet::const_iterator StateEdgeSetIt;
  
  typedef std::set<StateEdgeSet> StateEdgeSetSet;
  typedef StateEdgeSetSet::const_iterator StateEdgeSetSetIt;
  typedef StateEdgeSetSet::iterator StateEdgeSetSetNonConstIt;

private:
  typedef Digraph::NodeMap<CnaTriple> CnaTripleNodeMap;
  
private:
  const int _max_x;
  
  Digraph _G;
  IntNodeMap _x;
  IntNodeMap _y;
  IntNodeMap _xbar;
  IntNodeMap _ybar;
  EdgeTypeArcMap _type;
  
  Node4Matrix _toNode;
  Node _root;
  
  StateEdgeSetSet _result;
  
public:
  static const StateEdgeSetSet& getStateTrees(const IntPairSet& L,
                                              int xy_max_c,
                                              bool includeMutationEdge);

  static void print(const StateEdgeSet& S,
                    std::ostream& out)
  {
    for (StateEdgeSetIt it = S.begin(); it != S.end(); ++it)
    {
      out << "( (" << it->first._x << "," << it->first._y << "," << it->first._z << ") , ("
          << it->second._x << "," << it->second._y << "," << it->second._z << ") )  ";
    }
  }
  
  static void print(const StateEdgeSetSet& setS,
                    std::ostream& out)
  {
    for (StateEdgeSetSetIt it1 = setS.begin(); it1 != setS.end(); ++it1)
    {
      const StateEdgeSet& S = *it1;
      print(S, std::cout);
    }
  }
  
private:
  typedef std::map<IntPairSet, StateEdgeSetSet> Dictionary;
  static Dictionary _dict;
//  static std::mapStateGraph* _pG;
};
  
bool operator<(const StateGraph::CnaTriple& lhs, const StateGraph::CnaTriple& rhs);
  
} // namespace gm

#endif // STATEGRAPH_H
