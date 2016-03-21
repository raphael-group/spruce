/*
 * rootedcladisticsampler.h
 *
 *  Created on: 19-oct-2015
 *      Author: M. El-Kebir
 */

#ifndef ROOTEDCLADISTICSAMPLER_H
#define ROOTEDCLADISTICSAMPLER_H

#include <lemon/adaptors.h>
#include <random>
#include "utils.h"
#include "rootedcladisticsampler.h"
#include "rootedcladisticfullancestrygraph.h"

namespace gm {

class RootedCladisticSampler
{
public:
  DIGRAPH_TYPEDEFS(Digraph);
  
  RootedCladisticSampler(const RootedCladisticFullAncestryGraph& G,
                         int seed);
  
  void sample();
  
  void writeEdgeList(std::ostream& out) const;
  
protected:
  typedef std::list<Arc> ArcList;
  typedef ArcList::const_iterator ArcListIt;
  typedef ArcList::iterator ArcListNonConstIt;
  typedef ArcList::const_reverse_iterator ArcListRevIt;
  typedef std::vector<ArcList> ArcListVector;
  typedef std::list<ArcList> ArcListList;
  typedef lemon::SubDigraph<const Digraph> SubDigraph;
  typedef SubDigraph::ArcIt SubArcIt;
  typedef SubDigraph::NodeIt SubNodeIt;
  typedef SubDigraph::OutArcIt SubOutArcIt;
  typedef SubDigraph::InArcIt SubInArcIt;
  typedef Digraph::NodeMap<IntSet> IntSetNodeMap;
  
  void init(SubDigraph& subG, SubDigraph& T, ArcList& F);
  
  bool isFirstAncestor(const SubDigraph& T,
                       int c,
                       Node v_ci,
                       Node v_dj) const;
  
  bool isAncestor(const SubDigraph& T,
                  Node v_dj,
                  Node v_ci) const;
  
  bool isArborescence(const SubDigraph& T) const;
  
  void addArc(SubDigraph& T,
              Arc a_cidj) const;
  void removeArc(SubDigraph& T,
                 Arc a_cidj) const;
  
private:
  void grow(SubDigraph& G,
            SubDigraph& T,
            ArcList& F);
  
  virtual bool isValid(const SubDigraph& T) const;
  bool isValid(const SubDigraph& T, Arc a_ciel) const;
  
  void writeDOT(std::ostream& out,
                const SubDigraph& T) const;
  
  void generateResult(const SubDigraph& T, Node v_ci);
  
  int rand(int a, int b)
  {
    std::uniform_int_distribution<> dis(a, b);
    return dis(_generator);
  }
  
protected:
  const RootedCladisticFullAncestryGraph& _G;
  ArcList _result;
  std::mt19937 _generator;
};
  
} // namespace gm

#endif // ROOTEDCLADISTICSAMPLER