/*
 * character.h
 *
 *  Created on: 13-jan-2016
 *      Author: M. El-Kebir
 */

#ifndef CHARACTER_H
#define CHARACTER_H

#include <list>
#include "utils.h"
#include "stategraph.h"

namespace gm {
  
class Character
{
public:
  Character();
  
  struct CopyState
  {
  public:
    CopyState(int x, int y, double mu)
      : _x(x)
      , _y(y)
      , _mu(mu)
    {
    }
    
    CopyState()
      : _x(-1)
      , _y(-1)
      , _mu(0)
    {
    }
    
    int x() const { return _x; }
    int y() const { return _y; }
    double mu() const { return _mu; }
    
  private:
    int _x;
    int _y;
    double _mu;
  };
  
  typedef std::list<CopyState> CopyStateList;
  typedef CopyStateList::const_iterator CopyStateListIt;
  typedef CopyStateList::iterator CopyStateListNonConstIt;
  typedef std::map<StateGraph::StateEdgeSet, StlRealIntervalTensor> FrequencyMap;
  typedef FrequencyMap::const_iterator FrequencyMapIt;
  typedef FrequencyMap::iterator FrequencyMapNonConstIt;
  typedef StateGraph::IntPairSet IntPairSet;
  
  double vafLB() const { return _vafLB; }
  double vafUB() const { return _vafUB; }
  const CopyStateList& L() const { return _L; }
  
  IntPairSet L0() const
  {
    IntPairSet res;
    for (CopyStateListIt it = _L.begin(); it != _L.end(); ++it)
    {
      if (it->mu() == 0)
        res.insert(std::make_pair(it->x(), it->y()));
    }
    return res;
  }
  
  double mu(int x, int y) const
  {
    assert(0 <= x && x < _M.size());
    assert(0 <= y && y < _M[x].size());
    
    return _M[x][y];
  }
  
  void remove(const IntPairSet& L0);
  
  int sampleIndex() const { return _sampleIndex; }
  int characterIndex() const { return _characterIndex; }
  
  const std::string& sampleLabel() const { return _sampleLabel; }
  const std::string& characterLabel() const { return _characterLabel; }
  
  void solve(int xy_max_c,
             int xy_max_overall,
             bool includeMutationEdge,
             FrequencyMap& mapToF);
  
private:
  friend std::ostream& operator<<(std::ostream& out, const Character& c);
  friend std::istream& operator>>(std::istream& in, Character& c);
  
  static RealInterval intersect(double lb1, double ub1,
                                double lb2, double ub2)
  {
    RealInterval res;
    res.first = std::max(lb1, lb2);
    res.second = std::min(ub1, ub2);
    return res;
  }
  
private:
  double _vafLB;
  double _vaf;
  double _vafUB;
  CopyStateList _L;
  StlDoubleMatrix _M;

  int _sampleIndex;
  int _characterIndex;
  
  std::string _sampleLabel;
  std::string _characterLabel;
};
  
std::ostream& operator<<(std::ostream& out, const Character& c);
std::istream& operator>>(std::istream& in, Character& c);

} // namespace gm

#endif // CHARACTER_H
