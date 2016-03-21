/*
 * charactermatrix.h
 *
 *  Created on: 15-jan-2016
 *      Author: M. El-Kebir
 */

#ifndef CHARACTERMATRIX_H
#define CHARACTERMATRIX_H

#include "utils.h"
#include "character.h"
#include "statetree.h"

namespace gm {

class CharacterMatrix
{
public:
  typedef StateGraph::CnaTriple CnaTriple;
  typedef std::vector<IntSet> IntSetVector;
  typedef std::set<IntSet> IntSetSet;
  typedef IntSetSet::const_iterator IntSetSetIt;
  typedef IntSetSet::iterator IntSetSetNonConstIt;
  typedef std::vector<IntSetSet> IntSetSetVector;

  CharacterMatrix();
  
  void init();
  
  void writeCompatibleStateTrees(const StlDoubleVector& purityValues,
                                 std::ostream& out);
  
  void writeCompatibleStateTrees(std::ostream& out, int c);
  
  void setIntervals(std::istream& in);
  
  const IntSet& interval(int c) const
  {
    assert(0 <= c && c < n());
    return _charToInterval[c];
  }
  
  const IntSetSet& intervals() const
  {
    return _intervals;
  }
  
  bool isFeasible(int c) const
  {
    assert(0 <= 0 && c < _n);
    return !_F[0][c].empty();
  }
  
  int maxX(int c) const
  {
    assert(0 <= c && c < _n);
    return _maxX[c];
  }
  
  int m() const
  {
    return _m;
  }
  
  int n() const
  {
    return _n;
  }
  
  int k() const
  {
    return _stateToTriple.size();
  }
  
  const std::string& label(int i) const
  {
    assert(0 <= i && i < k());
    return _stateLabel[i];
  }
  
  int tripleToState(const CnaTriple& triple) const
  {
    CnaTripleMapIt it = _tripleToState.find(triple);
    if (it == _tripleToState.end())
      return -1;
    else
      return it->second;
  }
  
  const CnaTriple& stateToTriple(int i) const
  {
    assert(0 <= i && i < k());
    return _stateToTriple[i];
  }
  
  int numStateTrees(int c) const
  {
    assert(0 <= c && c < _n);
    
    return _F[0][c].size();
  }
  
  const Character& operator()(int p, int c) const
  {
    assert(0 <= p && p < _m);
    assert(0 <= c && c < _n);
    
    return _M[p][c];
  }
  
  StateTree stateTree(int s, int c) const;
  
  StlRealIntervalVector ccf(const StlDoubleVector& purityValues, int s, int c) const;
  
  typedef std::pair<IntPair, IntPair> IntPairPair;
  typedef std::set<IntPairPair> IntPairPairSet;
  
  IntPairPairSet copyTree(int s, int c) const;
  
  void applyHeuristic();
  
private:
  typedef std::vector<Character> StlCharacterVector;
  typedef std::vector<StlCharacterVector> StlCharacterMatrix;
  
  typedef std::vector<Character::FrequencyMap> FrequencyMapVector;
  typedef std::vector<FrequencyMapVector> FrequencyMapMatrix;
  
  typedef std::map<StateGraph::CnaTriple, int> CnaTripleMap;
  typedef CnaTripleMap::const_iterator CnaTripleMapIt;
  typedef StateGraph::CnaTripleVector CnaTripleVector;
  
  typedef std::vector<std::string> StringVector;

public:
  double get_f_lb(int s, int p, int c, int i) const
  {
    assert(0 <= i && i < k());
    return get_f_lb(s, p, c, _stateToTriple[i]._x, _stateToTriple[i]._y, _stateToTriple[i]._z);
  }
  
  double get_f_lb(int s, int p, int c, int x, int y, int z) const
  {
    assert(0 <= p && p < _m);
    assert(0 <= c && c < _n);
    assert(0 <= s && s < _F[p][c].size());
    assert(0 <= x && x <= _maxXY);
    assert(0 <= y && y <= _maxXY);
    assert(0 <= z && z <= x);
    
    Character::FrequencyMapIt it = _F[p][c].begin();
    for (int i = 0; i < s; ++i) ++it;
    
    const StlRealIntervalTensor& freqMap = it->second;
    
    return freqMap[x][y][z].first;
  }
  
  double get_f_ub(int s, int p, int c, int i) const
  {
    assert(0 <= i && i < k());
    return get_f_ub(s, p, c, _stateToTriple[i]._x, _stateToTriple[i]._y, _stateToTriple[i]._z);
  }
  
  double get_f_ub(int s, int p, int c, int x, int y, int z) const
  {
    assert(0 <= p && p < _m);
    assert(0 <= c && c < _n);
    assert(0 <= s && s < _F[p][c].size());
    assert(0 <= x && x <= _maxXY);
    assert(0 <= y && y <= _maxXY);
    assert(0 <= z && z <= x);
    
    Character::FrequencyMapIt it = _F[p][c].begin();
    for (int i = 0; i < s; ++i) ++it;
    
    return it->second[x][y][z].second;
  }
  
private:
  int _m;
  int _n;
  StlCharacterMatrix _M;
  FrequencyMapMatrix _F;
  StlIntVector _maxX;
  int _maxXY;
  
  CnaTripleMap _tripleToState;
  CnaTripleVector _stateToTriple;
  StringVector _stateLabel;
  
  IntSetVector _charToInterval;
  IntSetSet _intervals;
  
  friend std::ostream& operator<<(std::ostream& out, const CharacterMatrix& M);
  friend std::istream& operator>>(std::istream& in, CharacterMatrix& M);
};

std::ostream& operator<<(std::ostream& out, const CharacterMatrix& M);
std::istream& operator>>(std::istream& in, CharacterMatrix& M);

} // namespace gm

#endif // CHARACTERMATRIX_H
