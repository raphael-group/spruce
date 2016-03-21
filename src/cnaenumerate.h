/*
 * cnaenumerate.h
 *
 *  Created on: 19-oct-2015
 *      Author: M. El-Kebir
 */

#ifndef CNAENUMERATE_H
#define CNAENUMERATE_H

#include "charactermatrix.h"
#include "solutionset.h"

namespace gm {
  
class CnaEnumerate
{
public:
  typedef std::vector<StateTree> StateTreeVector;
  
  CnaEnumerate(const CharacterMatrix& M);
  
  const CharacterMatrix& M() const
  {
    return _M;
  }
  
  void enumerate(int limit,
                 int timeLimit,
                 int threads,
                 bool monoclonal,
                 const IntSet& whiteList);
  
  const SolutionSet& sols() const
  {
    return _sols;
  }
  
  void init();
  
  int combinations() const
  {
    return _combinations;
  }
  
  int real_n() const
  {
    return _real_n;
  }
  
  void get(int combination,
           RealTensor& F,
           StateTreeVector& S) const;
  
private:
  bool next(StlIntVector& pi) const;
  
  void solve(const StlIntVector& pi,
             int limit,
             int timeLimit,
             int threads,
             bool monoclonal,
             const IntSet& whiteList);
  
  void init(const StlIntVector& pi,
            RealTensor& F,
            StateTreeVector& S) const;
  
private:
  const CharacterMatrix& _M;
  SolutionSet _sols;
  int _treeSize;

  int _combinations;
  int _real_n;
};
  
} // namespace gm

#endif // CNAENUMERATE_H
