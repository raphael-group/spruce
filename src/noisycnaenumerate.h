/*
 * noisycnaenumerate.h
 *
 *  Created on: 19-oct-2015
 *      Author: M. El-Kebir
 */

#ifndef NOISYCNAENUMERATE_H
#define NOISYCNAENUMERATE_H

#include "charactermatrix.h"
#include "solutionset.h"
#include "compatibilitygraph.h"
#include "rootedcladisticnoisyancestrygraph.h"

namespace gm {
  
class NoisyCnaEnumerate
{
public:
  typedef std::vector<StateTree> StateTreeVector;
  
  NoisyCnaEnumerate(const CharacterMatrix& M,
                    const StlDoubleVector& purityValues,
                    const CompatibilityGraph& comp,
                    int lowerbound);
  
  const CharacterMatrix& M() const
  {
    return _M;
  }
  
  void enumerate(int limit,
                 int timeLimit,
                 int threads,
                 int state_tree_limit,
                 bool monoclonal,
                 int offset,
                 const IntSet& whiteList);
  
  const SolutionSet& sols() const
  {
    return _sols;
  }
  
  void init(int state_tree_limit);
  
  int combinations() const
  {
    return _combinations;
  }
  
  int real_n() const
  {
    return _real_n;
  }
  
  void get(int state_tree_limit,
           int combination,
           RealTensor& F,
           StateTreeVector& S,
           RealTensor& F_lb,
           RealTensor& F_ub) const;
  
private:
  bool next(int state_tree_limit,
            int offset,
            StlIntVector& pi) const;
  
  void solve(const StlIntVector& pi,
             int limit,
             int timeLimit,
             int threads,
             int state_tree_limit,
             bool monoclonal,
             const IntSet& whiteList);
    
  void collapse(const StlIntVector& mapNewCharToOldChar,
                const StlIntVector& mapOldCharToNewChar,
                RootedCladisticNoisyAncestryGraph& G);
  
  void fixTrunk(const RealTensor& F_ub,
                const StateTreeVector& S,
                const StlIntVector& mapNewCharToOldChar,
                const StlIntVector& mapOldCharToNewChar,
                RootedCladisticNoisyAncestryGraph& G);
  
private:
  const CharacterMatrix& _M;
  const StlDoubleVector& _purityValues;
  const CompatibilityGraph& _comp;
  
  SolutionSet _sols;
  int _treeSize;
  
  unsigned long _combinations;
  int _real_n;
};
  
} // namespace gm

#endif // CNAENUMERATE_H
