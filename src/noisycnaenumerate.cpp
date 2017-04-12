/*
 * noisycnaenumerate.cpp
 *
 *  Created on: 19-oct-2015
 *      Author: M. El-Kebir
 */

#include "noisycnaenumerate.h"
#include "realtensor.h"
#include "rootedcladisticnoisyenumeration.h"
#include <stdio.h>

namespace gm {

NoisyCnaEnumerate::NoisyCnaEnumerate(const CharacterMatrix& M,
                                     const StlDoubleVector& purityValues,
                                     const CompatibilityGraph& comp,
                                     int lowerbound)
  : _M(M)
  , _purityValues(purityValues)
  , _comp(comp)
  , _sols()
  , _treeSize(lowerbound)
  , _combinations(-1)
  , _real_n(-1)
{
}
  
void NoisyCnaEnumerate::init(int state_tree_limit)
{
  const int m = _M.m();
  const int n = _M.n();
  assert(m > 0 && n > 0);
  
  _combinations = _comp.combinations();
}
  
void NoisyCnaEnumerate::enumerate(int limit,
                                  int timeLimit,
                                  int threads,
                                  int state_tree_limit,
                                  bool monoclonal,
                                  int offset,
                                  const IntSet& whiteList)
{
  const int n = _M.n();
  
  _sols.clear();
  
  StlIntVector pi(n, 0);
  pi[0] = offset;
  int count = 0;
  do {
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      if (state_tree_limit != -1)
      {
        std::cerr << std::endl << "State tree combination " << ++count << "/" << state_tree_limit << " ..." << std::endl;
      }
      else
      {
        std::cerr << std::endl << "State tree combination " << ++count << "/" << _combinations << " ..." << std::endl;
      }
    }
    solve(pi, limit, timeLimit, threads, state_tree_limit, monoclonal, whiteList);
  } while (next(state_tree_limit, offset, pi));
}
  
bool NoisyCnaEnumerate::next(int state_tree_limit,
                             int offset,
                             StlIntVector& pi) const
{
  if (pi[0] < _comp.combinations() - 1 &&
      (pi[0] - offset < state_tree_limit - 1 || state_tree_limit == -1))
  {
    ++pi[0];
    return true;
  }
  else
  {
    return false;
  }
}

void NoisyCnaEnumerate::collapse(const StlIntVector& mapNewCharToOldChar,
                                 const StlIntVector& mapOldCharToNewChar,
                                 RootedCladisticNoisyAncestryGraph& G)
{
  typedef std::set<IntPair> IntPairSet;
  
  int k = _M.k();
  const auto& intervals = _M.intervals();
  for (const IntSet& interval : intervals)
  {
    IntSet remappedInterval;
    for (int c : interval)
    {
      int cc = mapOldCharToNewChar[c];
      if (cc != -1)
        remappedInterval.insert(cc);
    }
    if (remappedInterval.size() > 1)
    {
      // get the copy states
      IntPairSet XY;
      const StateTree& S = G.S(*remappedInterval.begin());
      for (int i = 0; i < k; ++i)
      {
        if (S.isPresent(i))
        {
          const auto& xyz = _M.stateToTriple(i);
          // skip state 1,1
          if (xyz._x != 1 || xyz._y != 1)
          {
            XY.insert(IntPair(xyz._x, xyz._y));
          }
        }
      }
      
      for (const IntPair& xy : XY)
      {
        assert(xy.first != 1 || xy.second != 1);
        
        // collect all char-state pairs correspond to CNAs
        IntPairSet toCollapse;
        for (int c : remappedInterval)
        {
          const StateTree& S_c = G.S(c);
          for (int i = 0; i < k; ++i)
          {
            if (S_c.isPresent(i) && _M.stateToTriple(i)._x == xy.first && _M.stateToTriple(i)._y == xy.second)
            {
              int pi_i = S_c.parent(i);
              assert(0 <= pi_i && pi_i < k);
              
              if (_M.stateToTriple(pi_i)._x != xy.first || _M.stateToTriple(pi_i)._y != xy.second)
              {
                // we got a CNA state
                toCollapse.insert(IntPair(c, i));
              }
            }
          }
        }
        
        G.collapse(toCollapse);
      }
    }
  }
}
  
void NoisyCnaEnumerate::fixTrunk(const RealTensor& F_ub,
                                 const StateTreeVector& S,
                                 const StlIntVector& mapNewCharToOldChar,
                                 const StlIntVector& mapOldCharToNewChar,
                                 RootedCladisticNoisyAncestryGraph& G)
{
  const int m = F_ub.m();
  const int k = F_ub.k();
  const int n = F_ub.n();
  
  typedef std::set<IntPair> IntPairSet;
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Truncal characters (CCF_lb >= 1, across all samples):" << std::endl;
  }
  
  IntPairSet truncalCharacters;
  for (int c = 0; c < n; ++c)
  {
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cerr << "Character " << _M(0, mapNewCharToOldChar[c]).characterLabel() << " (" << mapNewCharToOldChar[c] << ") : ";
      S[c].writeEdgeList(std::cerr);
      std::cerr << std::endl;
    }
    
    bool truncal = true;
    int mutation_i = -1;
    for (int p = 0; p < m; ++p)
    {
      double ccf = 0;
      for (int i = 0; i < k; ++i)
      {
        if (!S[c].isPresent(i))
          continue;
        
        const auto& xyz = _M.stateToTriple(i);
        if (xyz._z >= 1)
        {
          ccf += F_ub(i, p, c);
          
          const auto& pi_xyz = _M.stateToTriple(S[c].parent(i));
          if (pi_xyz._z == 0 && xyz._z == 1)
          {
            mutation_i = i;
          }
        }
      }
      ccf /= _purityValues[p];
      
      if (g_tol.less(ccf, 1))
      {
        truncal = false;
        break;
      }
    }
    
    if (truncal)
    {
      assert(mutation_i != -1);
      const StateTree& S_c = S[c];
      truncalCharacters.insert(IntPair(c, mutation_i));
      while ((mutation_i = S_c.parent(mutation_i)) != 0)
      {
        truncalCharacters.insert(IntPair(c, mutation_i));
      }
    }
  }
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Truncal character-state pairs:" << std::endl;
    for (auto ci : truncalCharacters)
    {
      auto xyz = _M.stateToTriple(ci.second);
      std::cerr << _M(0, mapNewCharToOldChar[ci.first]).characterLabel()
                << " , (" << xyz._x << "," << xyz._y << "," << xyz._z << ")" << std::endl;
    }
  }
  
  Digraph::Node monoClonalRoot = G.collapse(truncalCharacters);
  G.makeMonoclonal(monoClonalRoot);
}
  
void NoisyCnaEnumerate::solve(const StlIntVector& pi,
                              int limit,
                              int timeLimit,
                              int threads,
                              int state_tree_limit,
                              bool monoclonal,
                              const IntSet& whiteList)
{
  RealTensor F, F_lb, F_ub;
  StateTreeVector S;
  StlIntVector mapNewCharToOldChar;
  StlIntVector mapOldCharToNewChar;
  
  _comp.init(pi[0], F, S, F_lb, F_ub, mapNewCharToOldChar, mapOldCharToNewChar);
  
  RootedCladisticNoisyAncestryGraph G(F, S, F_lb, F_ub);
  G.init();
  
  collapse(mapNewCharToOldChar, mapOldCharToNewChar, G);
  if (monoclonal && !_purityValues.empty())
  {
    fixTrunk(F_ub, S, mapNewCharToOldChar, mapOldCharToNewChar, G);
  }
    
  G.setLabels(F);
//  G.RootedCladisticAncestryGraph::writeDOT(std::cerr);
  
  IntSet remappedWhiteList;
  for (const int c : whiteList)
  {
    int cc = mapOldCharToNewChar[c];
    if (cc != -1)
    {
      remappedWhiteList.insert(cc);
    }
  }
  
  RootedCladisticNoisyEnumeration enumerate(G,
                                            limit,
                                            timeLimit,
                                            threads,
                                            _treeSize,
                                            monoclonal,
                                            monoclonal && !_purityValues.empty(),
                                            remappedWhiteList);
  enumerate.run();
  
  if (enumerate.objectiveValue() >= _treeSize)
  {
    if (enumerate.objectiveValue() > _treeSize)
    {
      _sols.clear();
      _treeSize = enumerate.objectiveValue();
    }
    enumerate.populateSolutionSet(_sols);
  }
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
    {
      std::cerr << std::endl;
    }
    std::cerr << "Tree size: " << enumerate.objectiveValue() << "/" << _treeSize
              << " (" << enumerate.solutionCount() << " trees)" << std::endl;
  }
}
  
void NoisyCnaEnumerate::get(int state_tree_limit,
                            int combination,
                            RealTensor& F,
                            StateTreeVector& S,
                            RealTensor& F_lb,
                            RealTensor& F_ub) const
{
  assert(0 <= combination && combination < _combinations);

  StlIntVector mapNewCharToOldChar, mapOldCharToNewChar;
  const int n = _M.n();
  
  StlIntVector pi(n, 0);
  for (int count = 0; count < combination; ++count)
  {
    next(state_tree_limit, 0, pi);
  }
  
  _comp.init(pi[0], F, S, F_lb, F_ub, mapNewCharToOldChar, mapOldCharToNewChar);
}
  
} // namespace gm
