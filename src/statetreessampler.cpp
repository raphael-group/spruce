/*
 * statetreessampler.cpp
 *
 *  Created on: 6-feb-2016
 *      Author: M. El-Kebir
 */

#include "statetreessampler.h"
#include "utils.h"
#include <algorithm>

namespace gm {

StateTreesSampler::StateTreesSampler(int seed,
                                     int maxXY,
                                     int n,
                                     int maxCopyEvents)
  : _maxXY(maxXY)
  , _n(n)
  , _maxCopyEvents(maxCopyEvents)
  , _generator(seed)
  , _sampleSpace()
{
  // initialize sample space
  initSampleSpace();
}

IntPair StateTreesSampler::sampleCopyState()
{
  std::uniform_int_distribution<> dis(0, _sampleSpace.size() - 1);
  
  int i = dis(_generator);
  return _sampleSpace[i];
}
  
StateTreesSampler::StateTreeVector StateTreesSampler::sample()
{
  CnaTripleMap tripleToState;
  CnaTripleVector stateToTriple;
  StringVector stateLabel;
  
  StateTreeVector S;
  
  // now sample state trees
  StateEdgeSetVector stateTreePerCharacter(_n);
  
  std::uniform_int_distribution<> dis_02(0, _maxCopyEvents);
  for (int c = 0; c < _n; ++c)
  {
    int numberOfNonNormalCopyStates = dis_02(_generator);
    
    StateGraph::IntPairSet L;
    L.insert(std::make_pair(1, 1));
    for (int i = 0; i < numberOfNonNormalCopyStates; ++i)
    {
      IntPair xy = sampleCopyState();
      L.insert(xy);
//      std::cout << xy.first << "," << xy.second << std::endl;
    }
    
    StateGraph::StateEdgeSetSet setS = StateGraph::getStateTrees(L, _maxXY, true);

    // remove sets that don't have a mutation state
    for (StateGraph::StateEdgeSetSetNonConstIt it2 = setS.begin(); it2 != setS.end();)
    {
      bool remove = true;
      
      for (StateGraph::StateEdgeSetIt it3 = it2->begin(); it3 != it2->end(); ++it3)
      {
        if (it3->second._z > 0)
        {
          remove = false;
          break;
        }
      }
      
      if (remove)
        it2 = setS.erase(it2);
      else
        ++it2;
    }
    
    std::uniform_int_distribution<> dis_0s(0, setS.size() - 1);
    int n = dis_0s(_generator);
    
    StateGraph::StateEdgeSetSetIt it = setS.begin();
    for (int i = 0; i < n; ++i)
    {
      ++it;
    }
    
    stateTreePerCharacter[c] = *it;
  }
  
  // now construct state trees
  StateGraph::CnaTriple triple110(1,1,0);
  stateToTriple.push_back(triple110);
  tripleToState[triple110] = 0;
  
  for (int c = 0; c < _n; ++c)
  {
    const StateGraph::StateEdgeSet& edgeSet = stateTreePerCharacter[c];
    for (StateGraph::StateEdgeSetIt it2 = edgeSet.begin(); it2 != edgeSet.end(); ++it2)
    {
      if (tripleToState.find(it2->first) == tripleToState.end())
      {
        tripleToState[it2->first] = stateToTriple.size();
        stateToTriple.push_back(it2->first);
      }
      if (tripleToState.find(it2->second) == tripleToState.end())
      {
        tripleToState[it2->second] = stateToTriple.size();
        stateToTriple.push_back(it2->second);
      }
    }
  }
  
  const int kk = stateToTriple.size();
  stateLabel = StringVector(kk);
  char buf[1024];
  for (int i = 0; i < kk; ++i)
  {
    const StateGraph::CnaTriple& triple = stateToTriple[i];
    snprintf(buf, 1024, "(%d,%d,%d)", triple._x, triple._y, triple._z);
    stateLabel[i] = buf;
  }
  
  // construct state tree vector
  S.clear();
  for (int c = 0; c < _n; ++c)
  {
    StlIntVector pi(stateToTriple.size(), -2);
    pi[0] = -1;
    
    const StateGraph::StateEdgeSet& edgeSet = stateTreePerCharacter[c];
    for (StateGraph::StateEdgeSetIt edgeIt = edgeSet.begin(); edgeIt != edgeSet.end(); ++edgeIt)
    {
      assert(tripleToState.find(edgeIt->first) != tripleToState.end());
      assert(tripleToState.find(edgeIt->second) != tripleToState.end());
      pi[tripleToState.find(edgeIt->second)->second] = tripleToState.find(edgeIt->first)->second;
    }
    
    StateTree SS(pi);
    for (int i = 0; i < kk; ++i)
    {
      SS.setLabel(i, stateLabel[i]);
    }
    
    S.push_back(SS);
  }
  
  return S;
}
  
void StateTreesSampler::initSampleSpace()
{
  for (int x = 0; x <= _maxXY; ++x)
  {
    for (int y = 0; y <= x; ++y)
    {
      // skip state (0,0)
      if (x == 0 && y == 0) continue;
      
      IntPair xy = std::make_pair(x, y);
      
      int d = 2 * _maxXY - (abs(x - 1) + abs(y - 1));
//      std::cout << xy.first << "," << xy.second << " " << d << std::endl;
      for (int i = 0; i < d; ++i)
      {
        _sampleSpace.push_back(xy);
      }
    }
  }
}

} // namespace gm
