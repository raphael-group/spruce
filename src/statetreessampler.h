/*
 * statetreessampler.h
 *
 *  Created on: 6-feb-2016
 *      Author: M. El-Kebir
 */

#ifndef STATETREESSAMPLER_H
#define STATETREESSAMPLER_H

#include "utils.h"
#include "stategraph.h"
#include "statetree.h"
#include <random>

namespace gm {

class StateTreesSampler
{
public:
  StateTreesSampler(int seed,
                    int maxXY,
                    int n,
                    int maxCopyEvents);
  
  typedef std::vector<StateTree> StateTreeVector;
  typedef std::vector<IntPair> IntPairVector;
  typedef std::map<StateGraph::CnaTriple, int> CnaTripleMap;
  typedef CnaTripleMap::const_iterator CnaTripleMapIt;
  typedef StateGraph::CnaTripleVector CnaTripleVector;
  typedef std::vector<std::string> StringVector;
  typedef std::vector<StateGraph::StateEdgeSet> StateEdgeSetVector;

  StateTreeVector sample();

private:
  void initSampleSpace();
  IntPair sampleCopyState();

private:
  const int _maxXY;
  const int _n;
  const int _maxCopyEvents;
  
  std::mt19937 _generator;
  IntPairVector _sampleSpace;
};

} // namespace gm

#endif // STATETREESSAMPLER_H
