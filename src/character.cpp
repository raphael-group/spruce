/*
 * character.cpp
 *
 *  Created on: 13-jan-2016
 *      Author: M. El-Kebir
 */

#include "character.h"
#include "utils.h"
#include "stategraph.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace gm {
  
Character::Character()
  : _vafLB(0)
  , _vaf(0)
  , _vafUB(1)
  , _L()
  , _M()
  , _sampleIndex(-1)
  , _characterIndex(-1)
  , _sampleLabel()
  , _characterLabel()
{
}
  
void Character::remove(const IntPairSet& L0)
{
  for (StateGraph::IntPairSetIt it = L0.begin(); it != L0.end(); ++it)
  {
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cerr << "Removing copy-number state (" << it->first << "," << it->second
                << ") from character " << _characterLabel
                << " (" << _characterIndex
                << ") in sample " << _sampleLabel
                << " (" << _sampleIndex << "), mu = 0 across all samples"
                << std::endl;
    }
    
    assert(_M[it->first][it->second] == 0);
    for (CopyStateListNonConstIt it2 = _L.begin(); it2 != _L.end();)
    {
      if (it2->x() == it->first && it2->y() == it->second)
      {
        it2 = _L.erase(it2);
      }
      else
      {
        ++it2;
      }
    }
  }
}
  
void Character::solve(int xy_max_c,
                      int xy_max_overall,
                      bool includeMutationEdge,
                      FrequencyMap& mapToF)
{
  // 0. obtain L
  StateGraph::IntPairSet L;
  for (CopyStateListIt it = _L.begin(); it != _L.end(); ++it)
  {
    L.insert(std::make_pair(it->x(), it->y()));
  }
  
  // 1. enumerate all state trees
  const StateGraph::StateEdgeSetSet& setS = StateGraph::getStateTrees(L, xy_max_c, includeMutationEdge);
  
  // 2. now check for each S in setS whether S is feasible
  for (StateGraph::StateEdgeSetSetIt it = setS.begin(); it != setS.end(); ++it)
  {
    const StateGraph::StateEdgeSet& S = *it;
    if (S.empty()) continue;
    
    StlRealIntervalTensor f(xy_max_overall + 1,
                            StlRealIntervalMatrix(xy_max_overall + 1,
                                                  StlRealIntervalVector(xy_max_overall + 1, std::make_pair(0, 0))));
    
    // 2a. get vertices of S
    StateGraph::CnaTripleSet verticesS;
    verticesS.insert(StateGraph::CnaTriple(1,1,0));
    
    int mut_x = -1, mut_y = -1;

    for (StateGraph::StateEdgeSetIt it2 = S.begin(); it2 != S.end(); ++it2)
    {
      const StateGraph::CnaTriple& s = it2->first;
      const StateGraph::CnaTriple& t = it2->second;
      
      verticesS.insert(s);
      verticesS.insert(t);
      
      if (s._x == t._x && s._y == t._y)
      {
        mut_x = s._x;
        mut_y = s._y;
      }
    }
    
    // 2b. determine LB and UB for VAF
    double numerator = 0;
    double denominator = 0;
    for (StateGraph::CnaTripleSetIt it2 = verticesS.begin(); it2 != verticesS.end(); ++it2)
    {
      const StateGraph::CnaTriple& triple = *it2;
      
      numerator += triple._z * _M[triple._x][triple._y];
      
      // count mutation state at most once
      if (!(triple._x == mut_x && triple._y == mut_y) || triple._z == 0)
        denominator += (triple._x + triple._y) * _M[triple._x][triple._y];
    }
    
    double LB = mut_x != -1 ? (numerator - _M[mut_x][mut_y]) / denominator : numerator / denominator;
    double UB = numerator / denominator;
    
    // 2c. stuff
    RealInterval vaf_restricted = intersect(_vafLB, _vafUB, LB, UB);
    if (!g_tol.less(vaf_restricted.second, vaf_restricted.first))
    {
      for (StateGraph::CnaTripleSetIt it2 = verticesS.begin(); it2 != verticesS.end(); ++it2)
      {
        const StateGraph::CnaTriple& triple = *it2;
        if (triple._x == mut_x && triple._y == mut_y && triple._z == 1)
        {
          f[triple._x][triple._y][triple._z].first = vaf_restricted.first * denominator - numerator + _M[mut_x][mut_y];
          f[triple._x][triple._y][triple._z].second = vaf_restricted.second * denominator - numerator + _M[mut_x][mut_y];
          
          if (!g_tol.nonZero(f[triple._x][triple._y][triple._z].first))
          {
            f[triple._x][triple._y][triple._z].first = 0;
          }
          if (!g_tol.nonZero(f[triple._x][triple._y][triple._z].second))
          {
            f[triple._x][triple._y][triple._z].second = 0;
          }
          if (!g_tol.different(f[triple._x][triple._y][triple._z].first, f[triple._x][triple._y][triple._z].second))
          {
            f[triple._x][triple._y][triple._z].second = f[triple._x][triple._y][triple._z].first;
          }
          
          assert(f[triple._x][triple._y][triple._z].first <= f[triple._x][triple._y][triple._z].second);
          
          f[triple._x][triple._y][0].first = _M[mut_x][mut_y] - f[triple._x][triple._y][triple._z].second;
          f[triple._x][triple._y][0].second = _M[mut_x][mut_y] - f[triple._x][triple._y][triple._z].first;

          assert(f[triple._x][triple._y][0].first <= f[triple._x][triple._y][0].second);
        }
        else if (triple._x != mut_x || triple._y != mut_y)
        {
          f[triple._x][triple._y][triple._z] = std::make_pair(_M[triple._x][triple._y], _M[triple._x][triple._y]);
        }
        
        if (!g_tol.nonZero(f[triple._x][triple._y][triple._z].first))
        {
          f[triple._x][triple._y][triple._z].first = 0;
        }
        if (!g_tol.nonZero(f[triple._x][triple._y][triple._z].second))
        {
          f[triple._x][triple._y][triple._z].second = 0;
        }
      }
      
      mapToF[S] = f;
    }
  }
}
  
std::ostream& operator<<(std::ostream& out, const Character& c)
{
  out << c._sampleIndex << "\t" << c._sampleLabel << "\t"
      << c._characterIndex << "\t" << c._characterLabel << "\t"
      << c._vafLB << "\t" << c._vaf << "\t" << c._vafUB;
  
  for (Character::CopyStateListIt it = c._L.begin(); it != c._L.end(); ++it)
  {
    out << "\t" << it->x() << "\t" << it->y() << "\t" << it->mu();
  }
  
  return out;
}
  
std::istream& operator>>(std::istream& in, Character& c)
{ 
  std::string line;
  gm::getline(in, line);
  
  StringVector s;
  boost::split(s, line, boost::is_any_of("\t "));
  
  try
  {
    c._sampleIndex = boost::lexical_cast<int>(s[0]);
    c._sampleLabel = s[1];
    c._characterIndex = boost::lexical_cast<int>(s[2]);
    c._characterLabel = s[3];
    c._vafLB = boost::lexical_cast<double>(s[4]);
    c._vaf = boost::lexical_cast<double>(s[5]);
    c._vafUB = boost::lexical_cast<double>(s[6]);
  }
  catch (boost::bad_lexical_cast& e)
  {
    throw std::runtime_error(getLineNumber() + "Error: " + e.what());
  }
  
  if (std::isnan(c._vafLB))
  {
    throw std::runtime_error(getLineNumber() + "Error: vafLB should not be 'nan'");
  }
  
  if (std::isnan(c._vafUB))
  {
    throw std::runtime_error(getLineNumber() + "Error: vafUB should not be 'nan'");
  }
  
  c._L.clear();
  
  int offset = 7;
  int t = (s.size() - offset) / 3;
  int max_x = -1, max_y = -1;
  for (int i = 0; i < t; ++i)
  {
    int x = boost::lexical_cast<int>(s[offset + 3*i]);
    int y = boost::lexical_cast<int>(s[offset + 3*i + 1]);
    double mu = boost::lexical_cast<double>(s[offset + 3*i + 2]);
    
    if (x > max_x) max_x = x;
    if (y > max_y) max_y = y;
    
    Character::CopyState cpyState(x, y, mu);
    
    c._L.push_back(cpyState);
  }
  
  c._M = StlDoubleMatrix(max_x + 1, StlDoubleVector(max_y + 1, 0));
  for (Character::CopyStateListIt it = c._L.begin(); it != c._L.end(); ++it)
  {
    const Character::CopyState& cpyState = *it;
    c._M[cpyState.x()][cpyState.y()] = cpyState.mu();
  }
  
  return in;
}
  
} // gm
