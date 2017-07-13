/*
 * charactermatrix.cpp
 *
 *  Created on: 15-jan-2016
 *      Author: M. El-Kebir
 */

#include "charactermatrix.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace gm {

CharacterMatrix::CharacterMatrix()
  : _m(0)
  , _n(0)
  , _M()
  , _F()
  , _maxX()
  , _maxXY(0)
  , _tripleToState()
  , _stateToTriple()
  , _stateLabel()
  , _charToInterval()
  , _intervals()
{
}
  
void CharacterMatrix::applyHeuristic()
{
  // make groups
  typedef std::set<CnaTriple> CnaTripleSet;
  typedef std::map<CnaTripleSet, IntSet> CnaTripleSetMap;
  typedef std::vector<CnaTripleSetMap> CnaTripleSetMapVector;
  
  CnaTripleSetMapVector stateTreeVertexSet(_n);
  for (int c = 0; c < _n; ++c)
  {
    int idx = 0;
    int minL1 = std::numeric_limits<int>::max();
    for (Character::FrequencyMapIt it = _F[0][c].begin(); it != _F[0][c].end(); ++it, ++idx)
    {
      const StateGraph::StateEdgeSet& edgeSet = it->first;

      CnaTripleSet V;
      int L1 = 0;
      for (const StateGraph::StateEdge& e : edgeSet)
      {
        V.insert(e.first);
        V.insert(e.second);
        
        L1 += abs(e.first._x - e.second._x) + abs(e.first._y - e.second._y);
      }
      if (L1 < minL1)
      {
        stateTreeVertexSet[c][V].clear();
        stateTreeVertexSet[c][V].insert(idx);
        minL1 = L1;
      }
      else if (L1 == minL1)
      {
        stateTreeVertexSet[c][V].insert(idx);
      }
    }
    
    IntSet retain;
    for (const auto& V : stateTreeVertexSet[c])
    {
      for (auto kv : V.second)
      {
        retain.insert(kv);
      }
    }
    
    for (int p = 0; p < _m; ++p)
    {
      idx = 0;
      Character::FrequencyMap newF;
      for (Character::FrequencyMapIt it = _F[p][c].begin(); it != _F[p][c].end(); ++idx, ++it)
      {
        if (retain.count(idx))
        {
          newF[it->first] = it->second;
        }
        else
        {
          if (p == 0)
          {
            if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
            {
              std::cerr << "Removing character " <<  _M[0][c].characterLabel() << " (" << c << ") : ";
              StateGraph::print(it->first, std::cerr);
              std::cerr << std::endl;
            }
          }
        }
      }
      _F[p][c] = newF;
    }
  }
}
  
void CharacterMatrix::setIntervals(std::istream& in)
{
  typedef std::vector<std::string> StringVector;

  std::string line;
  
  IntSetSetVector X(m());
  for (int p = 0; p < m(); ++p)
  {
    gm::getline(in, line);
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t"));
    
    for (int i = 0; i < s.size(); ++i)
    {
      StringVector ss;
      boost::split(ss, s[i], boost::is_any_of(" "));
      
      IntSet S;
      for (int j = 0; j < ss.size(); ++j)
      {
        try
        {
          int c = boost::lexical_cast<int>(ss[j]);
          S.insert(c);
        }
        catch (boost::bad_lexical_cast& e)
        {
          throw std::runtime_error(getLineNumber() + "Error: " + e.what());
        }
      }
      X[p].insert(S);
    }
  }
  
  IntSetVector charToSet(n());
  IntSetSet XX = X[0];
  
  // partition refinement
  for (int p = 1; p < m(); ++p)
  {
    const IntSetSet& X_p = X[p];
    
    // this will be the new partition, i.e. XX refined with X_p
    IntSetSet XXX;
    for (IntSetSetIt it1 = X_p.begin(); it1 != X_p.end(); ++it1)
    {
      const IntSet& S = *it1;
      
      for (IntSetSetIt it2 = XX.begin(); it2 != XX.end(); ++it2)
      {
        const IntSet& SS = *it2;

        /// SSS = SS intersection S
        IntSet SSS;
        std::set_intersection(SS.begin(), SS.end(), S.begin(), S.end(), std::inserter(SSS, SSS.begin()));

        // is S contained in SS?
        if (!SSS.empty())
        {
          XXX.insert(SSS);
        }
      }
    }
    XX = XXX;
  }
  
  _charToInterval = IntSetVector(n());
  for (IntSetSetIt it2 = XX.begin(); it2 != XX.end(); ++it2)
  {
    const IntSet& SS = *it2;
    for (IntSetIt it3 = SS.begin(); it3 != SS.end(); ++it3)
    {
      _charToInterval[*it3] = SS;
    }
  }
  
  _intervals = XX;
}
  
StateTree CharacterMatrix::stateTree(int s, int c) const
{
  assert(0 <= c && c < _n);
  assert(0 <= s && s < _F[0][c].size());

  const int k = _stateToTriple.size();
  
  Character::FrequencyMapIt it = _F[0][c].begin();
  for (int i = 0; i < s; ++i) ++it;
  const StateGraph::StateEdgeSet& edgeSet = it->first;
  
  StlIntVector pi(_stateToTriple.size(), -2);
  pi[0] = -1;
  
  for (StateGraph::StateEdgeSetIt edgeIt = edgeSet.begin(); edgeIt != edgeSet.end(); ++edgeIt)
  {
    assert(_tripleToState.find(edgeIt->first) != _tripleToState.end());
    assert(_tripleToState.find(edgeIt->second) != _tripleToState.end());
    pi[_tripleToState.find(edgeIt->second)->second] = _tripleToState.find(edgeIt->first)->second;
  }
  
  StateTree S(pi);
  for (int i = 0; i < k; ++i)
  {
    S.setLabel(i, _stateLabel[i]);
  }
  return S;
}
  
StlRealIntervalVector CharacterMatrix::ccf(const StlDoubleVector& purityValues,
                                           int s, int c) const
{
  StlRealIntervalVector res(_m, RealInterval(0, 0));
  
  const int kk = k();
  for (int i = 0; i < kk; ++i)
  {
    const auto& xyz = _stateToTriple[i];
    if (xyz._z > 0)
    {
      for (int p = 0; p < _m; ++p)
      {
        res[p].first += get_f_lb(s, p, c, i);
        res[p].second += get_f_ub(s, p, c, i);
      }
    }
  }
  
  for (int p = 0; p < _m; ++p)
  {
    res[p].first /= purityValues[p];
    res[p].second /= purityValues[p];
  }
  
  return res;
}
  
CharacterMatrix::IntPairPairSet CharacterMatrix::copyTree(int s, int c) const
{
  IntPairPairSet res;
  
  assert(0 <= c && c < _n);
  assert(0 <= s && s < _F[0][c].size());
  
  Character::FrequencyMapIt it = _F[0][c].begin();
  for (int i = 0; i < s; ++i) ++it;
  const StateGraph::StateEdgeSet& edgeSet = it->first;

  for (StateGraph::StateEdgeSetIt edgeIt = edgeSet.begin(); edgeIt != edgeSet.end(); ++edgeIt)
  {
    const StateGraph::CnaTriple& xyz_1 = edgeIt->first;
    const StateGraph::CnaTriple& xyz_2 = edgeIt->second;
    
    if (xyz_1._x == xyz_2._x && xyz_1._y == xyz_2._y)
    {
      assert(xyz_1._z == 0);
      assert(xyz_2._z == 1);
      
      // skip
    }
    else
    {
      res.insert(std::make_pair(std::make_pair(xyz_1._x, xyz_1._y), std::make_pair(xyz_2._x, xyz_2._y)));
    }
  }
  
  return res;
}
  
void CharacterMatrix::writeCompatibleStateTrees(std::ostream& out,
                                                int c)
{
  _F = FrequencyMapMatrix(_m, FrequencyMapVector(_n));
  _maxX = StlIntVector(_n, 0);
  _maxXY = 0;
  bool includeMutationEdge_c = false;
  
  // remove (x,y) if mu_(x,y) = 0 across all samples p
  StateGraph::IntPairSet L0 = _M[0][c].L0();
  for (StateGraph::IntPairSetNonConstIt it = L0.begin(); it != L0.end();)
  {
    bool remove = false;
    for (int p = 1; p < _m; ++p)
    {
      if (_M[p][c].mu(it->first, it->second) != 0)
      {
        remove = true;
      }
    }
    
    // don't remove (1,1)
    if (remove || *it == std::make_pair(1, 1))
    {
      it = L0.erase(it);
    }
    else
    {
      ++it;
    }
  }
  
  // determine maxX
  for (int p = 0; p < _m; ++p)
  {
    includeMutationEdge_c = includeMutationEdge_c || (_M[p][c].vafUB() > 0);
    
    _M[p][c].remove(L0);
    const Character::CopyStateList & L = _M[p][c].L();
    for (Character::CopyStateListIt it = L.begin(); it != L.end(); ++it)
    {
      if (_maxX[c] < it->x())
      {
        _maxX[c] = it->x();
      }
      if (_maxXY < it->x())
      {
        _maxXY = it->x();
      }
    }
  }
  
  // compute feasible state trees
  for (int p = 0; p < _m; ++p)
  {
    _M[p][c].solve(_maxX[c], _maxXY, includeMutationEdge_c, _F[p][c]);
  }
  
  typedef std::map<StateGraph::StateEdgeSet, IntSet> CountMap;
  typedef CountMap::const_iterator CountMapIt;
  
  CountMap X;
  for (int p = 0; p < _m; ++p)
  {
    const Character::FrequencyMap& F_pc = _F[p][c];
    for (Character::FrequencyMapIt it = F_pc.begin(); it != F_pc.end(); ++it)
    {
      const StateGraph::StateEdgeSet& S = it->first;
      if (X.find(S) == X.end())
      {
        X[S] = IntSet();
        X[S].insert(p);
      }
      else
      {
        X[S].insert(p);
      }
    }
  }
  
  for (CountMapIt it = X.begin(); it != X.end(); ++it)
  {
    const StateGraph::StateEdgeSet& S = it->first;
    const IntSet& P = it->second;
    
    out << P.size() << "\t[ ";
    for (IntSetIt it2 = P.begin(); it2 != P.end(); ++it2)
    {
      out << *it2 << " ";
    }
    out << "]\t";
    
    StateGraph::print(S, out);
    
    out << std::endl;
  }
}
  
void CharacterMatrix::writeCompatibleStateTrees(const StlDoubleVector& purityValues,
                                                std::ostream& out)
{
  out << "c\tlabel\tstate tree";
  for (int p = 0; p < _m; ++p)
  {
    out << "\t" << _M[p][0].sampleLabel() << "\t" << _M[p][0].sampleLabel();
  }
  out << "\t" << "#samples clonal" << std::endl;
  
  for (int c = 0; c < _n; ++c)
  {
    for (int s = 0; s < numStateTrees(c); ++s)
    {
      Character::FrequencyMapIt it = _F[0][c].begin();
      for (int i = 0; i < s; ++i) ++it;
      
      out << c << "\t" << _M[0][c].characterLabel() << "\t";
      StateGraph::print(it->first , out);
      
      int clonal = 0;
      StlRealIntervalVector ccfs = ccf(purityValues, s, c);
      for (int p = 0; p < _m; ++p)
      {
        out << "\t" << ccfs[p].first << "\t" << ccfs[p].second;
        if (ccfs[p].second >= 1)
        {
          ++clonal;
        }
      }
      
      out << "\t" << clonal << std::endl;
    }
  }
}
  
void CharacterMatrix::init()
{
  _F = FrequencyMapMatrix(_m, FrequencyMapVector(_n));
  _maxX = StlIntVector(_n, 0);
  _maxXY = 0;

  StlBoolVector includeMutationEdge(_n, false);
  
  for (int c = 0; c < _n; ++c)
  {
    // remove (x,y) if mu_(x,y) = 0 across all samples p
    StateGraph::IntPairSet L0 = _M[0][c].L0();
    for (StateGraph::IntPairSetNonConstIt it = L0.begin(); it != L0.end();)
    {
      bool remove = false;
      for (int p = 1; p < _m; ++p)
      {
        if (_M[p][c].mu(it->first, it->second) != 0)
        {
          remove = true;
        }
      }
      
      // don't remove (1,1)
      if (remove || *it == std::make_pair(1, 1))
      {
        it = L0.erase(it);
      }
      else
      {
        ++it;
      }
    }

    // determine maxX
    for (int p = 0; p < _m; ++p)
    {
      includeMutationEdge[c] = includeMutationEdge[c] || (_M[p][c].vafUB() > 0);
      
      _M[p][c].remove(L0);
      const Character::CopyStateList & L = _M[p][c].L();
      for (Character::CopyStateListIt it = L.begin(); it != L.end(); ++it)
      {
        if (_maxX[c] < it->x())
        {
          _maxX[c] = it->x();
        }
        if (_maxXY < it->x())
        {
          _maxXY = it->x();
        }
      }
    }
  }
  
  // compute feasible state trees
  if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
  {
    std::cerr << std::endl << "Enumerating compatible state trees for each character ..." << std::endl;
  }
  for (int p = 0; p < _m; ++p)
  {
    for (int c = 0; c < _n; ++c)
    {
      if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
      {
        std::cerr << "Generating compatible state trees for character " << _M[p][c].characterLabel() << " (" << c
                  << ") in sample " << _M[p][c].sampleLabel() << " (" << p << ") ..." << std::flush;
      }

      _M[p][c].solve(_maxX[c], _maxXY, includeMutationEdge[c], _F[p][c]);
      
      if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
      {
        std::cerr << " Done: " << _F[p][c].size() << " state trees" << std::endl;
      }
    }
  }
  
  // compute intersection across samples for each character
  for (int c = 0; c < _n; ++c)
  {
    StateGraph::StateEdgeSetSet intersection;
    const Character::FrequencyMap& F_0c = _F[0][c];
    for (Character::FrequencyMapIt it = F_0c.begin(); it != F_0c.end(); ++it)
    {
      intersection.insert(it->first);
    }
    
    for (int p = 1; p < _m; ++p)
    {
      const Character::FrequencyMap& F_pc = _F[p][c];
      for (StateGraph::StateEdgeSetSetNonConstIt it = intersection.begin(); it != intersection.end();)
      {
        const StateGraph::StateEdgeSet& S = *it;
        if (F_pc.find(S) == F_pc.end())
        {
          it = intersection.erase(it);
        }
        else
        {
          ++it;
        }
      }
    }
    
    // now update state trees
    for (int p = 0; p < _m; ++p)
    {
      Character::FrequencyMap& F_pc = _F[p][c];
      for (Character::FrequencyMapNonConstIt it = F_pc.begin(); it != F_pc.end();)
      {
        const StateGraph::StateEdgeSet& S = it->first;
        if (intersection.find(S) == intersection.end())
        {
          it = F_pc.erase(it);
        }
        else
        {
          ++it;
        }
      }
    }
  }
  
//  for (int p = 0; p < _m; ++p)
//  {
//    for (int c = 0; c < _n; ++c)
//    {
//      std::cout << "State trees for (" << p << "," << c << "):" << std::endl;
//      for (Character::FrequencyMapIt it = _F[p][c].begin(); it != _F[p][c].end(); ++it)
//      {
//        StateGraph::print(it->first, std::cout);
//      }
//      std::cout << std::endl;
//    }
//  }
  
  // update state map
  _tripleToState.clear();
  _stateToTriple.clear();
  
  StateGraph::CnaTriple triple110(1,1,0);
  _stateToTriple.push_back(triple110);
  _tripleToState[triple110] = 0;
  
  for (int c = 0; c < _n; ++c)
  {
    const Character::FrequencyMap& f_0c = _F[0][c];
    for (Character::FrequencyMapIt it = f_0c.begin(); it != f_0c.end(); ++it)
    {
      const StateGraph::StateEdgeSet& S = it->first;
      for (StateGraph::StateEdgeSetIt it2 = S.begin(); it2 != S.end(); ++it2)
      {
        if (_tripleToState.find(it2->first) == _tripleToState.end())
        {
          _tripleToState[it2->first] = _stateToTriple.size();
          _stateToTriple.push_back(it2->first);
        }
        if (_tripleToState.find(it2->second) == _tripleToState.end())
        {
          _tripleToState[it2->second] = _stateToTriple.size();
          _stateToTriple.push_back(it2->second);
        }
      }
    }
  }
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    for (int c = 0; c < _n; ++c)
    {
      std::cerr << "Number of state trees for character " << _M[0][c].characterLabel() << " (" << c << ") : " << numStateTrees(c) << std::endl;
    }
  }
  
  const int kk = k();
  _stateLabel = StringVector(kk);
  char buf[1024];
  for (int i = 0; i < kk; ++i)
  {
    const StateGraph::CnaTriple& triple = _stateToTriple[i];
    snprintf(buf, 1024, "(%d,%d,%d)", triple._x, triple._y, triple._z);
    _stateLabel[i] = buf;
  }
}
  
std::ostream& operator<<(std::ostream& out, const CharacterMatrix& M)
{
  out << M._m << " # m" << std::endl;
  out << M._n << " # n" << std::endl;
  
  for (int p = 0; p < M._m; ++p)
  {
    for (int c = 0; c < M._n; ++c)
    {
      out << M._M[p][c] << std::endl;
    }
  }
  
  return out;
}
  
std::istream& operator>>(std::istream& in, CharacterMatrix& M)
{
  g_lineNumber = 0;
  
  std::string line;
  gm::getline(in, line);
  
  int m = -1;
  int n = -1;
  
  std::stringstream ss(line);
  ss >> m;
  
  if (m <= 0)
  {
    throw std::runtime_error(getLineNumber() + "Error: m should be nonnegative");
  }
  
  gm::getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n;
  
  if (n <= 0)
  {
    throw std::runtime_error(getLineNumber() + "Error: n should be nonnegative");
  }
  
  M._m = m;
  M._n = n;
  M._M = CharacterMatrix::StlCharacterMatrix(m, CharacterMatrix::StlCharacterVector(n));
  StlBoolMatrix present(m, StlBoolVector(n, false));
  
  while (in.good())
  {
    gm::getline(in, line);
    if (line == "" || line[0] == '#')
      continue;
    
    ss.clear();
    ss.str(line);
    
    Character character;
    ss >> character;
    
    int c = character.characterIndex();
    int p = character.sampleIndex();
    
    if (!(0 <= p && p < m) || !(0 <= c && c < n))
    {
      throw std::runtime_error(getLineNumber() + "Error: invalid character ("
                               + boost::lexical_cast<std::string>(p)
                               + ","
                               + boost::lexical_cast<std::string>(c)
                               + ")");
    }
    
    if (present[p][c])
    {
      throw std::runtime_error(getLineNumber() + "Error: duplicate character ("
                               + boost::lexical_cast<std::string>(p)
                               + ","
                               + boost::lexical_cast<std::string>(c)
                               + ")");
    }
    
    M._M[p][c] = character;
    present[p][c] = true;
  }
  
  for (int p = 0; p < m; ++p)
  {
    for (int c = 0; c < n; ++c)
    {
      if (!present[p][c])
      {
        throw std::runtime_error(getLineNumber() + "Error: missing character ("
                                 + boost::lexical_cast<std::string>(p)
                                 + ","
                                 + boost::lexical_cast<std::string>(c)
                                 + ")");
      }
    }
  }
  
  return in;
}
  
} // namespace gm
