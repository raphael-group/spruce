/* 
 * utils.h
 *
 *  Created on: 18-jun-2015
 *      Author: M. El-Kebir
 */

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <string>
#include <map>
#include <set>
#include <istream>
#include <lemon/tolerance.h>
#include <lemon/list_graph.h>
#include <stdlib.h>
#include <stdio.h>
#include <random>

#undef DEBUG

namespace gm {
  
typedef lemon::ListGraph Graph;
typedef lemon::ListDigraph Digraph;
typedef Digraph::NodeMap<std::string> StringNodeMap;
typedef std::map<std::string, Digraph::Node> StringToNodeMap;

typedef enum
{
  VERBOSE_NONE,
  VERBOSE_ESSENTIAL,
  VERBOSE_NON_ESSENTIAL,
  VERBOSE_DEBUG
} VerbosityLevel;

/// Verbosity level
extern VerbosityLevel g_verbosity;

/// Line number used for error handling in file I/O
extern int g_lineNumber;
  
/// Return current line number
std::string getLineNumber();
  
/// Random number generator
extern std::mt19937 g_rng;
  
typedef std::vector<bool> StlBoolVector;
typedef std::vector<StlBoolVector> StlBoolMatrix;
  
typedef std::vector<int> StlIntVector;
typedef StlIntVector::const_iterator StlIntVectorIt;
typedef std::vector<StlIntVector> StlIntMatrix;
typedef StlIntMatrix::const_iterator StlIntMatrixIt;
  
typedef std::vector<double> StlDoubleVector;
typedef std::vector<StlDoubleVector> StlDoubleMatrix;
typedef std::vector<StlDoubleMatrix> StlDoubleTensor;
  
typedef std::vector<std::string> StringVector;
  
typedef std::pair<double, double> RealInterval;
typedef std::vector<RealInterval> StlRealIntervalVector;
typedef std::vector<StlRealIntervalVector> StlRealIntervalMatrix;
typedef std::vector<StlRealIntervalMatrix> StlRealIntervalTensor;
  
typedef std::pair<int, int> IntPair;
typedef std::set<IntPair> IntPairSet;
  
bool intPairCompare(const IntPair& a, const IntPair& b);
  
typedef std::set<int> IntSet;
typedef IntSet::const_iterator IntSetIt;
  
std::ostream& operator<<(std::ostream& out, const StlBoolMatrix& M);
std::istream& operator>>(std::istream& in, StlBoolMatrix& M);

std::ostream& operator<<(std::ostream& out, const StlDoubleMatrix& M);
std::istream& operator>>(std::istream& in, StlDoubleMatrix& M);
  
std::ostream& operator<<(std::ostream& out, const StlDoubleTensor& M);
std::istream& operator>>(std::istream& in, StlDoubleTensor& M);
  
std::ostream& operator<<(std::ostream& out, const StlRealIntervalMatrix& M);
std::istream& operator>>(std::istream& in, StlRealIntervalMatrix& M);

StlBoolVector discretize(const StlDoubleVector& v);

StlBoolMatrix discretize(const StlDoubleMatrix& M);
  
std::istream& getline(std::istream& is, std::string& t);
  
extern lemon::Tolerance<double> g_tol;

} // namespace gm

#endif /* UTILS_H */
