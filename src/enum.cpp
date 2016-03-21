/*
 * enum.cpp
 *
 *  Created on: 29-sep-2015
 *      Author: M. El-Kebir
 */

#include <lemon/arg_parser.h>
#include <iostream>
#include <fstream>
#include "utils.h"
#include "config.h"
#include "rootedcladisticancestrygraph.h"
#include "rootedcladisticnoisyancestrygraph.h"
#include "rootedcladisticenumeration.h"
#include "rootedcladisticnoisyenumeration.h"
#include "cnaenumerate.h"
#include "noisycnaenumerate.h"
#include "character.h"
#include "charactermatrix.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

using namespace gm;

typedef RootedCladisticAncestryGraph::StateTreeVector StateTreeVector;

bool parse_purity_vector(const std::string& purityString,
                         int m,
                         StlDoubleVector& purityValues)
{
  // parse purity vector
  purityValues = StlDoubleVector(m, 0);
  StringVector s;
  boost::split(s, purityString, boost::is_any_of(" "));
  if (s.size() != m)
  {
    std::cerr << "Invalid purity vector size. Expected " << m << " values, got " << s.size() << std::endl;
    return false;
  }
  
  for (int p = 0; p < m; ++p)
  {
    purityValues[p] = boost::lexical_cast<double>(s[p]);
  }
  return true;
}

void read_cladistic_tensor(std::ifstream& inFile,
                           RealTensor& F,
                           StateTreeVector& S)
{
  inFile >> F;
  F.setLabels(inFile);
  std::string line;
  gm::getline(inFile, line);
  
  S = StateTreeVector(F.n(), StateTree(F.k()));
  
  for (int c = 0; c < F.n(); ++c)
  {
    inFile >> S[c];
  }
}

void read_cladistic_tensor_noisy(std::ifstream& inFile,
                                 RealTensor& F,
                                 StateTreeVector& S,
                                 RealTensor& F_lb,
                                 RealTensor& F_ub)
{
  read_cladistic_tensor(inFile, F, S);
  
  std::string line;
  gm::getline(inFile, line);
  inFile >> F_lb;
  gm::getline(inFile, line);
  inFile >> F_ub;
}

void enumerate(int limit,
               int timeLimit,
               int threads,
               int state_tree_limit,
               std::ifstream& inFile,
               const std::string& intervalFile,
               int lowerbound,
               bool monoclonal,
               const std::string& purityString,
               bool writeCliqueFile,
               bool readCliqueFile,
               const std::string& cliqueFile,
               bool shuffle,
               int offset,
               const IntSet& whiteList,
               SolutionSet& sols)
{
  CharacterMatrix M;
  inFile >> M;
  
  if (!intervalFile.empty())
  {
    std::ifstream inFile2(intervalFile.c_str());
    if (inFile2.good())
    {
      M.setIntervals(inFile2);
    }
  }
  
  StlDoubleVector purityValues;
  if (purityString != "" && !parse_purity_vector(purityString, M.m(), purityValues))
  {
    return;
  }
  if (monoclonal && purityString == "")
  {
    std::cerr << "Expected purity values" << std::endl;
    return;
  }
  
  std::cerr << "Intializing copy-state matrix ..." << std::endl;
  M.init();
  
  if (monoclonal)
  {
    M.applyHeuristic();
  }
  
  CompatibilityGraph* pComp = NULL;
  if (!readCliqueFile)
  {
    std::cerr << "Initializing compatibility graph ..." << std::endl;
    pComp = new CompatibilityGraph(M, shuffle);
  }
  else
  {
    std::ifstream inFile(cliqueFile);
    pComp = new CompatibilityGraph(M, shuffle, inFile);
  }
  
  if (writeCliqueFile)
  {
    std::ofstream outFile(cliqueFile);
    if (outFile.good())
    {
      pComp->write(outFile);
    }
  }
  
  NoisyCnaEnumerate alg(M, purityValues, *pComp, lowerbound);
  alg.init(state_tree_limit);
  alg.enumerate(limit, timeLimit, threads, state_tree_limit, monoclonal, offset, whiteList);

  sols = alg.sols();
  delete pComp;
  
  std::cerr << "Generated " << sols.solutionCount() << " solutions" << std::endl;
  inFile.close();
}

int main(int argc, char** argv)
{
  int limit = -1;
  int timeLimit = -1;
  int threads = 2;
  int state_tree_limit = -1;
  int random_seed = 0;
  int lowerbound = 0;
  bool monoclonal = false;
  bool deterministic = false;
  std::string purityString;
  std::string cliqueFile;
  int offset = 0;
  int verbosityLevel = 1;
  std::string whiteListString;
  
  lemon::ArgParser ap(argc, argv);
  ap.boolOption("-version", "Show version number")
    .refOption("v", "Verbosity level (default: 1)", verbosityLevel)
    .boolOption("p", "Perfect data")
    .boolOption("c", "Cladistic characters")
    .refOption("m", "Monoclonal", monoclonal)
    .refOption("purity", "Purity values (used for fixing trunk)",  purityString)
    .refOption("clique", "Clique file", cliqueFile)
    .refOption("d", "Deterministic enumeration of state tree combinations", deterministic)
    .refOption("t", "Number of threads (default: 2)", threads)
    .refOption("l", "Maximum number of trees to enumerate (default: -1)", limit)
    .refOption("ll", "Time limit in seconds (default: -1)", timeLimit)
    .refOption("s", "State tree combination limit (default: -1)", state_tree_limit)
    .refOption("o", "State tree combination offset (default: 0)", offset)
    .refOption("r", "Seed", random_seed)
    .refOption("lb", "Lowerbound (default: 0)", lowerbound)
    .refOption("w", "Characters that must be present in the solution trees", whiteListString)
    .other("input_1", "Input file")
    .other("input_2", "Interval file relating SNVs affected by the same CNA");
  ap.parse();

  g_rng = std::mt19937(random_seed);
  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);
  
  if (ap.given("-version"))
  {
    std::cout << "Version number: " << SPRUCE_VERSION << std::endl;
    return 0;
  }
  
  if (ap.files().size() == 0)
  {
    std::cerr << "Error: missing input file" << std::endl;
    return 1;
  }
  
  std::ifstream inFile(ap.files()[0].c_str());
  if (!inFile.good())
  {
    std::cerr << "Unable to open '" << ap.files()[0].c_str()
              << "' for reading" << std::endl;
    return 1;
  }
  
  IntSet whiteList;
  if (!whiteListString.empty())
  {
    StringVector s;
    boost::split(s, whiteListString, boost::is_any_of(", ;"));
    for (const std::string& str : s)
    {
      whiteList.insert(boost::lexical_cast<int>(str));
    }
  }
  
  bool readCliqueFile = false;
  bool writeCliqueFile = false;
  if (cliqueFile != "" && boost::filesystem::exists(cliqueFile))
  {
    readCliqueFile = true;
  }
  else if (cliqueFile != "")
  {
    writeCliqueFile = true;
  }
  
  SolutionSet sols;
  enumerate(limit, timeLimit, threads,
            state_tree_limit, inFile,
            (ap.files().size() > 1 ? ap.files()[1] : ""),
            lowerbound,
            monoclonal,
            purityString,
            writeCliqueFile,
            readCliqueFile,
            cliqueFile,
            !deterministic,
            offset,
            whiteList,
            sols);
  std::cout << sols;
  
  return 0;
}
