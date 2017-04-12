/*
 * cliques.cpp
 *
 *  Created on: 27-jun-2016
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "config.h"
#include "charactermatrix.h"
#include "compatibilitygraph.h"
#include <fstream>
#include <lemon/arg_parser.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace gm;

int main(int argc, char** argv)
{
  int size = -1;
  std::string filterString;
  int verbosityLevel = 1;
  int cliqueLimit;
  
  lemon::ArgParser ap(argc, argv);
  ap.boolOption("-version", "Show version number")
    .refOption("v", "Verbosity level (default: 1)", verbosityLevel)
    .refOption("s", "Maximal clique size (default: -1 (maximum))", size)
    .refOption("f", "Filter (default : \"\")", filterString)
    .refOption("l", "Clique limit (default : -1 (unlimited))", cliqueLimit)
    .other("input", "Input file");
  ap.parse();
  
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
  
  IntPairSet filter;
  IntSet whiteList;
  if (!filterString.empty())
  {
    StringVector s;
    boost::split(s, filterString, boost::is_any_of(";"));
    for (const std::string& str : s)
    {
      StringVector ss;
      boost::split(ss, str, boost::is_any_of(","));
      filter.insert(std::make_pair(boost::lexical_cast<int>(ss[0]), boost::lexical_cast<int>(ss[1])));
    }
  }

  CharacterMatrix M;
  inFile >> M;
  
  M.init();
  M.applyHeuristic();
  
  CompatibilityGraph G(M);
  G.setCliqueLimit(cliqueLimit);
  G.init(filter, size);
  
  G.write(std::cout);
  
  return 0;
}
