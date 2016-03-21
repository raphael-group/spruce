/*
 * compatiblestatetrees.cpp
 *
 *  Created on: 18-feb-2016
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
#include "character.h"
#include "charactermatrix.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace gm;

int main(int argc, char** argv)
{
  int c = -1;
  lemon::ArgParser ap(argc, argv);
  std::string purity;
  bool heuristic;
  ap.boolOption("-version", "Show version number")
    .synonym("v", "-version")
    .refOption("m", "Apply heuristic", heuristic)
    .refOption("c", "Character (default: -1)", c)
    .refOption("p", "Purity vector", purity)
    .other("input", "Input file");
  ap.parse();
  
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
  
  CharacterMatrix M;
  inFile >> M;
  
  // parse purity vector
  StlDoubleVector purityValues(M.m(), 0);
  StringVector s;
  boost::split(s, purity, boost::is_any_of(" "));
  if (s.size() != M.m())
  {
    std::cerr << "Invalid purity vector size. Expected " << M.m() << " values, got " << s.size() << std::endl;
    return 1;
  }
  for (int p = 0; p < M.m(); ++p)
  {
    purityValues[p] = boost::lexical_cast<double>(s[p]);
  }

  if (c == -1)
  {
    M.init();
    if (heuristic)
    {
      M.applyHeuristic();
    }
    M.writeCompatibleStateTrees(purityValues, std::cout);
  }
  else if (!(0 <= c && c < M.n()))
  {
    std::cerr << "Incorrect c" << std::endl;
    return 1;
  }
  else
  {
    M.writeCompatibleStateTrees(std::cout, c);
  }
  
  return 0;
}