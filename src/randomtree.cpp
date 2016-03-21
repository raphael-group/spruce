/*
 *  randomtree.cpp
 *
 *   Created on: 18-oct-2015
 *       Author: M. El-Kebir
 */

#include <lemon/arg_parser.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <random>
#include "utils.h"
#include "statetree.h"
#include "stategraph.h"
#include "rootedcladisticfullancestrygraph.h"
#include "rootedcladisticsampler.h"
#include "statetreessampler.h"
#include "config.h"

using namespace gm;

typedef std::vector<StateTree> StateTreeVector;

void parse(std::istream& in, StateTreeVector& S)
{
  // parse n
  int n;
  std::string line;
  gm::getline(in, line);
  
  std::stringstream ss(line);
  ss >> n;
  
  S = StateTreeVector(n, StateTree(1));
  for (int c = 0; c < n; ++c)
  {
    in >> S[c];
  }
}

int main(int argc, char** argv)
{
  int maxXY = 4;
  int n = 10;
  int maxCopyEvents = 2;
  int seed = 0;
  
  lemon::ArgParser ap(argc, argv);
  ap.boolOption("-version", "Show version number")
    .synonym("v", "-version")
    .refOption("s", "Seed", seed, true)
    .refOption("m", "Maximum number of copy events", maxCopyEvents, true)
    .refOption("x", "Maximum #maternal copies", maxXY, true)
    .refOption("n", "Number of characters", n, true);
  ap.parse();
  
  if (ap.given("-version"))
  {
    std::cout << "Version number: " << SPRUCE_VERSION << std::endl;
    return 0;
  }
  
  StateTreesSampler ssampler(seed, maxXY, n, maxCopyEvents);
  StateTreeVector S = ssampler.sample();
  
  RootedCladisticFullAncestryGraph G(S);
  G.init();
//  G.writeDOT(std::cerr);

  RootedCladisticSampler sampler(G, seed);
  sampler.sample();
  
  std::cout << n << " #n" << std::endl;
  for (int c = 0; c < S.size(); ++c)
  {
    std::cout << S[c];
  }
  sampler.writeEdgeList(std::cout);
  
  return 0;
}
