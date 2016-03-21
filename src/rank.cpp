/*
 *  rank.cpp
 *
 *   Created on: 15-oct-2015
 *       Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <cstring>
#include "utils.h"
#include "solutionset.h"

using namespace gm;

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0]
              << " <SOLUTION>" << std::endl;
    return 1;
  }
  
  SolutionSet sols;
  if (strcmp(argv[1], "-") == 0)
  {
    std::cin >> sols;
  }
  else
  {
    std::ifstream inFile(argv[1]);
    if (!inFile.good())
    {
      std::cerr << "Unable to open '" << argv[1]
                << "' for reading" << std::endl;
      return 1;
    }
    
    inFile >> sols;
    inFile.close();
  }
  
  sols.assignDistancesByOccurenceCounts();
  sols.sort();
  
  std::cout << sols;
  
  return 0;
}
