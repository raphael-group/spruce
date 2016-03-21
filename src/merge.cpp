/*
 *  merge.cpp
 *
 *   Created on: 21-feb-2016
 *       Author: M. El-Kebir
 */

#include <lemon/arg_parser.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "config.h"
#include "utils.h"
#include "solutionset.h"
#include "solutiongraph.h"
#include "perfectphylotree.h"

using namespace gm;

bool read_solutions(const std::string& filename,
                    SolutionSet& sols)
{
  if (filename == "-")
  {
    std::cin >> sols;
  }
  else
  {
    std::ifstream inFile(filename.c_str());
    if (!inFile.good())
    {
      std::cerr << "Unable to open '" << filename
                << "' for reading" << std::endl;
      return false;
    }
    
    inFile >> sols;
    inFile.close();
  }
  
  return true;
}

int main(int argc, char** argv)
{
  bool retainAll = false;
  lemon::ArgParser ap(argc, argv);
  ap.boolOption("-version", "Show version number")
    .synonym("v", "-version")
    .refOption("m", "Retain all trees (not only the largest)", retainAll)
    .other("solution", "Solution files");
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
  
  int treeSize = 0;
  SolutionSet allSols;
  for (const std::string& filename : ap.files())
  {
    SolutionSet newSols;
    read_solutions(filename, newSols);
    if (newSols.solutionCount() > 0)
    {
      const Solution& sol0 = newSols.solution(0);
      PerfectPhyloTree T(sol0.A(), sol0.S());
      
      int size_sol0 = T.numOfVertices();

      if (retainAll)
      {
        allSols.add(newSols);
      }
      else
      {
        if (size_sol0 >= treeSize)
        {
          if (size_sol0 > treeSize)
          {
            treeSize = size_sol0;
            allSols.clear();
          }

          allSols.add(newSols);
        }
      }
    }
  }
  
  allSols.unique();
  
  std::cout << allSols;
  
  return 0;
}
