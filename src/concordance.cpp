/*
 *  overlap.cpp
 *
 *   Created on: 6-oct-2015
 *       Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <cstring>
#include "utils.h"
#include "solutionset.h"
#include "solutiongraph.h"
#include "perfectphylotree.h"

using namespace gm;

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0]
              << " <SOL_1> <TRUE_SOL>" << std::endl;
    return 1;
  }
  
 
  SolutionSet sols1;

  if (strcmp(argv[1], "-") == 0)
  {
    std::cin >> sols1;
  }
  else
  {
    std::ifstream inFile1(argv[1]);
    if (!inFile1.good())
    {
      std::cerr << "Unable to open '" << argv[1]
                << "' for reading" << std::endl;
      return 1;
    }
    inFile1 >> sols1;
    inFile1.close();
  }
  
  std::ifstream inFile2(argv[2]);
  if (!inFile2.good())
  {
    std::cerr << "Unable to open '" << argv[2]
              << "' for reading" << std::endl;
    return 1;
  }
  
  SolutionSet sols2;
  inFile2 >> sols2;
  
  inFile2.close();
  
  int n1 = sols1.solutionCount();

  PerfectPhyloTree trueTree = PerfectPhyloTree(sols2.solution(0).A(), sols2.solution(0).S());
  trueTree.setLabels(sols2.solution(0).inferredF());
  int edgeCount = trueTree.numOfEdges();
  
  StlDoubleVector recallValues(n1, 0);
  
  for (int idx1 = 0; idx1 < n1; ++idx1)
  {
    const Solution& sol1 = sols1.solution(idx1);
    PerfectPhyloTree T1(sol1.A(), sol1.S());
    T1.setLabels(sol1.inferredF());
    
    int recall = trueTree.edgeRecall(T1);
    recallValues[idx1]  = (double)recall / edgeCount;
    
    std::cout << recallValues[idx1] << std::endl;
  }
  
  std::nth_element(recallValues.begin(), recallValues.begin() + n1 / 2, recallValues.end());
  
  std::cerr << recallValues[n1/2] << std::endl;
  
  return 0;
}
