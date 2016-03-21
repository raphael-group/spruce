/*
 *  ccf.cpp
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
  int sampleIdx;
  double purity = 0;
  
  lemon::ArgParser ap(argc, argv);
  ap.boolOption("-version", "Show version number")
    .synonym("v", "-version")
    .refOption("p", "Purity", purity)
    .refOption("s", "Sample index", sampleIdx, true)
    .other("solution", "Solution input file");
  ap.parse();

  if (ap.given("-version"))
  {
    std::cout << "Version number: " << ANCESTREE_VERSION << std::endl;
    return 0;
  }
  
  if (ap.files().size() == 0)
  {
    std::cerr << "Error: missing input file" << std::endl;
    return 1;
  }
  
  SolutionSet sols;
  if (!read_solutions(ap.files()[0], sols))
  {
    return 1;
  }
  
  for (int idx = 0; idx < sols.solutionCount(); ++idx)
  {
    const Solution& sol = sols.solution(idx);
    PerfectPhyloTree T(sol.A(), sol.S());
    T.setLabels(sol.inferredF());
    
    if (!ap.given("p"))
    {
      purity = 1 - sol.U()(sampleIdx, 0);
    }
    
    const auto& TT = T.T();
    
    const int n = sol.inferredF().n();
    StlDoubleVector ccf(n, 0);
    StlBoolVector present(n, false);
    for (PerfectPhyloTree::NodeIt v_ci(TT); v_ci != lemon::INVALID; ++v_ci)
    {
      if (v_ci == T.root())
        continue;
      
      const IntPair& ci = T.nodeToCharState(v_ci);
      const StateTree& S_c = T.S(ci.first);
      
//      const std::string& label_c = sol.inferredF().getColLabel(ci.first);
      const std::string& label_i = S_c.label(ci.second);

      present[ci.first] = true;
      
      int x = -1, y = -1, z = -1;
      sscanf(label_i.c_str(), "(%d,%d,%d)", &x, &y, &z);
      
      if (z > 0)
      {
        ccf[ci.first] += sol.inferredF()(ci.second, sampleIdx, ci.first);
      }
      
//      std::cout << T.label(v_ci) << " " << ci.first << " " << label_c << " " << label_i << " " << z << std::endl;
    }
    
    for (int c = 0; c < n; ++c)
    {
      if (present[c])
      {
        std::cout << sol.inferredF().getColLabel(c) << "\t" << sol.inferredF().getRowLabel(sampleIdx) << "\t" << idx << "\t" << ccf[c] / purity << "\t";
        T.S(c).writeEdgeList(std::cout);
        std::cout << std::endl;
      }
    }
  }
  
  return 0;
}
