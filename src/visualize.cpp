/*
 *  visualize.cpp
 *
 *   Created on: 1-oct-2015
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
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace gm;

typedef PerfectPhyloTree::StringToStringMap StringToStringMap;

bool parse_color_mapping(const std::string& filename,
                         StringToStringMap& label2color)
{
  std::ifstream inFile(filename.c_str());
  if (!inFile.good())
  {
    std::cerr << "Unable to open '" << filename
              << "' for reading" << std::endl;
    return false;
  }
  
  std::string line;
  StringVector s;
  while (inFile.good())
  {
    gm::getline(inFile, line);
    boost::split(s, line, boost::is_any_of("\t"));
    
    if (s.size() == 2)
    {
      label2color[s[0]] = s[1];
    }
  }
  
  return true;
}

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
  bool showStateVectors = false;
  bool showCumFreqs = false;
  bool showEdgeLabels = false;
  bool summarize = false;
  bool json = false;
  bool ascii = false;
  int solIdx = -1;
  bool usageMatrix = false;
  bool simple = false;
  std::string labelToColorFilename;
  
  lemon::ArgParser ap(argc, argv);
  ap.boolOption("-version", "Show version number")
    .synonym("v", "-version")
    .refOption("-showStateVectors", "Show state vectors", showStateVectors)
    .refOption("-showCumFreqs", "Show cumulative frequencies", showCumFreqs)
    .refOption("-showEdgeLabels", "Show edge labels", showEdgeLabels)
    .refOption("c", "Label to color file", labelToColorFilename)
    .refOption("u", "Print usage matrix", usageMatrix)
    .refOption("i", "Solution index (default: -1)", solIdx)
    .refOption("s", "Summarize", summarize)
    .refOption("simple", "Simple visualization", simple)
    .refOption("j", "Output JSON", json)
    .refOption("a", "ASCII", ascii)
    .other("solution", "Solution input file")
    .other("true_solution", "True solution file");
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
  
  SolutionSet sols, trueSols;
  if (!read_solutions(ap.files()[0], sols))
  {
    return 1;
  }
  if (ap.files().size() == 2)
  {
    if (!read_solutions(ap.files()[1], trueSols))
    {
      return 1;
    }
  }
  
  StringToStringMap label2color;
  if (!labelToColorFilename.empty() && !parse_color_mapping(labelToColorFilename, label2color))
  {
    return 1;
  }
  
  if (summarize && trueSols.solutionCount() > 0)
  {
    PerfectPhyloTree trueT(trueSols.solution(0).A(), trueSols.solution(0).S());
    trueT.setLabels(trueSols.solution(0).inferredF());
    sols.writeSummaryDOT(std::cout, trueT);
  }
  else if (summarize || json)
  {
    if (json)
      sols.writeSummaryJSON(std::cout, label2color);
    else
      sols.writeSummaryDOT(std::cout, label2color);
  }
  else
  {
    bool all = false;
    if (!(0 <= solIdx && solIdx < sols.solutionCount()))
    {
      std::cerr << "Invalid solution index, will visualize all" << std::endl;
      all = true;
    }
    
    if (all)
    {
      for (int idx = 0; idx < sols.solutionCount(); ++idx)
      {
        const Solution& sol = sols.solution(idx);
        PerfectPhyloTree T(sol.A(), sol.S());
        
        SolutionGraph solGraph(sol.inferredF(), sol.S(), sol, 0,
                               showStateVectors, showCumFreqs, showEdgeLabels);
        
        char buf[1024];
//        snprintf(buf, 1024, "%s.sol%d.dot", fs::basename(fs::path(ap.files()[0].c_str())).c_str(), idx);
        snprintf(buf, 1024, "sol%d.dot", idx);
        
        std::ofstream out(buf);
        solGraph.writeDOT(out);
      }
    }
    else
    {
      const Solution& sol = sols.solution(solIdx);
      PerfectPhyloTree T(sol.A(), sol.S());

      SolutionGraph solGraph(sol.inferredF(), sol.S(), sol, 0.05,
                             showStateVectors, showCumFreqs, showEdgeLabels);
      
      if (usageMatrix)
      {
        solGraph.writeUsageMatrix(std::cout);
      }
      else if (!ascii)
      {
        if (simple)
        {
          T.writeDOT(sol.inferredF(), sol.U(), label2color, std::cout);
        }
        else
        {
          solGraph.writeDOT(std::cout);
        }
      }
      else
      {
        solGraph.writeASCII(std::cout);
      }
    }
  }
  
  return 0;
}
