/*
 * dimacsparser.h
 *
 *  Created on: 19-nov-2014
 *      Author: M. El-Kebir
 */

#ifndef DIMACSPARSER_H
#define DIMACSPARSER_H

#include <assert.h>
#include <string.h>
#include <map>
#include <string>
#include <limits>
#include <lemon/core.h>
#include "parser.h"
#include "utils.h"

namespace nina {
namespace mwcs {

template<typename GR>
class DimacsParser : public Parser<GR>
{
public:
  /// Graph type
  typedef GR Graph;
  /// Base class type
  typedef Parser<GR> Parent;

private:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

public:
  typedef typename Parent::InvIdNodeMap InvIdNodeMap;
  typedef typename Parent::IdNodeMap IdNodeMap;
  typedef typename Parent::WeightNodeMap WeightNodeMap;
  typedef typename Parent::WeightEdgeMap WeightEdgeMap;

  using Parent::safeGetline;
  using Parent::_filename;
  using Parent::_pG;
  using Parent::_pIdNodeMap;
  using Parent::_pInvIdNodeMap;
  using Parent::_pWeightNodeMap;
  using Parent::_pWeightEdgeMap;
  using Parent::_nNodes;
  using Parent::_nEdges;

private:
  bool parseGraph(std::istream& in, int& lineNumber);
  bool parseNrNodes(std::istream& in, int& lineNumber);
  bool parseNode(std::istream& in, int& lineNumber);
  bool parseNrEdges(std::istream& in, int& lineNumber);
  bool parseEdge(std::istream& in, int& lineNumber);

public:
  DimacsParser(const std::string& filename);
  bool parse();
  const std::string& getName() const
  {
    return _name;
  }
  
protected:
  std::string _name;
};

template<typename GR>
inline DimacsParser<GR>::DimacsParser(const std::string& filename)
  : Parent(filename)
{
}

template<typename GR>
inline bool DimacsParser<GR>::parseNrNodes(std::istream& in, int& lineNumber)
{
  std::string line;

  if (safeGetline(in, line))
  {
    std::string text;
    std::stringstream ss(line);
    lineNumber++;
    ss >> text;

    if (text != "Vertices")
    {
      std::cerr << "Error at line " << lineNumber << ": expected 'Nodes'" << std::endl;
      return false;
    }

    ss >> _nNodes;
    _pG->reserveNode(_nNodes);
    return true;
  }
  else
  {
    std::cerr << "Premature end-of-file; expected node count" << std::endl;
    return false;
  }
}

template<typename GR>
inline bool DimacsParser<GR>::parseNrEdges(std::istream& in, int& lineNumber)
{
  std::string line;

  if (safeGetline(in, line))
  {
    std::string text;
    std::stringstream ss(line);
    lineNumber++;
    ss >> text;

    if (text != "Edges")
    {
      std::cerr << "Error at line " << lineNumber << ": expected 'Edges'" << std::endl;
      return false;
    }

    ss >> _nEdges;
    _pG->reserveEdge(_nEdges);
    
    return true;
  }
  else
  {
    std::cerr << "Premature end-of-file; expected edge count" << std::endl;
    return false;
  }
}

template<typename GR>
inline bool DimacsParser<GR>::parseGraph(std::istream& in, int& lineNumber)
{
  std::string line;
  
  // skip until "Name"
  while (safeGetline(in, line) && line.substr(0, 4) != "Name")
    lineNumber++;
  
  if (line.size() < 6)
  {
    std::cerr << "Error: missing 'Name'" << std::endl;
    return false;
  }
  _name = line.substr(5);

  // skip until "SECTION Finalsolution"
  while (safeGetline(in, line) && line != "SECTION Finalsolution")
    lineNumber++;

  if (line != "SECTION Finalsolution")
  {
    std::cerr << "Error: missing 'SECTION Finalsolution'" << std::endl;
    return false;
  }

  if (!parseNrNodes(in, lineNumber))
  {
    return false;
  }

  // add nodes
  for (int i = 0; i < _nNodes; i++)
  {
    if (!parseNode(in, lineNumber))
    {
      return false;
    }
  }

  // add edges
  if(!parseNrEdges(in, lineNumber))
  {
    return false;
  }

  for (int i = 0;i < _nEdges; i++)
  {
    if (!parseEdge(in, lineNumber))
    {
      return false;
    }
  }

  return true;
}
  
template<typename GR>
inline bool DimacsParser<GR>::parseNode(std::istream& in, int& lineNumber)
{
  std::string line;
  char buf[1024];

  if (safeGetline(in, line))
  {
    std::string text;
    std::stringstream ss(line);
    lineNumber++;
    ss >> text;

    if (text != "V")
    {
      std::cerr << "Error at line " << lineNumber << ": expected 'E'" << std::endl;
      return false;
    }

    int idU = -1;
    ss >> idU;

    if (idU < 1)
    {
      std::cerr << "Error at line " << lineNumber << ": expected node id in [1, "
                << "+inf]" << std::endl;
      return false;
    }
    
    Node u = _pG->addNode();
    _pWeightNodeMap->set(u, 0);
    
    snprintf(buf, 1024, "%d", idU);
    _pIdNodeMap->set(u, buf);
    (*_pInvIdNodeMap)[buf] = u;

    return true;
  }
  else
  {
    std::cerr << "Premature end-of-file; expected edge" << std::endl;
    return false;
  }
}

template<typename GR>
inline bool DimacsParser<GR>::parseEdge(std::istream& in, int& lineNumber)
{
  std::string line;

  if (safeGetline(in, line))
  {
    std::string text;
    std::stringstream ss(line);
    lineNumber++;
    ss >> text;

    if (text != "E")
    {
      std::cerr << "Error at line " << lineNumber << ": expected 'E'" << std::endl;
      return false;
    }

    std::string idU, idV ;
    ss >> idU >> idV;

    Node u = (*_pInvIdNodeMap)[idU];
    Node v = (*_pInvIdNodeMap)[idV];
    
    if (u == lemon::INVALID)
    {
      std::cerr << "Error at line " << lineNumber << ": invalid node id '" << idU << "'" << std::endl;
      return false;
    }
    if (v == lemon::INVALID)
    {
      std::cerr << "Error at line " << lineNumber << ": invalid node id '" << idV << "'" << std::endl;
      return false;
    }
    
    _pG->addEdge(u, v);

    return true;
  }
  else
  {
    std::cerr << "Premature end-of-file; expected edge" << std::endl;
    return false;
  }
}

template<typename GR>
inline bool DimacsParser<GR>::parse()
{
  if (!_pG)
    return false;

  std::ifstream in(_filename.c_str());
  if (!in.good())
  {
    std::cerr << "Error: could not open file "
              << _filename << " for reading" << std::endl;
    return false;
  }

  _pG->clear();

  int lineNumber = 0;
  return parseGraph(in, lineNumber);
}

} // namespace mwcs
} // namespace nina

#endif // DIMACSPARSER_H
