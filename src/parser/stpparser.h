/*
 * stpparser.h
 *
 *  Created on: 11-oct-2013
 *      Author: M. El-Kebir
 */

#ifndef STPPARSER_H
#define STPPARSER_H

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
class StpParser : public Parser<GR>
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
  bool parseHeader(std::istream& in, int& lineNumber);
  bool parseGraph(std::istream& in, int& lineNumber);
  bool parseNrNodes(std::istream& in, int& lineNumber);
  bool parseNrEdges(std::istream& in, int& lineNumber);
  bool parseNrTerminals(std::istream& in, int& lineNumber);
  bool parseEdge(std::istream& in, int& lineNumber);
  bool parseTerminal(std::istream& in, int& lineNumber);

public:
  StpParser(const std::string& filename);
  bool parse();
  const std::string& getName() const
  {
    return _name;
  }
  
protected:
  std::string _name;
};

template<typename GR>
inline StpParser<GR>::StpParser(const std::string& filename)
  : Parent(filename)
{
}

template<typename GR>
inline bool StpParser<GR>::parseHeader(std::istream& in, int& lineNumber)
{
  std::string line;
  if (safeGetline(in, line))
  {
    lineNumber++;
    return line.substr(0, 8) == "33D32945";
  }
  else
  {
    return false;
  }
}

template<typename GR>
inline bool StpParser<GR>::parseNrTerminals(std::istream& in, int& lineNumber)
{
  std::string line;

  if (safeGetline(in, line))
  {
    std::string text;
    std::stringstream ss(line);
    lineNumber++;
    ss >> text;

    if (text != "Terminals")
    {
      std::cerr << "Error at line " << lineNumber << ": expected 'Terminals'" << std::endl;
      return false;
    }

    int nTerminals = -1;
    ss >> nTerminals;
    if (nTerminals != _nNodes)
    {
      std::cerr << "Error at line " << lineNumber << ": terminal count must match node count" << std::endl;
      return false;
    }
    return true;
  }
  else
  {
    std::cerr << "Premature end-of-file; expected terminal count" << std::endl;
    return false;
  }
}

template<typename GR>
inline bool StpParser<GR>::parseNrNodes(std::istream& in, int& lineNumber)
{
  std::string line;

  if (safeGetline(in, line))
  {
    std::string text;
    std::stringstream ss(line);
    lineNumber++;
    ss >> text;

    if (text != "Nodes")
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
inline bool StpParser<GR>::parseNrEdges(std::istream& in, int& lineNumber)
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
inline bool StpParser<GR>::parseGraph(std::istream& in, int& lineNumber)
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

  // skip until "SECTION Graph"
  while (safeGetline(in, line) && line != "SECTION Graph")
    lineNumber++;

  if (line != "SECTION Graph")
  {
    std::cerr << "Error: missing 'SECTION Graph'" << std::endl;
    return false;
  }

  if (!parseNrNodes(in, lineNumber) || !parseNrEdges(in, lineNumber))
  {
    return false;
  }

  // add nodes
  for (int i = 0; i < _nNodes; i++)
  {
    _pG->addNode();
  }

  // add edges
  for (int i = 0;i < _nEdges; i++)
  {
    if (!parseEdge(in, lineNumber))
    {
      return false;
    }
  }

  // skip until "SECTION Terminals"
  while (safeGetline(in, line) && line != "SECTION Terminals")
    lineNumber++;

  if (line != "SECTION Terminals")
  {
    std::cerr << "Error: missing 'SECTION Terminals'" << std::endl;
    return false;
  }

  if (!parseNrTerminals(in, lineNumber))
  {
    return false;
  }

  for (int i = 0; i < _nNodes; i++)
  {
    if (!parseTerminal(in, lineNumber))
    {
      return false;
    }
  }

  return true;
}

template<typename GR>
inline bool StpParser<GR>::parseEdge(std::istream& in, int& lineNumber)
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

    int idU = -1, idV = -1;
    ss >> idU >> idV;

    if (!(0 < idU && idU <= _nNodes) || !(0 < idV && idV <= _nNodes))
    {
      std::cerr << "Error at line " << lineNumber << ": expected node id in [1, "
                << _nNodes << "]" << std::endl;
      return false;
    }

    Node u = _pG->nodeFromId(idU - 1);
    Node v = _pG->nodeFromId(idV - 1);
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
inline bool StpParser<GR>::parseTerminal(std::istream& in, int& lineNumber)
{
  std::string line;
  char buf[1024];

  if (safeGetline(in, line))
  {
    std::string text;
    std::stringstream ss(line);
    lineNumber++;
    ss >> text;

    if (text != "T")
    {
      std::cerr << "Error at line " << lineNumber << ": expected 'T'" << std::endl;
      return false;
    }

    int idU = -1;
    double weightU = -std::numeric_limits<double>::max();

    ss >> idU >> weightU;

    if (!(0 < idU && idU <= _nNodes))
    {
      std::cerr << "Error at line " << lineNumber << ": expected node id in [1, "
                << _nNodes << "]" << std::endl;
      return false;
    }

    if (weightU == -std::numeric_limits<double>::max())
    {
      std::cerr << "Error at line " << lineNumber << ": expected real-valued node weight"
                << std::endl;
      return false;
    }

    Node u = _pG->nodeFromId(idU - 1);
    _pWeightNodeMap->set(u, weightU);
    
    snprintf(buf, 1024, "%d", idU);
    _pIdNodeMap->set(u, buf);
    (*_pInvIdNodeMap)[buf] = u;

    return true;
  }
  else
  {
    std::cerr << "Premature end-of-file; expected terminal" << std::endl;
    return false;
  }
}

template<typename GR>
inline bool StpParser<GR>::parse()
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
  return parseHeader(in, lineNumber) && parseGraph(in, lineNumber);
}

} // namespace mwcs
} // namespace nina

#endif // STPPARSER_H
