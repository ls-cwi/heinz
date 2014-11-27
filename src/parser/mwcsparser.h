/*
 * mwcsparser.h
 *
 *  Created on: 7-aug-2012
 *     Authors: C.I. Bucur, M. El-Kebir
 */

#ifndef MWCSPARSER_H
#define MWCSPARSER_H

#include <assert.h>
#include <string.h>
#include <map>
#include <string>
#include <lemon/core.h>
#include "parser.h"
#include "utils.h"

namespace nina {
namespace mwcs {

template<typename GR>
class MwcsParser : public Parser<GR>
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
  const std::string _filenameEdges;
  bool parseNodes();
  bool parseEdges();

public:
  MwcsParser(const std::string &fNodes, const std::string &fEdges);
  bool parse();
};

template<typename GR>
inline MwcsParser<GR>::MwcsParser(const std::string &fNodes,
                                  const std::string &fEdges)
  : Parent(fNodes)
  , _filenameEdges(fEdges)
{
}

template<typename GR>
inline bool MwcsParser<GR>::parse()
{
  if (!_pG)
    return false;

  _pG->clear();

  return parseNodes() && parseEdges();
}

template<typename GR>
inline bool MwcsParser<GR>::parseNodes()
{
  assert(_pG);

  std::ifstream inFile(_filename.c_str());
  if (!inFile.good())
  {
    std::cerr << "Error: could not open file "
              << _filename << " for reading" << std::endl;
    return false;
  }

  std::string line;

  int lineNumber = 0;

  while (safeGetline(inFile, line))
  {
    ++lineNumber;
    if ((!line.empty() && line[0] == '#') || line.empty())
      continue;

    std::stringstream lineStream(line);

    std::string label;
    double score;

    lineStream >> label;
    if (!lineStream.good())
    {
      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        std::cout << "File error: label is wrong at line " << lineNumber << std::endl;
      }
      return false;
    }
    lineStream >> score;
    if (lineStream.bad() || lineStream.fail())
    {
      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        std::cout << "File error: score is wrong at line " << lineNumber << std::endl;
      }
      return false;
    }

    if (!lineStream.eof())
    {
      if (g_verbosity >= VERBOSE_DEBUG)
      {
        std::cout << "Warning: trailing characters at line " << lineNumber << std::endl;
      }
    }

    if (_pInvIdNodeMap->find(label) != _pInvIdNodeMap->end())
    {
      if (g_verbosity >= VERBOSE_DEBUG)
      {
        std::cout << "Warning: duplicate node with label " << label << " at line " << lineNumber << std::endl;
      }
      continue;
    }

    Node x = _pG->addNode();

    if (_pIdNodeMap) _pIdNodeMap->set(x, label);
    if (_pWeightNodeMap) _pWeightNodeMap->set(x, score);
    (*_pInvIdNodeMap)[label] = x;

    _nNodes++;
  }

  return true;
}

template<typename GR>
inline bool MwcsParser<GR>::parseEdges()
{
  assert(_pG);

  std::ifstream inFile(_filenameEdges.c_str());
  if (!inFile.good())
  {
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cerr << "Error: could not open file "
              << _filenameEdges << " for reading" << std::endl;
    }
    return false;
  }

  std::string line;

  int lineNumber = 0;

  while (safeGetline(inFile, line))
  {
    ++lineNumber;
    if ((!line.empty() && line[0] == '#') || line.empty())
      continue;

    std::stringstream lineStream(line);

    std::string label1, label2, trail;

    lineStream >> label1;
    if (!lineStream.good())
    {
      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        std::cout << "File error: first label is wrong at line " << lineNumber << std::endl;
      }
      return false;
    }
    lineStream >> label2;
    if (lineStream.bad() || lineStream.fail())
    {
      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        std::cout << "File error: second label is wrong at line " << lineNumber << std::endl;
      }
      return false;
    }

    lineStream >> trail;
    if (!trail.empty())
    {
      if (g_verbosity >= VERBOSE_DEBUG)
      {
        std::cout << "Warning: trailing characters at line "
                << lineNumber << std::endl;
      }
    }

    if (_pInvIdNodeMap->find(label1) == _pInvIdNodeMap->end())
    {
      if (g_verbosity >= VERBOSE_DEBUG)
      {
        std::cout << "Warning: node with label " << label1
                << " at line " << lineNumber
                << " not specified. The edge is skipped." << std::endl;
      }
      continue;
    }

    if (_pInvIdNodeMap->find(label2) == _pInvIdNodeMap->end())
    {
      if (g_verbosity >= VERBOSE_DEBUG)
      {
        std::cout << "Warning: node with label " << label2
                << " at line " << lineNumber
                << " not specified. The edge is skipped." << std::endl;
      }
      continue;
    }

    Node node1 = (*_pInvIdNodeMap)[label1];
    Node node2 = (*_pInvIdNodeMap)[label2];

    if (node1 == node2)
    {
      if (g_verbosity >= VERBOSE_DEBUG)
      {
        std::cout << "Warning: node with label " << label1
                << " at line " << lineNumber
                << " has a self-loop. The edge is skipped." << std::endl;
      }
      continue;
    }

    _pG->addEdge((*_pInvIdNodeMap)[label1], (*_pInvIdNodeMap)[label2]);
    _nEdges++;
  }

  return true;
}

} // namespace mwcs
} // namespace nina

#endif // MWCSPARSER_H
