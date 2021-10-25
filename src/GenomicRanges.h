/*
 * Authors: Gustavo V. Barroso 
 * Created: 22/10/2019
 * Last modified: 03/08/2020
 *
 */


#ifndef _GENOMICRANGES_
#define _GENOMICRANGES_

#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cmath>
#include <utility>

typedef std::vector<std::pair<unsigned int, unsigned int>> Windows;

// chr label to start & end coordinates 
typedef std::map<std::string, Windows> Ranges; 

// chr label to values
typedef std::map<std::string, std::vector<double>> Values;

class GenomicRanges {

protected:
  Ranges ranges_;
  Values values_;
  Values cdf_; 
  
  // to compare 2 bedgraphs faster
  std::vector<std::string> chrList_;

public:
  GenomicRanges(const Ranges& ranges, const Values& vals):
  ranges_(ranges),
  values_(vals),
  cdf_(),
  chrList_(0)
  { 
    genCdf_();    
  }

  GenomicRanges(const std::string& bedgraphPath):
  ranges_(),
  values_(),
  cdf_(),
  chrList_(0)
  {
    parseBedgraph_(bedgraphPath);
    genCdf_();
  }
  
  GenomicRanges():
  ranges_(),
  values_(),
  cdf_(),
  chrList_(0)
  { }
  
  GenomicRanges(GenomicRanges&& other):
  ranges_(),
  values_(),
  cdf_(),
  chrList_(0)
  {
    ranges_ = other.ranges_;
    values_ = other.values_;
    cdf_ = other.cdf_;
    chrList_ = other.chrList_;
  }
  
public:
  const Ranges& getRanges() const
  {
    return ranges_;
  }

  const Values& getValues() const
  {
    return values_;
  }
  
  const Values& getCdf() const
  {
    return cdf_;
  }
  
  const std::vector<std::string>& getChrList() const
  {
    return chrList_;
  }
  
  const size_t getNumChr() const 
  {
    return ranges_.size();
  }

  void setRanges(const Ranges& ranges)
  {
    ranges_ = ranges;
  }
  
  void setValues(const Values& values)
  {
    values_ = values;
  }
  
  void setChrList(const std::vector<std::string>& list)
  {
    chrList_ = list;
  }
  
  void writeToFile(const std::string& name);

  double fetchValue(const std::string& chr, unsigned int pos);

  const double fetchValue(const std::string& chr, unsigned int pos) const;

private:
  void genCdf_();

  void parseBedgraph_(const std::string& bedgraphPath);

};

#endif
