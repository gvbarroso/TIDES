/*
 * Authors: Gustavo V. Barroso 
 * Created: 12/11/2020
 * Last modified: 16/11/2020
 *
 */


#ifndef _RECOMBINATIONMAP_
#define _RECOMBINATIONMAP_

#include <chrono>

#include "GenomicRanges.h"

class RecombinationMap:
  public GenomicRanges {

private:
  std::string sex_; // "f" => female (mothers); "m" => male (fathers)
  std::map<std::string, double> lambdas_; // total rec. rate per chr (for Poisson draw)

public:
  RecombinationMap(const Ranges& ranges, const Values& vals, const std::string& sex = "*"):
  GenomicRanges(ranges, vals),
  sex_(sex),
  lambdas_()
  { 
    computeRates_();
  }

  RecombinationMap(const std::string& bedgraphPath, const std::string& sex = "*"):
  GenomicRanges(bedgraphPath),
  sex_(sex),
  lambdas_()
  { 
    computeRates_();
  }
  
  RecombinationMap():
  GenomicRanges(),
  sex_(),
  lambdas_()
  { }
  
  RecombinationMap(RecombinationMap&& other):
  GenomicRanges(std::move(other)),
  sex_(),
  lambdas_()
  {
    sex_ = other.sex_;
    lambdas_ = other.lambdas_;
  }
  
public:
  const std::string& getSex() const
  {
    return sex_;
  }
  
  const std::map<std::string, double>& getRates() const
  {
    return lambdas_;
  }
  
  std::map<std::string, size_t> drawNumberOfCrossoversPerChr() const;
  
private:  
  void computeRates_();
  
};
  
#endif
