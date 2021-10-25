/*
 * Authors: Gustavo V. Barroso 
 * Created: 22/10/2019
 * Last modified: 25/05/2020
 *
 */


#ifndef _DIPLOID_
#define _DIPLOID_

#include <vector>
#include <string>
#include <utility>

#include "HaploidGenome.h"

class Diploid {

protected:
  std::string sex_; // f: female | m: male
  std::pair<HaploidGenome, HaploidGenome> haplotypePair_; // 0: mother | 1: father
  size_t numHetSites_; // stores for fast access 
  size_t numHomoSites_; // stores for fast access 

public:
  Diploid(const std::pair<HaploidGenome, HaploidGenome>& p, const std::string& sex = "*"):
  sex_(sex),
  haplotypePair_(p),
  numHetSites_(0),
  numHomoSites_(0)
  { 
    countSites();
  }

  Diploid():
  sex_("*"),
  haplotypePair_(),
  numHetSites_(0),
  numHomoSites_(0)
  { }

public:
  const std::string& getSex() 
  {
    return sex_;
  }
  
  const std::pair<HaploidGenome, HaploidGenome>& getHaplotypePair() const
  {
    return haplotypePair_;
  }
  
  std::pair<HaploidGenome, HaploidGenome>& getHaplotypePair()
  {
    return haplotypePair_;
  }
  
  size_t getNumHetSites() const 
  {
    return numHetSites_;
  }

  size_t getNumHomoSites() const 
  {
    return numHomoSites_;
  }
  
  void countSites(); // computes numHetSites_ and numHomoSites_

};

#endif
