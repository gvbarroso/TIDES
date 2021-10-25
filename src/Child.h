/*
 * Authors: Gustavo V. Barroso 
 * Created: 14/11/2019
 * Last modified: 12/11/2020
 *
 */


#ifndef _CHILD_
#define _CHILD_

#include <utility>
#include <iostream>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "RealDiploid.h"
#include "Zygote.h"
#include "RecombinationMap.h"

class Child:
  public RealDiploid {

private:
  // the distribution of could-have-been siblings is invariant to model params,
  // so we simulate it only once and store / re-use it throughout the ABC simulations:
  std::vector<Zygote> chbs_; 
  
  std::shared_ptr<RealDiploid> mother_;
  std::shared_ptr<RealDiploid> father_;
  
public:
  Child(const std::pair<HaploidGenome, HaploidGenome>& hapPair,
        std::shared_ptr<RealDiploid> mom, std::shared_ptr<RealDiploid> dad,
        const std::string& name, const std::string& sex = "*"):
  RealDiploid(hapPair, name, sex),
  chbs_(0),
  mother_(mom),
  father_(dad)
  { }

public:
  const std::vector<Zygote>& getSiblings() const
  {
    return chbs_;
  }
  
  size_t getNumSiblings()
  {
    return chbs_.size();
  }

  std::shared_ptr<RealDiploid> getMother()
  {
    return mother_;
  }

  std::shared_ptr<RealDiploid> getFather()
  {
    return father_;
  }

  void genPotentialSiblings(size_t n, const RecombinationMap& femMap, const RecombinationMap& maleMap);
  
  void addDeNovoMutations(double lambda);

private:
  std::vector<Zygote> genBatchOfSiblings_(size_t number, const RecombinationMap& femMap, const RecombinationMap& maleMap);
  
};

#endif
