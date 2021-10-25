/*
 * Authors: Gustavo V. Barroso 
 * Created: 21/10/2019
 * Last modified: 23/11/2020
 *
 */


#ifndef _REALDIPLOID_
#define _REALDIPLOID_

#include <random>
#include <array>
#include <chrono>

#include "Diploid.h"
#include "RecombinationMap.h"
#include "HaploidGenome.h"


class RealDiploid:
  public Diploid {

private:
  std::mt19937 gen_;
  std::uniform_int_distribution<int> boolean_;
  std::uniform_real_distribution<double> distValues_;
  std::string name_;
  
public:
  RealDiploid(const std::pair<HaploidGenome, HaploidGenome>& hapPair,
              const std::string& name, const std::string& sex = "*"):
  Diploid(hapPair, sex),
  gen_(),
  boolean_(0, 1),
  distValues_(0., 1.),
  name_(name)
  { 
    std::array<int, 624> seedData;
    unsigned sem = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re(sem);
    std::generate_n(seedData.data(), seedData.size(), std::ref(re));
    std::seed_seq seq(std::begin(seedData), std::end(seedData));
    gen_.seed(seq);
  }

  RealDiploid():
  Diploid(),
  gen_(),
  boolean_(),
  distValues_(),
  name_("")
  { }
  
public:
  const std::string& getName() const
  {
    return name_;
  }
  
  void setSeed(size_t seed)
  {
    gen_.seed(seed);
  }
  
  void resetSeed()
  {
    std::array<int, 624> seedData;
    unsigned sem = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re(sem);
    std::generate_n(seedData.data(), seedData.size(), std::ref(re));
    std::seed_seq seq(std::begin(seedData), std::end(seedData));
    gen_.seed(seq);
  }
  
  // generates a single gamete via recombination of haploidPair_
  Haplotype generateGamete(const RecombinationMap& map);

private:
  std::map<std::string, std::vector<size_t>> drawBreakpoints_(const RecombinationMap& m);

};

#endif
