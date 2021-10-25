/*
 * Authors: Gustavo V. Barroso 
 * Created: 12/11/2020
 * Last modified: 16/11/2020
 *
 */

#include <cmath>
#include <random>
#include <iostream>

#include "RecombinationMap.h"
 
std::map<std::string, size_t> RecombinationMap::drawNumberOfCrossoversPerChr() const
{
  std::map<std::string, size_t> crossoverPerChr; //ret 
  
  std::mt19937 gen;
  std::array<int, 624> seedData;
  unsigned sem = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine re(sem);
  std::generate_n(seedData.data(), seedData.size(), std::ref(re));
  std::seed_seq seq(std::begin(seedData), std::end(seedData));
  gen.seed(seq);
  
  for(size_t i = 0; i < chrList_.size(); ++i)
  {
    double lambda = lambdas_.at(chrList_[i]);
    std::poisson_distribution<size_t> dist(lambda);
    
    size_t numCrossover = dist(gen);
    crossoverPerChr.emplace(chrList_[i], numCrossover);
    
    //std::cout << chrList_[i] << "\t" << lambda << "\t" << numCrossover << std::endl;
  }
  
  return crossoverPerChr;
}
 

void RecombinationMap::computeRates_()
{  
  for(size_t i = 0; i < chrList_.size(); ++i)
  {
    std::string chr = chrList_[i];
    
    std::vector<double> reparRates = values_.at(chr);
    auto itR = std::begin(ranges_.at(chr));
    
    for(size_t j = 0; j < reparRates.size(); ++j)
      reparRates[j] = -(std::log(1. - reparRates[j]));
    
    double lambda = 0.; // total rec. rate for focal chr
    for(size_t j = 0; j < reparRates.size(); ++j)
      lambda += reparRates[j] * (itR[j].second - itR[j].first);
    
    lambdas_.emplace(std::move(chr), lambda);
  }
}
