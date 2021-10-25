/*
 * Authors: Gustavo V. Barroso 
 * Created: 14/11/2019
 * Last modified: 12/11/2020
 *
 */


#include <cmath>
#include <thread>
#include <random>

#include "Child.h"
#include "Tools.h"
#include "Global.h"


void Child::genPotentialSiblings(size_t n, const RecombinationMap& femMap, const RecombinationMap& maleMap) 
{
  chbs_.resize(n);
  
  size_t numThreads = std::min(n, NUM_AVAIL_THREADS);
  size_t batchSize = (n + numThreads - 1) / numThreads; // ceiling
  
  auto threadVec = std::vector<std::thread>(0);
  threadVec.reserve(numThreads);

  for(size_t i = 0; i < n; i += batchSize)
  {  
    size_t rest = n - i; // when we are near the end, just do the rest
    size_t focalBatchSize = std::min(batchSize, rest);
    
    // spawn a new thread for items [i..i + focalBatchSize)
    auto newThread = std::thread([&, i, focalBatchSize]() // capture only recMaps by ref
    {
      auto zygotes = genBatchOfSiblings_(focalBatchSize, femMap, maleMap);
      std::move(std::begin(zygotes), std::end(zygotes), std::begin(chbs_) + i); 
    });

    threadVec.emplace_back(std::move(newThread));
  }
  
  // wait for every thread to finish
  for(auto& t : threadVec)
    t.join();

}

void Child::addDeNovoMutations(double lambda)
{
  std::array<int, 624> seedData;
  unsigned sem = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine re(sem);
  std::generate_n(seedData.data(), seedData.size(), std::ref(re));
  std::seed_seq seq(std::begin(seedData), std::end(seedData));
  
  std::mt19937 gen; 
  gen.seed(seq);
  
  std::poisson_distribution<size_t> dist(lambda);
  
  std::vector<size_t> perSibling(0); // num de novo snps per sibling
  perSibling.reserve(chbs_.size());
  
  std::generate_n(std::back_inserter(perSibling), chbs_.size(),[&] { return dist(gen); });
  
  for(size_t i = 0; i < chbs_.size(); ++i)
    chbs_[i].setDeNovoCount(perSibling[i]);
}

std::vector<Zygote> Child::genBatchOfSiblings_(size_t numSibs, const RecombinationMap& femMap, const RecombinationMap& maleMap)
{
  auto genNewSibling = [&]()
  {
    Diploid d(std::make_pair(HaploidGenome(mother_->generateGamete(femMap)),
                             HaploidGenome(father_->generateGamete(maleMap))));
    
    return Zygote(d);
  };

  auto vec = std::vector<Zygote>(0);
  vec.reserve(numSibs); 
  
  std::generate_n(std::back_inserter(vec), numSibs, genNewSibling);

  return vec;
}
