/*
 * Authors: Gustavo V. Barroso 
 * Created: 22/10/2019
 * Last modified: 23/11/2020
 *
 */

#include <iterator>
#include <algorithm>

#include "Diploid.h"


void Diploid::countSites() {
    
  size_t sumHomoPos = 0; // sum over chromosomes
  size_t sumHetPos = 0; // sum over chromosomes

  auto gen1 = std::begin(haplotypePair_.first.getGenome()); // "maternal"
  auto gen2 = std::begin(haplotypePair_.second.getGenome()); // "paternal"

  // for each chr (gen1 and gen2 share the same chr's)
  for(; gen1 != std::end(haplotypePair_.first.getGenome()); ++gen1, ++gen2) 
  {
    auto startMom = std::begin(gen1->second);
    auto endMom = std::end(gen1->second);
    auto startDad = std::begin(gen2->second);
    auto endDad = std::end(gen2->second);

    std::vector<size_t> overlap(0);
    overlap.reserve(std::min(gen1->second.size(), gen2->second.size()));
  
    std::set_intersection(startMom, endMom, startDad, endDad,std::back_inserter(overlap));

    sumHomoPos += overlap.size();
    sumHetPos += gen1->second.size() + gen2->second.size() - 2 * overlap.size();  
  }

  numHomoSites_ = sumHomoPos;
  numHetSites_ = sumHetPos;
}
