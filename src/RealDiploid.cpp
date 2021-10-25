/*
 * Authors: Gustavo V. Barroso 
 * Created: 21/10/2019
 * Last modified: 26/11/2020
 *
 */


#include <algorithm>
#include <utility>
#include <random>
#include <cmath>

#include "RealDiploid.h"
#include "Tools.h"
#include "Global.h"


Haplotype RealDiploid::generateGamete(const RecombinationMap& map) 
{
  // chr -> SNPs
  Haplotype gameteGenome; // ret
  
  // chr -> crossover positions
  std::map<std::string, std::vector<size_t>> recBps = drawBreakpoints_(map); 
  
  for(auto itC = std::begin(recBps); itC != std::end(recBps); ++itC) // for all chr
  { 
    std::string chr = itC->first;
    std::vector<size_t> recChr(0);
    
    auto p0begin = std::begin(haplotypePair_.first.getChromosome(chr));
    auto p0end = std::end(haplotypePair_.first.getChromosome(chr));
    
    auto p1begin = std::begin(haplotypePair_.second.getChromosome(chr));
    auto p1end = std::end(haplotypePair_.second.getChromosome(chr));
    
    size_t numCrossovers = itC->second.size();
    bool focalParent = boolean_(gen_); // which parent gives 1st part of haplotype
    //bool rec = boolean_(gen_); // if a chromatid with crossover is chosen
    bool rec = true; // WARNING for SLiM simulations (active above line real data)
    
    if(numCrossovers == 0 || !rec) // if no cross-over in chr
    {
      if(focalParent)
        recChr.insert(std::end(recChr), p1begin, p1end);
      else
        recChr.insert(std::end(recChr), p0begin, p0end);
    }
    else
    {
      // vectors of interators pointing crossover positions  
      std::vector<std::vector<size_t>::const_iterator> bpsItsZero(0);
      std::vector<std::vector<size_t>::const_iterator> bpsItsOne(0);
      
      bpsItsZero.reserve(numCrossovers);
      bpsItsOne.reserve(numCrossovers);
      
      for(size_t i = 0; i < numCrossovers; ++i) 
      {
        size_t pos = itC->second[i]; 
        
        bpsItsZero.emplace_back(std::lower_bound(p0begin, p0end, pos));
        bpsItsOne.emplace_back(std::lower_bound(p1begin, p1end, pos));
      }
        
      if(focalParent) // if copying starts with parent w/ index 1 in the pair
      {
        bpsItsOne.insert(std::begin(bpsItsOne), p1begin);
        
        if(numCrossovers % 2 == 0) // even num of crossovers
          bpsItsOne.push_back(p1end);
        else 
          bpsItsZero.push_back(p0end);
      }
      else
      {
        bpsItsZero.insert(std::begin(bpsItsZero), p0begin);   

        if(numCrossovers % 2 == 0) // even num of crossovers
          bpsItsZero.push_back(p0end);
        else 
          bpsItsOne.push_back(p1end);
      }
      
      // put the two parts of crossover segments in the recombining chr
      // two separate loops because bpsItsZero & bpsItsOne may differ in size
      // no problem because recChr will be sorted downstream
      
      for(size_t i = 0; i < bpsItsZero.size(); i+=2) 
        recChr.insert(std::end(recChr), bpsItsZero[i], bpsItsZero[i + 1]);
      
      for(size_t i = 0; i < bpsItsOne.size(); i+=2) 
        recChr.insert(std::end(recChr), bpsItsOne[i], bpsItsOne[i + 1]);
      
      // sorts and deletes duplicates (in case of overlapping regions near the crossover)
      std::sort(std::begin(recChr), std::end(recChr));
      recChr.erase(std::unique(std::begin(recChr), std::end(recChr)), std::end(recChr));
    }

    gameteGenome.emplace(std::make_pair(std::move(chr), std::move(recChr)));
  } 

  return gameteGenome;
}

std::map<std::string, std::vector<size_t>> RealDiploid::drawBreakpoints_(const RecombinationMap& map) 
{ 
  std::map<std::string, std::vector<size_t>> recPosPerChr; // ret
  std::map<std::string, size_t> crossPerChr = map.drawNumberOfCrossoversPerChr(); 
  
  auto itN = std::begin(crossPerChr); // number of crossovers
  auto itR = std::begin(map.getRanges()); // genomic windows
  auto itD = std::begin(map.getCdf()); // rec map CDF
  
  for(; itR != std::end(map.getRanges()); ++itR, ++itD, ++itN) // for each chr
  { 
    size_t numCrossovers = itN->second; // number of crossovers drawn for focal chr
    std::vector<size_t> recPos(0);
    recPos.reserve(numCrossovers);
    
    for(size_t i = 0; i < numCrossovers; ++i)
    {
      double sumRates = itD->second.back(); // last pos. of CDF
      double x = distValues_(gen_); // runif [0, 1]
      double y = x * sumRates;
    
      auto samp = std::upper_bound(std::begin(itD->second), std::end(itD->second), y);
      auto bin = samp - std::begin(itD->second); // genomic bin where crossover happens
    
      std::pair<size_t, size_t> window = itR->second[bin];
      std::uniform_int_distribution<size_t> distPos(window.first, window.second);

      recPos.emplace_back(distPos(gen_));
    }
    std::sort(std::begin(recPos), std::end(recPos));
    recPos.erase(std::unique(std::begin(recPos), std::end(recPos)), std::end(recPos));
    
    recPosPerChr.emplace(itR->first, std::move(recPos));
  }
  
  return recPosPerChr; 
}

