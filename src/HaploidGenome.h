/*
 * Authors: Gustavo V. Barroso 
 * Created: 18/10/2019
 * Last modified: 25/05/2020
 *
 */


#ifndef _HAPLOIDGENOME_H_
#define _HAPLOIDGENOME_H_

#include <vector>
#include <string>
#include <map>
#include <iterator>
#include <utility>
#include <iostream>

typedef std::map<std::string, std::vector<size_t>> Haplotype;

class HaploidGenome {

private:
  // we just store the positions of non-synonymous derived sites
  Haplotype chromosomes_; 

public:
  HaploidGenome(Haplotype& chromosomes):
  chromosomes_()
  { 
    chromosomes_.merge(chromosomes);
  }

  HaploidGenome(Haplotype&& chromosomes):
  chromosomes_()
  { 
    chromosomes_.merge(chromosomes);
  }
  
  HaploidGenome():
  chromosomes_()
  { }
  
public:
  const Haplotype& getGenome() const
  {
    return chromosomes_;
  }
  
  Haplotype& getGenome()
  {
    return chromosomes_;
  }
  
  const size_t getNumChromosomes() const 
  {
    return chromosomes_.size();
  }
  
  const std::vector<size_t>& getChromosome(const std::string& chrLabel) const 
  {
    return chromosomes_.at(chrLabel);
  }

  size_t fetchNumNonSynDerived() 
  {
    size_t number = 0;
    for(auto it = std::begin(chromosomes_); it != std::end(chromosomes_); ++it)
      number += it->second.size();

    return number;
  }
  
  void printSites();

};

#endif
