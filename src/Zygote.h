/*
 * Authors: Gustavo V. Barroso 
 * Created: 14/11/2019
 * Last modified: 25/05/2020
 *
 */


#ifndef _ZYGOTE_
#define _ZYGOTE_

#include <utility>
#include <algorithm>


class Zygote {
private:
  size_t numHomoSites_;    
  size_t numHetSites_;
  size_t deNovoCount_; // num de novo mutations 

  // which of the haplotypes (in each of the parents) is copied first:
  std::vector<bool> momHap_; // 1 per chr
  std::vector<bool> dadHap_; // 1 per chr
  
  // DSB coordinates:
  std::vector<size_t> momBps_; // 1 per chr
  std::vector<size_t> dadBps_; // 1 per chr
  
public:
  Zygote():
  numHomoSites_(0),
  numHetSites_(0),
  deNovoCount_(0),
  momHap_(0),
  dadHap_(0),
  momBps_(0),
  dadBps_(0)
  { }
  
  Zygote(const Diploid& d):
  numHomoSites_(0),
  numHetSites_(0),
  deNovoCount_(0),
  momHap_(0),
  dadHap_(0),
  momBps_(0),
  dadBps_(0)
  {
    numHomoSites_ = d.getNumHomoSites();
    numHetSites_ = d.getNumHetSites();
  }
  
  Zygote(Zygote&& other):
  numHomoSites_(0),
  numHetSites_(0),
  deNovoCount_(0),
  momHap_(0),
  dadHap_(0),
  momBps_(0),
  dadBps_(0)
  { 
    numHomoSites_ = other.numHomoSites_;
    numHetSites_ = other.numHetSites_;
    deNovoCount_ = other.deNovoCount_;
    std::swap(momHap_, other.momHap_);
    std::swap(dadHap_, other.dadHap_);
    std::swap(momBps_, other.momBps_);
    std::swap(dadBps_, other.dadBps_);
    
    other.numHomoSites_ = 0;
    other.numHetSites_ = 0;
    other.deNovoCount_ = 0;
  }
  
  Zygote& operator=(const Zygote& other)
  {
    if(this == &other)
      return *this;
    else
    {
      numHomoSites_ = other.numHomoSites_;
      numHetSites_ = other.numHetSites_;
      deNovoCount_ = other.deNovoCount_;
      momHap_ = other.momHap_;
      dadHap_ = other.dadHap_;
      momBps_ = other.momBps_;
      dadBps_ = other.dadBps_;
    }
    
    return *this;
  }
  
  Zygote& operator=(Zygote& other)
  {
    if(this == &other)
      return *this;
    else
    {
      numHomoSites_ = other.numHomoSites_;
      numHetSites_ = other.numHetSites_;
      deNovoCount_ = other.deNovoCount_;
      momHap_ = other.momHap_;
      dadHap_ = other.dadHap_;
      momBps_ = other.momBps_;
      dadBps_ = other.dadBps_;
    }
    
    return *this;
  }
  
public:
  size_t getNumHetSites() const 
  {
    return numHetSites_ + deNovoCount_;
  }

  size_t getNumHomoSites() const 
  {
    return numHomoSites_;
  }
  
  const std::vector<bool>& getMomHap() const
  {
    return momHap_;
  }
  
  const std::vector<bool>& getDadHap() const
  {
    return dadHap_;
  }
  
  const std::vector<size_t>& getMomBps() const
  {
    return momBps_;
  }
  
  const std::vector<size_t>& getDadBps() const
  {
    return dadBps_;
  }
  
  size_t getDeNovoCount() const 
  {
    return deNovoCount_;
  }
  
  void setDeNovoCount(size_t count)
  {
    deNovoCount_ = count;
  }
  
  void setNumHetSites(size_t count)
  {
    numHetSites_ = count;
  }
  
  void setNumHomoSites(size_t count)
  {
    numHomoSites_ = count;
  }
  
};

#endif
