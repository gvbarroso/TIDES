/*
 * Authors: Gustavo V. Barroso 
 * Created: 24/10/2019
 * Last modified: 24/07/2020
 *
 */


#ifndef _SIMULATOR_
#define _SIMULATOR_

#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <map>
#include <utility>
#include <memory>

#include "HashTable.h"
#include "RealDiploid.h"
#include "Child.h"
#include "Tools.h"



class Simulator {

private:
  std::map<std::string, std::vector<double>> newPrior_; // name->values 
  std::map<std::string, std::vector<double>> simParams_; // name->values 
  std::map<std::string, std::vector<double>> simStats_; /// name->values 
  std::map<std::string, double> obsStats_; // name->value
  
  std::pair<double, double> sUnif_;
  std::pair<double, double> hUnif_;
  
  std::pair<double, double> parentMeans_; // 1st: mean het; 2nd: mean homo in parents
    
  bool isPilot_; // tells whether we are in pilot mode
  
public:
  Simulator(const std::vector<double>& selInterval,
            const std::vector<double>& domInterval):
  newPrior_(),
  simParams_(),
  simStats_(),
  obsStats_(),
  sUnif_(std::make_pair(selInterval[0], selInterval[1])),
  hUnif_(std::make_pair(domInterval[0], domInterval[1])),
  parentMeans_(),
  isPilot_()
  { }
  
public:
  const std::map<std::string, std::vector<double>>& getAllSimParams() const 
  {
    return simParams_;
  }

  const std::vector<double>& getSimParamValues(const std::string& param) const 
  {
    return simParams_.at(param);
  }  
  
  void setPilot(bool mode)
  {
    isPilot_ = mode;
  }
  
  void setNewPrior(const std::string& paramName, const std::vector<double>& v)
  {
    newPrior_.at(paramName) = v;
  }
  
  void abcReject(size_t n, size_t f);  
      
  void init(size_t n);
  
  void reset(size_t n);
  
  void pilot(const std::vector<std::shared_ptr<Child>>& kids, size_t n, size_t f);
  
  void simulate(const std::vector<std::shared_ptr<Child>>& kids, size_t n);

  void writeSimParamsToFiles(size_t n);
  
  void writeStatsToFiles(size_t n);
  
  void computeObsStats(const std::vector<std::shared_ptr<Child>>& kids);

  void writeChbsSiteDistToFile(const std::vector<std::shared_ptr<Child>>& kids);
  
  void writeTriosSiteDistToFile(const std::vector<std::shared_ptr<Child>>& kids);

private:
  void evalParams_(const std::vector<std::shared_ptr<Child>>& kids, size_t s, size_t b);
  
  double drawH_();
  
  double drawS_();
};

#endif

