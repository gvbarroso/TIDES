/*
 * Authors: Gustavo V. Barroso 
 * Created: 28/10/2019
 * Last modified: 21/12/2020
 *
 */

#include <array>
#include <cmath>
#include <thread>
#include <chrono>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "BppTextTools.h"
#include "Simulator.h"
#include "Global.h"
#include "Tools.h"
#include "Utils.h"

using namespace boost::accumulators;

void Simulator::init(size_t numParamDraws)
{
  // initialize map of sim params
  std::vector<double> tmp(numParamDraws);
  simParams_.emplace(std::make_pair("h", tmp));  
  simParams_.emplace(std::make_pair("s", tmp));  
    
  // initialize map of summary stats
  simStats_.emplace(std::make_pair("diff_het", tmp));
  simStats_.emplace(std::make_pair("diff_homo", tmp));
}

void Simulator::reset(size_t numParamDraws)
{
  simParams_.clear();
  simStats_.clear();
  
  init(numParamDraws);
}

void Simulator::abcReject(size_t numPilot, size_t numAccepted)
{
  std::vector<double> obs(0); // vector of observed sum stats values
  obs.reserve(obsStats_.size());
  for(auto it = std::begin(obsStats_); it != std::end(obsStats_); ++it)
    obs.emplace_back(it->second);
  
  // standardizing sum stats
  size_t obsIndex = 0; // helps standardizing obs sum stat
  for(auto it = std::begin(simStats_); it != std::end(simStats_); ++it)
  {
    std::pair<double, double> p = Utils::mad(it->second);
    double median = p.first;
    double mad = p.second;
    
    for(size_t i = 0; i < it->second.size(); ++i)
      it->second[i] = (it->second[i] - median) / mad;

    obs[obsIndex] = (obs[obsIndex] - median) / mad;    
    ++obsIndex;
  } 
  
  std::vector<double> euclidean(0); // euclidean dist for every sim. i
  euclidean.reserve(numPilot);
  
  for(size_t i = 0; i < numPilot; ++i) // for each pilot simulation
  {
    std::vector<double> sim(0);
    sim.reserve(obs.size());
    
    for(auto it = std::begin(simStats_); it != std::end(simStats_); ++it)
      sim.emplace_back(it->second[i]); // standardized sum stats for sim. i
    
    double d = 0.;
    for(size_t j = 0; j < sim.size(); ++j)
      d += std::pow((sim[j] - obs[j]), 2.);
    
    euclidean[i] = std::sqrt(d); 
    //std::cout << euclidean[i] << std::endl; // NOTE DEBUG
  }

  // sorts simulation INDICES according to their euclidean distances
  std::vector<size_t> simIndices(numPilot);
  std::iota(std::begin(simIndices), std::end(simIndices), 0);
  
  std::sort(std::begin(simIndices), std::end(simIndices),
            [&](std::size_t i, std::size_t j){return euclidean[i] < euclidean[j];});
    
  std::cout << "\n\nKept " << numAccepted << " simulations ";
  std::cout << "after rejection algorithm in pilot run.\n\n";
  
  auto start = std::begin(simIndices);
  auto finish = std::begin(simIndices) + numAccepted;
  
  std::vector<size_t> acceptedSims(start, finish); // indices of accepted sims
  for(auto it = std::begin(simParams_); it != std::end(simParams_); ++it)
  { 
    std::vector<double> acceptedVals(numAccepted);
    
    for(size_t i = 0; i < numAccepted; ++i)
      acceptedVals[i] = it->second[acceptedSims[i]];
        
    // for the "step function" in drawParameters
    std::sort(std::begin(acceptedVals), std::end(acceptedVals)); 
    newPrior_.emplace(std::make_pair(it->first, std::move(acceptedVals)));
  }
}

void Simulator::pilot(const std::vector<std::shared_ptr<Child>>& kids, size_t numPilot, size_t numAccepted)
{
  simulate(kids, numPilot); // pilot sims  
  abcReject(numPilot, numAccepted); // fills in newPrior_
}
  
void Simulator::writeSimParamsToFiles(size_t numParamDraws)
{
  boost::iostreams::filtering_ostream paramStream; 
  
  std::ofstream paramFile;
  paramFile.open("sim_params_dist.gz", std::ios_base::out | std::ios_base::binary);
  
  paramStream.push(boost::iostreams::gzip_compressor());
  paramStream.push(paramFile);  
  
  for(auto it = std::begin(simParams_); it != std::end(simParams_); ++it)
    paramStream << it->first << "\t";
  
  paramStream << std::endl;
  
  // all the iterators for simParams_
  std::vector<std::map<std::string, std::vector<double>>::const_iterator> parIts(0);
  parIts.reserve(simParams_.size());
  
  for(size_t i = 0; i < simParams_.size(); ++i)
  {
    auto it = std::begin(simParams_);
    std::advance(it, i);
    parIts.emplace_back(it);
  }
  
  // writes in "matrix format" (all params values and same file)
  // columns = param; rows = simulation index
  for(size_t i = 0; i < numParamDraws; ++i) 
  {
    for(auto it = std::begin(parIts); it != std::end(parIts); ++it)
      paramStream << (*it)->second[i] << "\t";
    
    paramStream << std::endl;
  }

  boost::iostreams::close(paramStream);
}  

void Simulator::writeStatsToFiles(size_t numParamDraws)
{  
  boost::iostreams::filtering_ostream simStatsStream; 
  
  std::ofstream sumStatFile;
  sumStatFile.open("sum_stats_dist.gz", std::ios_base::out | std::ios_base::binary);
  
  simStatsStream.push(boost::iostreams::gzip_compressor());
  simStatsStream.push(sumStatFile);  
  
  for(auto it = std::begin(simStats_); it != std::end(simStats_); ++it) // header
    simStatsStream << it->first << "\t";
  
  simStatsStream << std::endl;
  
  // all the iterators for simStats_
  std::vector<std::map<std::string, std::vector<double>>::const_iterator> pIts(0);
  pIts.reserve(simStats_.size());
  
  for(size_t i = 0; i < simStats_.size(); ++i)
  {
    auto it = std::begin(simStats_);
    std::advance(it, i);
    pIts.emplace_back(it);
  }
  
  // writes in "matrix format" (all sum stats values and same file)
  // columns = sum stat; rows = simulation index
  for(size_t i = 0; i < numParamDraws; ++i) 
  {
    for(auto it = std::begin(pIts); it != std::end(pIts); ++it)
      simStatsStream << (*it)->second[i] << "\t";
    
    simStatsStream << std::endl;
  }

  boost::iostreams::close(simStatsStream);
  
  // observed summary stats:
  std::ofstream obsStatFile;
  obsStatFile.open("obs_sum_stats.txt");
  
  for(auto it = std::begin(obsStats_); it != std::end(obsStats_); ++it) // header
    obsStatFile << it->first << "\t";
  
  obsStatFile << std::endl;
  
  for(auto it = std::begin(obsStats_); it != std::end(obsStats_); ++it) // header
    obsStatFile << it->second << "\t";
  
  obsStatFile << std::endl;
  obsStatFile.close();
}

void Simulator::computeObsStats(const std::vector<std::shared_ptr<Child>>& kids)
{
  size_t n = kids.size();
  
  // for summary statistics in accumulators below
  std::vector<double> kidsHetSites(n);
  std::vector<double> kidsHomoSites(n);
  
  // (each parent is represented proportional to its number of children)
  std::vector<double> parentsHetSites(n);
  std::vector<double> parentsHomoSites(n);
  
  for(size_t i = 0; i < n; ++i)
  {
    kidsHetSites[i] = kids[i]->getNumHetSites();
    kidsHomoSites[i] = kids[i]->getNumHomoSites();
    
    parentsHetSites[i] = (kids[i]->getMother()->getNumHetSites() +
                          kids[i]->getFather()->getNumHetSites()) / 2.;
    
    parentsHomoSites[i] = (kids[i]->getMother()->getNumHomoSites() +
                           kids[i]->getFather()->getNumHomoSites()) / 2.;
  }
  
  accumulator_set<double, stats<tag::mean>> accHetParents;
  std::for_each(std::begin(parentsHetSites), std::end(parentsHetSites),
                boost::bind<void>(boost::ref(accHetParents), _1));
  
  accumulator_set<double, stats<tag::mean>> accHomoParents;
  std::for_each(std::begin(parentsHomoSites), std::end(parentsHomoSites),
                boost::bind<void>(boost::ref(accHomoParents), _1));
  
  // stores for fast access
  parentMeans_.first = mean(accHetParents);
  parentMeans_.second = mean(accHomoParents);
  
  accumulator_set<double, stats<tag::mean, tag::variance>> accHomo;
  std::for_each(std::begin(kidsHomoSites), std::end(kidsHomoSites),
                boost::bind<void>(boost::ref(accHomo), _1));

  accumulator_set<double, stats<tag::mean, tag::variance>> accHet;
  std::for_each(std::begin(kidsHetSites), std::end(kidsHetSites),
                boost::bind<void>(boost::ref(accHet), _1));
  
  obsStats_.emplace(std::make_pair("diff_het",
(mean(accHet) - parentMeans_.first) / parentMeans_.first));
  obsStats_.emplace(std::make_pair("diff_homo",
(mean(accHomo) - parentMeans_.second) / parentMeans_.second));
}

void Simulator::simulate(const std::vector<std::shared_ptr<Child>>& kids, size_t numDraws)
{ 
  reset(numDraws); // resets simParams_ and simStats_ to size numDraws
  
  size_t numThreads = std::min(numDraws, NUM_AVAIL_THREADS);
  size_t batchSize = (numDraws + numThreads - 1) / numThreads;
  
  auto threadVec = std::vector<std::thread>(0);
  threadVec.reserve(numThreads);
  
  for(size_t i = 0; i < numDraws; i += batchSize)
  {  
    size_t rest = numDraws - i; // when we are near the end, just do the rest
    size_t focalBatchSize = std::min(batchSize, rest);

    // spawn a new thread for items [i..i + focalBatchSize)
    auto newThread = std::thread([&, i, focalBatchSize]() // capture only kids by ref
    {
      evalParams_(kids, i, focalBatchSize);  
    });

    threadVec.emplace_back(std::move(newThread));
  }
  
  for(auto& t : threadVec)
    t.join();
}

void Simulator::evalParams_(const std::vector<std::shared_ptr<Child>>& children, size_t startIndex, size_t batchSize) 
{  
  std::vector<std::shared_ptr<Child>> kids(children); // copy shared ptrs
  std::vector<double> chbsHomo(kids.size());
  std::vector<double> chbsHet(kids.size());
  std::vector<double> chbsFit(kids.front()->getNumSiblings());
    
  std::uniform_real_distribution<double> runif(0., 1.); 
  std::array<int, 624> seedData;
  unsigned sem = std::chrono::system_clock::now().time_since_epoch().count() + startIndex;
  std::default_random_engine re(sem);
  std::generate_n(seedData.data(), seedData.size(), std::ref(re));
  std::seed_seq seq(std::begin(seedData), std::end(seedData));
  
  std::mt19937 gen; 
  gen.seed(seq);

  Tools pb; 
  for(size_t i = startIndex; i < (startIndex + batchSize); ++i)
  { 
    double h = drawH_();
    double s = drawS_();
    
    for(size_t j = 0; j < kids.size(); ++j) 
    {
      // computes zygote fitnesses based on param draw  
      for(size_t k = 0; k < chbsFit.size(); ++k)
      { 
        double fitHomo = std::exp(-(kids[j]->getSiblings()[k].getNumHomoSites() * s));
        double fitHet = std::exp(-(kids[j]->getSiblings()[k].getNumHetSites() * h * s));
        
        double fitTotal = fitHomo * fitHet;
        // it is possible that fitTotal > 1 for some parameter combinations (eg h < 0).
        // since we are interested in absolute fitness (viability selection), we limit:
        chbsFit[k] = std::min(1., fitTotal);
      }

      for(size_t k = 1; k < chbsFit.size(); ++k)
        chbsFit[k] += chbsFit[k - 1]; // accumulates
      
      double sumFit = chbsFit.back();
      double x = runif(gen); // runif [0, 1]
      double y = x * sumFit;
      
      auto samp = std::upper_bound(std::begin(chbsFit), std::end(chbsFit), y);
      auto idx = samp - std::begin(chbsFit);
            
      chbsHomo[j] = kids[j]->getSiblings()[idx].getNumHomoSites();
      chbsHet[j] = kids[j]->getSiblings()[idx].getNumHetSites();
    } 
    
    accumulator_set<double, stats<tag::mean, tag::variance>> accHomo;
    std::for_each(std::begin(chbsHomo), std::end(chbsHomo),
                  boost::bind<void>(boost::ref(accHomo), _1));

    accumulator_set<double, stats<tag::mean, tag::variance>> accHet;
    std::for_each(std::begin(chbsHet), std::end(chbsHet),
                  boost::bind<void>(boost::ref(accHet), _1));

    simStats_.at("diff_het")[i] = (mean(accHet) - parentMeans_.first)/parentMeans_.first;
    simStats_.at("diff_homo")[i] = (mean(accHomo) - parentMeans_.second)/parentMeans_.second;
    
    // stores param draw itself
    simParams_.at("h")[i] = h;
    simParams_.at("s")[i] = s;
    
    if(startIndex == 0) // hacky: only prints pb for the 1st thread
      pb.print_pb((i + 1) / static_cast<double>(startIndex + batchSize));
  } 
}

double Simulator::drawH_()
{ 
  unsigned sem = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine re(sem);
  std::mt19937 gen;
  gen.seed(re());
  
  if(isPilot_) // if pilot run mode
  {
    std::uniform_real_distribution<double> unifPrior(hUnif_.first, hUnif_.second);   
    return unifPrior(gen);
  }
  else
  {
    size_t numAccepted = newPrior_.at("h").size();
    
    // to draw a random value in the sorted vector of pilotPosteriors
    std::uniform_int_distribution<size_t> distVals(0, numAccepted - 1);
    size_t pos = distVals(gen);
    
    if(pos == 0) // if first pos in the vector there's no left neighbor
    {
      std::uniform_real_distribution<double> draw(newPrior_.at("h")[pos],
                                                  newPrior_.at("h")[pos + 1]);
      return draw(gen);
    }
    else if(pos == numAccepted - 1) // if last pos in the vector there's no right neighbor
    {
      std::uniform_real_distribution<double> draw(newPrior_.at("h")[pos],
                                                  newPrior_.at("h")[pos - 1]);
      return draw(gen);
    }
    else
    {
      std::uniform_real_distribution<double> draw(newPrior_.at("h")[pos - 1],
                                                  newPrior_.at("h")[pos + 1]);
      return draw(gen);
    }    
  }
}

double Simulator::drawS_()
{ 
  unsigned sem = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine re(sem);
  std::mt19937 gen;
  gen.seed(re());
    
  if(isPilot_) // if pilot run mode
  {
    std::uniform_real_distribution<double> unifPrior(sUnif_.first, sUnif_.second);  
    double exponent = unifPrior(gen);
    return std::pow(10., -exponent);
  }
  else
  {
    size_t numAccepted = newPrior_.at("s").size();
    
    // to draw a random value in the sorted vector of pilot pilotPosteriors
    std::uniform_int_distribution<size_t> distVals(0, numAccepted - 1);
    size_t pos = distVals(gen);
    
    if(pos == 0) // if first pos in the vector there's no left neighbor
    {
      std::uniform_real_distribution<double> draw(newPrior_.at("s")[pos],
                                                  newPrior_.at("s")[pos + 1]);
      return draw(gen);
    }
    else if(pos == numAccepted - 1) // if last pos in the vector there's no right neighbor
    {
      std::uniform_real_distribution<double> draw(newPrior_.at("s")[pos],
                                                  newPrior_.at("s")[pos - 1]);
      return draw(gen);
    }
    else
    {
      std::uniform_real_distribution<double> draw(newPrior_.at("s")[pos - 1],
                                                  newPrior_.at("s")[pos + 1]);
      return draw(gen);
    }    
  }
}
