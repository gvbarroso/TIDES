/*
 * Authors: Gustavo V. Barroso 
 * Created: 13/07/2020
 * Last modified: 18/11/2020
 *
 */


#ifndef _UTILS_
#define _UTILS_

#include <utility>
#include <chrono>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "BppExceptions.h"
#include "BppTextTools.h"
#include "Child.h"
#include "GenomicRanges.h"
#include "HaploidGenome.h"
#include "HashTable.h"
#include "Tools.h"

#define leftHap (*itIndv)->getHaplotypePair().first.getGenome()
#define rightHap (*itIndv)->getHaplotypePair().second.getGenome()

using namespace boost::accumulators;

// chr -> HashTable with SNP freqs.
typedef std::map<std::string, HashTable> SnpInfo;

class Utils {
  
public:
  // median absolute deviation
  static std::pair<double, double> mad(const std::vector<double>& vals) 
  {
    std::vector<double> cpy(0); // copy because of sorting
    std::copy(std::begin(vals), std::end(vals), std::back_inserter(cpy));
    
    size_t n = cpy.size();
    std::sort(std::begin(cpy), std::end(cpy));
    
    double median = -1.;
    if(n % 2 == 0)
      median = (cpy[n / 2] + cpy[(n / 2) - 1]) / 2.;
    else
      median = cpy[n / 2];
    
    std::vector<double> absDev(0);
    absDev.reserve(n);
    for(size_t i = 0; i < n; ++i)
      absDev.emplace_back(std::abs(cpy[i] - median));
 
    std::sort(std::begin(absDev), std::end(absDev));
  
    double mad = -1.;
    if(n % 2 == 0)
      mad = (absDev[n / 2] + absDev[(n / 2) - 1]) / 2.;
    else
      mad = absDev[n / 2]; 
    
    if(mad <= 0.)
      throw bpp::Exception("Utils::non-positive median absolute deviation!");
    
    return std::make_pair(median, mad);
  }
  
  static double pearson(const std::vector<double>& x, const std::vector<double>& y)
  {  
    accumulator_set<double, stats<tag::mean, tag::variance>> accX;
    std::for_each(std::begin(x), std::end(x), boost::bind<void>(boost::ref(accX), _1));
  
    accumulator_set<double, stats<tag::mean, tag::variance>> accY;
    std::for_each(std::begin(y), std::end(y), boost::bind<void>(boost::ref(accY), _1));
  
    double x_bar = mean(accX);
    double y_bar = mean(accY);
    
    double cov = 0.;
    for(size_t i = 0; i < x.size(); ++i)
      cov += (x[i] - x_bar) * (y[i] - y_bar);
    
    cov /= x.size();
    
    double rs = cov / (std::sqrt(variance(accX)) * std::sqrt(variance(accY)));
    return rs;
  }
  
  static double covar(const std::vector<double>& x, const std::vector<double>& y)
  {  
    accumulator_set<double, stats<tag::mean>> accX;
    std::for_each(std::begin(x), std::end(x), boost::bind<void>(boost::ref(accX), _1));
  
    accumulator_set<double, stats<tag::mean>> accY;
    std::for_each(std::begin(y), std::end(y), boost::bind<void>(boost::ref(accY), _1));
  
    double x_bar = mean(accX);
    double y_bar = mean(accY);
    
    double cov = 0.;
    for(size_t i = 0; i < x.size(); ++i)
      cov += (x[i] - x_bar) * (y[i] - y_bar);
    
    cov /= x.size();
    
    return cov;
  }
  
  static double watterson_theta(const std::vector<size_t>& sfs)
  {
    size_t num = 0;
    double denom = 0.;
    for(size_t i = 1; i < sfs.size(); ++i) // sfs[1] = snps with abs freq = 1
    {
      num += sfs[i];
      denom += 1. / i;
    }
    
    return static_cast<double>(num) / denom;
  }
  
  static double tajima_pi(const std::vector<size_t>& sfs)
  {
    size_t n = sfs.size();
    size_t denom = n * (n - 1);
    size_t num = 0;
    
    for(size_t i = 1; i < n; ++i) // sfs[1] = snps with abs freq = 1
      num += i * (n - i) * sfs[i];
    
    return (2. * static_cast<double>(num)) / static_cast<double>(denom);
  }
  
  static double tajima_norm(const std::vector<size_t>& sfs)
  {
    double n = static_cast<double>(sfs.size() - 1);
    
    double S = 0.;
    double a1 = 0.;
    double a2 = 0.;
    for(size_t i = 1; i < n - 1; ++i)
    {
      S += sfs[i];
      a1 += 1. / i;
      a2 += 1. / (i * i);
    }
    S += sfs[n - 1];
    
    double b1 = (n + 1.) / (3. * (n - 1.));
    double b2 = (2. * (n * n + n + 3.)) / (9. * n * (n - 1.));
    
    double c1 = b1 - (1. / a1);
    double c2 = b2 - (n + 2.) / (a1 * n) + a2 / (a1 * a1);
    
    double e1 = c1 / a1;
    double e2 = c2 / (a1 * a1 + a2);
    
    double var_d = e1 * S + e2 * S * (S - 1.);
    
    return std::sqrt(var_d);
  }
  
  static GenomicRanges merge_bedgraphs(const GenomicRanges& b1, const GenomicRanges& b2)
  {
    GenomicRanges ret;
    Ranges rangs;
    Values vals;
    
    // concatenates chr lists from both bedgraphs
    std::vector<std::string> chrList = b1.getChrList();
    chrList.insert(std::end(chrList), std::begin(b2.getChrList()),
std::end(b2.getChrList()));
    std::sort(std::begin(chrList), std::end(chrList));
    chrList.erase(std::unique(std::begin(chrList), std::end(chrList)),
std::end(chrList));
    
    ret.setChrList(chrList);
    
    for(auto itC = std::begin(chrList); itC != std::end(chrList); ++itC)
    { 
      // test existence of chr in both bedgraphs
      bool matchB1 = std::any_of(std::begin(b1.getChrList()),
std::end(b1.getChrList()), [&](const std::string& s) { return s == *itC; });
      
      bool matchB2 = std::any_of(std::begin(b2.getChrList()),
std::end(b2.getChrList()), [&](const std::string& s) { return s == *itC; });
      
      if(matchB1 && matchB2) // if chr is in both bedraphs 
      {
        Windows w(0); // for rangs
        std::vector<double> v(0); // for vals
        
        auto r1 = std::begin(b1.getRanges().at(*itC));
        auto v1 = std::begin(b1.getValues().at(*itC));
        auto r2 = std::begin(b2.getRanges().at(*itC));
        auto v2 = std::begin(b2.getValues().at(*itC));
        
        // (modifiable) copies of windows
        auto p1 = *r1;
        auto p2 = *r2; 
        
        // keeps track of advancing the iterators
        bool moved1 = true;
        bool moved2 = true;
        
        while((r1 != std::end(b1.getRanges().at(*itC))) || 
              (r2 != std::end(b2.getRanges().at(*itC))))
        {
          // checks whether windows have NOT been modified before copying
          if(moved1)
            p1 = *r1;
          if(moved2)
            p2 = *r2; 
        
          double avgVal = ((*v1) + (*v2)) / 2.; // avg between two bedgraphs
          
          // if bedgraph1 is exausted
          if(r1 == std::end(b1.getRanges().at(*itC)))
          {
            v.push_back(*v2);  
            w.push_back(std::make_pair(p2.first, p2.second));
            
            ++r2; ++v2;
            moved2 = true;
          }
          // if bedgraph2 is exausted
          else if(r2 == std::end(b2.getRanges().at(*itC)))
          {
            v.push_back(*v1);  
            w.push_back(std::make_pair(p1.first, p1.second));
            
            ++r1; ++v1;
            moved1 = true;
          }
          
          // otherwise there are windows to be merged
          
          else if(p1.first == p2.first) 
          {
            if(p1.second == p2.second) 
            {
              v.push_back(avgVal);  
              w.push_back(std::make_pair(p1.first, p1.second));
              
              ++r1; ++v1;
              ++r2; ++v2;
              moved1 = true; moved2 = true;
            }
            else if(p1.second < p2.second)
            {
              v.push_back(avgVal);
              w.push_back(std::make_pair(p1.first, p1.second));
              
              ++r1; ++v1;
              moved1 = true;
              if((*r1).first >= p2.second) 
              {
                v.push_back(*v2);  
                w.push_back(std::make_pair(p1.second, p2.second));
                
                ++r2; ++v2;
                moved2 = true;
              }
              else 
              {
                p2.first = p1.second;
                moved2 = false;
              }
            }
            else if(p1.second > p2.second)
            {
              v.push_back(avgVal); 
              w.push_back(std::make_pair(p1.first, p2.second));
              
              ++r2; ++v2;
              moved2 = true;
              
              if((*r2).first >= p1.second)
              {    
                v.push_back(*v1);
                w.push_back(std::make_pair(p2.second, p1.second));
                
                ++r1; ++v1; 
                moved1 = true;
              }
              else
              {
                p1.first = p2.second;
                moved1 = false;
              }
            }
          }
          else if(p1.first < p2.first)
          {
            if(p1.second == p2.second)
            {
              v.push_back(*v1); 
              v.push_back(avgVal);
              
              w.push_back(std::make_pair(p1.first, p2.first));
              w.push_back(std::make_pair(p2.first, p2.second));
              
              ++r1; ++v1;
              ++r2; ++v2;
              moved1 = true; moved2 = true;
            }
            else if(p1.second < p2.second)
            {
              if(p1.second < p2.first) // if windows don't overlap
              {
                v.push_back(*v1); 
                w.push_back(std::make_pair(p1.first, p1.second));
                
                ++r1; ++v1;
                moved1 = true;
                moved2 = false;
              }
              else
              {
                v.push_back(*v1); 
                v.push_back(avgVal);    
             
                w.push_back(std::make_pair(p1.first, p2.first));
                w.push_back(std::make_pair(p2.first, p1.second));
              
                ++r1; ++v1;
                moved1 = true;
              
                if((*r1).first >= p2.second)
                {
                  v.push_back(*v2);
                  w.push_back(std::make_pair(p1.second, p2.second));
                
                  ++r2; ++v2;
                  moved2 = true;
                }
                else
                {
                  p2.first = p1.second;
                  moved2 = false;
                }
              }
            }
            else if(p1.second > p2.second)
            {
              v.push_back(*v1); 
              v.push_back(avgVal);    
             
              w.push_back(std::make_pair(p1.first, p2.first));
              w.push_back(std::make_pair(p2.first, p2.second));
              
              ++r2; ++v2;
              moved2 = true;
              
              if((*r2).first >= p1.second)
              {
                v.push_back(*v1);
                w.push_back(std::make_pair(p2.second, p1.second));
                
                ++r1; ++v1;
                moved1 = true;
              }
              else
              {
                p1.first = p2.second;
                moved1 = false;
              }
            }
          }
          else if(p1.first > p2.first)
          { 
            if(p1.second == p2.second)
            {
              v.push_back(*v2);
              v.push_back(avgVal);
              
              w.push_back(std::make_pair(p2.first, p1.first));
              w.push_back(std::make_pair(p1.first, p1.second));
              
              ++r1; ++v1;
              ++r2; ++v2;
              moved1 = true; moved2 = true;
            }           
            else if(p1.second < p2.second)
            {
              v.push_back(*v2);    
              v.push_back(avgVal);
              
              w.push_back(std::make_pair(p2.first, p1.first));
              w.push_back(std::make_pair(p1.first, p1.second));
              
              ++r1; ++v1;
              moved1 = true;
              if((*r1).first >= p2.second)
              {
                v.push_back(*v1);  
                w.push_back(std::make_pair(p1.second, p2.second));
                
                ++r2; ++v2;
                moved2 = true;
              }
              else
              {
                p2.first = p1.second;
                moved2 = false;
              }
            }
            else if(p1.second > p2.second)
            {
              if(p1.first > p2.second) // if windows don't overlap
              {
                v.push_back(*v2); 
                w.push_back(std::make_pair(p2.first, p2.second));
                
                ++r2; ++v2;
                moved2 = true;
                moved1 = false;
              }
              else
              {
                v.push_back(*v2);    
                v.push_back(avgVal);
              
                w.push_back(std::make_pair(p2.first, p1.first));
                w.push_back(std::make_pair(p1.first, p2.second));
              
                ++r2; ++v2;
                moved2 = true;
                if((*r2).first >= p1.second)
                {
                  v.push_back(*v1); 
                  w.push_back(std::make_pair(p2.second, p1.second));
                
                  ++r1; ++v1;
                  moved1 = true;
                }
                else
                {
                  p1.first = p2.second;
                  moved1 = false;
                }
              }
            }
          }
        }
        rangs.emplace(std::make_pair(*itC, std::move(w)));
        vals.emplace(std::make_pair(*itC, std::move(v)));
      }
      else if(matchB1 && !matchB2)
      {
        // copy all info
        rangs.emplace(std::make_pair(*itC, b1.getRanges().at(*itC)));
        vals.emplace(std::make_pair(*itC, b1.getValues().at(*itC)));
      }
      else if(matchB2 && !matchB1)
      {
        // copy all info
        rangs.emplace(std::make_pair(*itC, b2.getRanges().at(*itC)));
        vals.emplace(std::make_pair(*itC, b2.getValues().at(*itC)));
      }
    }
    
    ret.setRanges(std::move(rangs));
    ret.setValues(std::move(vals));
    
    return ret;
  }
  
  static void write_sfs(const std::vector<size_t>& sfs, const std::string& name)
  {
    std::ofstream sfsFile;
    sfsFile.open(name);
  
    if(sfsFile.is_open())
    { 
      sfsFile << "Freq.\tCount" << std::endl; // header
    
      for(size_t i = 0; i < sfs.size(); ++i)
        sfsFile << i << "\t" << sfs[i] << std::endl;
    }  
    else
      throw bpp::Exception("Could not write SFS to file: " + name);
  
    sfsFile.close();
  }
  
  // sampleSize = number of individuals (diploids)
  static std::vector<size_t> fetch_sfs(size_t sampleSize, SnpInfo& freqs)
  {
    // SNPs with 0 freq. due to 1) loss from parents or 2) [removed] de novo mutations;
    std::vector<size_t> sfs(2 * sampleSize, 0); // there are 2N copies of each allele
  
    for(auto itChr = std::begin(freqs); itChr != std::end(freqs); ++itChr)
    {
      for(size_t i = 0; i < itChr->second.getTableSize(); ++i)
      {
        HashNode* entry = itChr->second.getTable()[i];
      
        while(entry != NULL)
        {  
          size_t f = entry->getValue(); 
          ++sfs[f];
        
          entry = entry->getNext(); // moves to next node in same bucket
        }
      }
    }
    return sfs;
  }
  
  static SnpInfo fetch_snp_freqs(const Haplotype& refSnps,
const std::vector<std::shared_ptr<Diploid>>& indvs)
  { 
    // build from refSnps, hence contains SNPs with 0 frequency in kids || moms || dads
    SnpInfo freqs;
  
    std::string chr = "";
    for(auto itChr = std::begin(refSnps); itChr != std::end(refSnps); ++itChr)
    {
      chr = itChr->first;
      size_t tableSize = itChr->second.size() + itChr->second.size() / 2; // 1.5X OG size
  
      HashTable table(tableSize, itChr->second.back());
      table.init(itChr->second); // built from all refSnps_ here
  
      for(auto itIndv = std::begin(indvs); itIndv != std::end(indvs); ++itIndv)
      { 
        for(auto itL = std::begin(leftHap.at(chr)); itL != std::end(leftHap.at(chr)); ++itL)
          table.increment(*itL, 1.); 
    
        for(auto itR = std::begin(rightHap.at(chr)); itR != std::end(rightHap.at(chr)); ++itR)
          table.increment(*itR, 1.);
      }
      // 
      freqs.emplace(std::make_pair(std::move(chr), std::move(table))); 
    } 
  
    return freqs;
  }

  static void compute_snp_freqs(SnpInfo& f, const std::vector<std::shared_ptr<Diploid>>& d)
  { 
    std::string chr = "";  
    for(auto itChr = std::begin(f); itChr != std::end(f); ++itChr)
    {
      chr = itChr->first;
      itChr->second.reset(); // sets all SNP frequencies to 0
  
      for(auto itIndv = std::begin(d); itIndv != std::end(d); ++itIndv)
      { 
        for(auto itL = std::begin(leftHap.at(chr)); itL != std::end(leftHap.at(chr)); ++itL)
        itChr->second.increment(*itL, 1.); 
    
        for(auto itR = std::begin(rightHap.at(chr)); itR != std::end(rightHap.at(chr)); ++itR)
          itChr->second.increment(*itR, 1.);
      }
    } 
  }
  
  // NOTE freqParents should be passed as const here
  static void remove_denovo_snps(std::vector<std::shared_ptr<Child>>& kids,
SnpInfo& freqParents)
  {
    Tools pb; // this function can be parallelized
  
    size_t cc = 0; // pb helper  
    size_t x = 0; // de novo mutations counter
  
    std::string chr = "";
    for(auto itChr = std::begin(freqParents); itChr != std::end(freqParents); ++itChr)
    {
      chr = itChr->first;
    
      for(size_t i = 0; i < itChr->second.getTableSize(); ++i) 
      {
        HashNode* entry = itChr->second.getTable()[i];
      
        while(entry != NULL)
        {  
          if(entry->getValue() == 0.) // if SNP has freq. 0 in parents
          {
            size_t coord = entry->getKey(); // accesses the SNP coordinate
          
            for(auto itIndv = std::begin(kids); itIndv != std::end(kids); ++itIndv)
            {
              auto itSnp = std::lower_bound(std::begin(leftHap.at(chr)),
                                            std::end(leftHap.at(chr)), coord);
            
              if(itSnp != std::end(leftHap.at(chr)) && *itSnp == coord)
              { 
                ++x;
                leftHap.at(chr).erase(itSnp); // removes SNP if present in focal child
                (*itIndv)->countSites(); // updates counts of sites
                break; // breaks since de novo mutations have abs. freq. = 1
              }
          
              // if we didn't find the SNP in 1st Haplotype, search 2nd
              itSnp = std::lower_bound(std::begin(rightHap.at(chr)),
                                       std::end(rightHap.at(chr)), coord);
          
              if(itSnp != std::end(rightHap.at(chr)) && *itSnp == coord)
              {   
                ++x;
                rightHap.at(chr).erase(itSnp); // removes SNP if present in focal child
                (*itIndv)->countSites(); // updates counts of sites
                break; // breaks since de novo mutations have abs. freq. = 1
              }
            }
          }
          // moves to next node in same bucket
          entry = entry->getNext(); 
        }
      }
    
      double prog = static_cast<double>(cc + 1) / static_cast<double>(freqParents.size());
      pb.print_pb(prog);
    
      ++cc;
    }
    //
    std::cout << "\n[removed " << x << " de novo SNPs.]\n";
  }
  
    // NOTE info should be passed as const here
  static void write_snp_info(const Haplotype& refSnps, SnpInfo& info,
const std::string& name)
  {
    std::ofstream freqsFile;
    freqsFile.open(name, std::ios_base::out | std::ios_base::binary);  
  
    boost::iostreams::filtering_ostream freqStream; 
    freqStream.push(boost::iostreams::gzip_compressor());
    freqStream.push(freqsFile);  
  
    freqStream << "chr\tpos\tfreq" << std::endl;

    for(auto itChr = std::begin(refSnps); itChr != std::end(refSnps); ++itChr)
    {
      for(size_t i = 0; i < itChr->second.size(); ++i)
        freqStream << itChr->first << "\t" << itChr->second[i] << "\t" << info.at(itChr->first).getValue(itChr->second[i]) << std::endl;
    }
    //
    boost::iostreams::close(freqStream);
  }

};

#endif
