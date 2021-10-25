/*
 * Authors: Gustavo V. Barroso 
 * Created: 23/10/2019
 * Last modified: 21/07/2020
 *
 */

#include <iostream>
#include <fstream>
#include <numeric>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "BppApplication.h"
#include "BppTextTools.h"
#include "GenomicRanges.h"


double GenomicRanges::fetchValue(const std::string& chr, unsigned int pos) 
{
  // naive search, not optimised
  auto itW = std::begin(ranges_.at(chr));
  auto itV = std::begin(values_.at(chr));

  for(; itW != std::end(ranges_.at(chr)) && itV != std::end(values_.at(chr)); ++itW, ++itV) 
  {
    if(pos >= itW -> first && pos <= itW -> second)
      return *itV;
  }

  throw bpp::Exception("TIDES::GenomicRanges::could not find value for pos: " + bpp::TextTools::toString(pos));
}

const double GenomicRanges::fetchValue(const std::string& chr, unsigned int pos) const 
{
  // naive search, not optimised
  auto itR = std::begin(ranges_.at(chr));
  auto itV = std::begin(values_.at(chr));

  for(; itR != std::end(ranges_.at(chr)) && itV != std::end(values_.at(chr)); ++itR, ++itV)
  {
    if(pos >= itR->first && pos <= itR->second)
      return *itV;
  }

  throw bpp::Exception("TIDES::GenomicRanges::could not find value for pos: " +
bpp::TextTools::toString(pos));
}

void GenomicRanges::writeToFile(const std::string& name)
{
  boost::iostreams::filtering_ostream stream; 
  
  std::ofstream file;
  file.open(name, std::ios_base::out | std::ios_base::binary);
  
  stream.push(boost::iostreams::gzip_compressor());
  stream.push(file);  
  
  auto itR = std::begin(ranges_);
  auto itV = std::begin(values_);
  
  for(; itR != std::end(ranges_); ++itR, ++itV)
  {
    auto itW = std::begin(itR->second);
    auto itD = std::begin(itV->second);
    
    for(; itW != std::end(itR->second); ++itW, ++itD)
    {
      stream << itR->first << "\t" << (*itW).first << "\t";
      stream << (*itW).second << "\t" << *itD << std::endl; 
    }
  }

  boost::iostreams::close(stream);
}

void GenomicRanges::parseBedgraph_(const std::string& file) 
{       
  std::ifstream bedgraphFile;
  bedgraphFile.open(file, std::ios::in);
  
  if(bedgraphFile.is_open()) 
  {  
    bool firstLine = true;
    unsigned int startCoord, endCoord = 0;
    double value = 0.;
    std::string focalChr, prevChr, line = "";
    std::vector<std::string> splitLine(0);
    std::vector<std::pair<unsigned int, unsigned int>> chrRanges(0);
    std::vector<double> chrValues(0);

    while(std::getline(bedgraphFile, line)) 
    {             
      boost::split(splitLine, line, [](char c) { return c == '\t'; });
      
      focalChr = splitLine[0];
      startCoord = std::stoul(splitLine[1]);
      endCoord = std::stoul(splitLine[2]);
      value = std::stod(splitLine[3]); 
      
      if(firstLine) 
      {
        prevChr = focalChr;
        firstLine = false;
      }
       
      if(focalChr == prevChr) // keep pushing
      { 
        chrRanges.push_back(std::make_pair(startCoord, endCoord));
        chrValues.push_back(value);
      }

      else // store chr info and reset
      { 
        chrList_.push_back(prevChr);
            
        ranges_.emplace(std::make_pair(prevChr, std::move(chrRanges)));
        values_.emplace(std::make_pair(prevChr, std::move(chrValues)));

        chrRanges.push_back(std::make_pair(startCoord, endCoord));
        chrValues.push_back(value);

        prevChr = focalChr;
      }
    }
    // last chr misses the comparisons above
    chrList_.push_back(prevChr);
    ranges_.emplace(std::make_pair(prevChr, std::move(chrRanges)));
    values_.emplace(std::make_pair(prevChr, std::move(chrValues)));

    bedgraphFile.close();
  }
    
  else
    throw bpp::Exception("TIDES::GenomicRanges::could not open bedgraph:" + file);
}

void GenomicRanges::genCdf_()
{
  for(auto itV = std::begin(values_); itV != std::end(values_); ++itV) // for all chr
  {
    std::vector<double> focalCdf(itV->second);
    
    for(size_t i = 1; i < focalCdf.size(); ++i) // for al windows in chr
      focalCdf[i] += focalCdf[i - 1];  // accumulates
    
    // stores
    cdf_.emplace(std::make_pair(itV->first, std::move(focalCdf)));
  } 
}

