/*
 * Authors: Gustavo V. Barroso
 * Created: 01/11/2019
 * Last modified: 12/11/2020
 *
 */


#include <stdlib.h> 
#include <algorithm>

#include "HaploidGenome.h"
#include "Vcf.h"
#include "RecombinationMap.h"


// called once for mothers, another for fathers and another for children
std::vector<std::pair<HaploidGenome, HaploidGenome>> Vcf::parseSeqs(
const RecombinationMap& femMap, const RecombinationMap& maleMap) 
{
  checkBasicInfo_();
  
  ///////////////////////////////////////////////////////////
  boost::iostreams::filtering_istream boostStream; 
  
  if(compressionType_ == "gzip")
    boostStream.push(boost::iostreams::gzip_decompressor());
  else if(compressionType_ == "zip") 
    boostStream.push(boost::iostreams::zlib_decompressor());
  else if(compressionType_ == "bzip2")
    boostStream.push(boost::iostreams::bzip2_decompressor());
  else if(compressionType_ != "none")
    throw bpp::Exception("TIDES::Vcf::mis-specified sequence compression type!");

  std::ifstream inFile(fileName_, std::ios_base::in | std::ios_base::binary);  

  if(!inFile.is_open())
    throw bpp::Exception("TIDES::Vcf::could not open seq. file: " + fileName_); 
  else
    boostStream.push(inFile);
  ///////////////////////////////////////////////////////////
  
  // Haplotype type is defined in HaploidGenome.h
  std::vector<Haplotype> leftHaplotypes(numDiploids_);
  std::vector<Haplotype> rightHaplotypes(numDiploids_);

  // for a given chr
  std::vector<std::vector<size_t>> leftChr(numDiploids_,
std::vector<size_t>(0));
  std::vector<std::vector<size_t>> rightChr(numDiploids_,
std::vector<size_t>(0));

  std::string line, chr, prevChr = "";
  std::vector<std::string> splitLine(0);
  std::vector<std::string> prevLine(0);

  while(std::getline(boostStream, line)) 
  {  
    if(line.at(0) != '#') // gets lines with data  
    {    
      boost::split(splitLine, line, [](char c) { return c == '\t'; });
      chr = splitLine[0];      
      
      if(prevLine.size() == 0) // if first pos. in VCF 
        prevChr = chr;
      else
        prevChr = prevLine[0];
            
      if(chr == prevChr)
        addSite_(splitLine, leftChr, rightChr);
      else
      {
        for(size_t i = 0; i < numDiploids_; ++i)
        {
          // stores info about prev chr (std::move 'cleans' vectors so addSite_ works)
          leftHaplotypes[i].emplace(std::make_pair(prevChr, std::move(leftChr[i])));
          rightHaplotypes[i].emplace(std::make_pair(prevChr, std::move(rightChr[i])));
        }
        // starts seq for new chr
        addSite_(splitLine, leftChr, rightChr);
      }
      // moves along
      prevLine = splitLine;
    }
  }
  // the last chr
  for(size_t i = 0; i < numDiploids_; ++i)
  {
    leftHaplotypes[i].emplace(std::make_pair(chr, std::move(leftChr[i])));
    rightHaplotypes[i].emplace(std::make_pair(chr, std::move(rightChr[i])));
  }

  // safety first
  checkChrlabels_(leftHaplotypes.front(), femMap.getRanges(), maleMap.getRanges());

  std::vector<std::pair<HaploidGenome, HaploidGenome>> hapList(0);
  hapList.reserve(numDiploids_);
  for(size_t i = 0; i < numDiploids_; ++i)
    hapList.emplace_back(std::make_pair(HaploidGenome(std::move(leftHaplotypes[i])),
                                        HaploidGenome(std::move(rightHaplotypes[i]))));

  return hapList;
}

void Vcf::addSite_(const std::vector<std::string>& splitLine,
                   std::vector<std::vector<size_t>>& leftChr,
                   std::vector<std::vector<size_t>>& rightChr)
{
  size_t siteCoord = static_cast<size_t>(std::stoi(splitLine[1]));
  
  for(size_t i = 0; i < numDiploids_; ++i)
  { 
    if(splitLine[i + numInfoCols_][0] == '1')
      leftChr[i].push_back(siteCoord);
 
    if(splitLine[i + numInfoCols_][2] == '1')
      rightChr[i].push_back(siteCoord);
  }
}

void Vcf::checkChrlabels_(const Haplotype& chromosomes, const Ranges& frmap, const Ranges& mrmap)
{
  auto itVcf = std::begin(chromosomes);
  auto itFrmap = std::begin(frmap);
  auto itMrmap = std::begin(mrmap);

  std::string labVcf, labFrmap, labMrmap = "";

  for(; itVcf != std::end(chromosomes) && itFrmap != std::end(frmap) && itMrmap
!= std::end(mrmap); ++itVcf, ++itFrmap, ++itMrmap)
  {
    labVcf = itVcf->first;
    labFrmap = itFrmap->first;
    labMrmap = itMrmap->first;

    if((labFrmap != labMrmap) || (labFrmap != labVcf))
    {
      std::string msg = "VCF: " + labVcf + " Fem. map: " + labFrmap + "; Male map: " + labMrmap;
      throw bpp::Exception("TIDES::Chromosome labels from VCF and recombination maps don't match: " + msg);
    }
  }
}

void Vcf::checkBasicInfo_() 
{  
  ///////////////////////////////////////////////////////////
  boost::iostreams::filtering_istream boostStream; 
  
  if(compressionType_ == "gzip")
    boostStream.push(boost::iostreams::gzip_decompressor());
  else if(compressionType_ == "zip") 
    boostStream.push(boost::iostreams::zlib_decompressor());
  else if(compressionType_ == "bzip2")
    boostStream.push(boost::iostreams::bzip2_decompressor());
  else if(compressionType_ != "none")
    throw bpp::Exception("TIDES::Vcf::mis-specified sequence compression type!");

  std::ifstream inFile(fileName_, std::ios_base::in | std::ios_base::binary);  

  if(!inFile.is_open())
    throw bpp::Exception("TIDES::Vcf::could not open seq. file: " + fileName_); 
  else
    boostStream.push(inFile);
  ///////////////////////////////////////////////////////////
  
  std::string line = "";  
  std::vector<std::string> splitLine(0);

  while(std::getline(boostStream, line)) 
  {
    if(line[0] == '#' && line[1] != '#') //header  
    { 
      boost::split(splitLine, line, [](char c){ return c == '\t'; }); 
      
      //checks VCF format to update number of columns before genotypes      
      if(std::any_of(std::begin(splitLine), std::end(splitLine),
[&](const std::string& s) { return s == "FORMAT"; }))
        numInfoCols_ = 9; 
      
      //assumes there are no cols right of genotypes:
      numDiploids_ = splitLine.size() - numInfoCols_; 
      
      // stores individual ids
      auto startIds = std::begin(splitLine) + static_cast<std::vector< 
std::string>::difference_type >(numInfoCols_);

      auto endIds = std::end(splitLine);
      
      names_.reserve(numDiploids_);
      names_.insert(std::end(names_), startIds, endIds);
    }
    else if(line[0] != '#') // first line with genotypes
    {
      boost::split(splitLine, line, [](char c){ return c == '\t'; });
      
      for(size_t i = 0; i < numDiploids_; ++i)
        if(splitLine[i + numInfoCols_][1] == '/')
          throw bpp::Exception("TIDES::Character '/' in genotypes indicates that VCF is not phased!");
      //
      break;
    }
  }
}
