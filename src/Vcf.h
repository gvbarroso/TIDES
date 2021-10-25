/*
 * Authors: Gustavo V. Barroso
 * Created: 01/11/2019
 * Last modified: 25/05/2020
 *
 */


#ifndef _VCF_
#define _VCF_

#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <memory>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include "BppTextTools.h"
#include "BppExceptions.h"
#include "HaploidGenome.h"
#include "RecombinationMap.h"


class Vcf {
private:
  size_t numDiploids_;
  size_t numInfoCols_;
  
  std::string fileName_;
  std::string compressionType_;
  
  std::vector<std::string> names_; // stores indv ids from VCF header
  
public:
  Vcf(const std::string& name = "NA", const std::string& type = "NA"):
  numDiploids_(0),
  numInfoCols_(8),
  fileName_(name),
  compressionType_(type),
  names_(0)
  { }

public:
  size_t getNumberOfInfoColumns() 
  {
    return numInfoCols_;
  }
  
  size_t getNumberOfDiploids() 
  {
    return numDiploids_;
  }
  
  const std::string& getFileName() const
  {
    return fileName_;
  }
  
  void setFileName(const std::string& name)
  {
    fileName_ = name;
  }
  
  const std::string& getCompressionType() const 
  {
    return compressionType_;
  }
  
  void setCompressionType(const std::string& type)
  {
    compressionType_ = type;
  }
  
  const std::vector<std::string>& getIndvNames() const
  {
    return names_;
  }
  
  std::vector< std::pair<HaploidGenome, HaploidGenome>> parseSeqs(const RecombinationMap& femMap,
                                                                  const RecombinationMap& maleMap);

private:
  void addSite_(const std::vector<std::string>& splitLine, 
                std::vector<std::vector<size_t>>& leftSnps,
                std::vector<std::vector<size_t>>& rightSnps);

  void checkChrlabels_(const Haplotype& chromosomes, const Ranges& frmap, const Ranges& mrmap);
  
  void checkBasicInfo_();

 };

#endif
