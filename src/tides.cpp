/*
 * Authors: Gustavo V. Barroso 
 * Created: 22/10/2019
 * Last modified: 10/11/2020
 *
 */


#include <string>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <utility>
#include <numeric>
#include <algorithm>
#include <memory>
#include <map>
#include <cmath>
#include <thread>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include "BppApplication.h"
#include "BppApplicationTools.h"
#include "BppTextTools.h"
#include "GenomicRanges.h"
#include "RealDiploid.h"
#include "Child.h"
#include "Vcf.h"
#include "Simulator.h"
#include "Tools.h"
#include "Utils.h"
#include "Global.h"

size_t NUM_AVAIL_THREADS; // extern variable in Global.h

int main(int argc, char *argv[]) { 
    
  
  std::cout << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                       TIDES version 0.0.1                      *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*         Inferring Natural Selection from Family Trios          *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "* Authors: Gustavo V. Barroso         Created:       22/Oct/2019 *" << std::endl;
  std::cout << "*          Kirk E. Lohmueller         Last Modified: 21/Dec/2020 *" << std::endl;
  std::cout << "*                                                                *" << std::endl;
  std::cout << "******************************************************************" << std::endl;
  std::cout << std::endl;

  if(argc == 1) {
    std::cout << "To use TIDES in Terminal => tides params=params_file.bpp\n\n";
    
    std::cout << "**********************************************************\n\n";
    
    std::cout << "params_file.bpp is a text file with the following options:\n\n";
    
    std::cout << "# the type of file compression used in the VCFs [OPTIONS]\n";
    std::cout << "seq_compression = ['zip', 'bzip2', 'gzip' OR 'none']\n\n"; 
    std::cout << "# the (relative) paths to VCF files, e.g.:\n";
    std::cout << "mothers_file = mothers.vcf.gz\n";
    std::cout << "fathers_file = fathers.vcf.gz\n";
    std::cout << "children_file = children.vcf.gz\n\n";

    std::cout << "# the (relative) path to a tab-separated file indicating\n";
    std::cout << "# the trio 'relatedness' among indv. IDs (as found in the VCF header):\n";
    std::cout << "trio_ids_file = trio_ids.txt\n\n";
    std::cout << "# for example, trio_ids_file may begin with:\n";
    std::cout << "# 'child_id\tmother_id\tfather_id => [MANDATORY HEADER]\n";
    std::cout << "# i0\ti11632\ti17883\n";
    std::cout << "# i1\ti23651\ti1848\n";
    std::cout << "# i2\ti11211\ti12911\n";
    std::cout << "# i3\ti19837\ti14224\n";
    std::cout << "# i4\ti21042\ti12911'\n";
    std::cout << "# note that trios may be overlapping, e.g.\n";
    std::cout << "# child i2 and child i4 share the same father (half-sibs).\n\n";
    
    std::cout << "# the (relative) path to (sex-specific) recombination maps, e.g.:\n";
    std::cout << "female_rmap_path = ../frmap.bedgraph\n";
    std::cout << "male_rmap_path = ../mrmap.bedgraph\n";
    std::cout << "# if only a sex-average map is available, provide same path twice\n\n";

    std::cout << "mu = # de novo mutation rate per nucleotide\n";
    std::cout << "seq_len = # total seq. length (number of sites for mu)\n\n";
    
    std::cout << "# range of uniform priors in the pilot run:\n";
    std::cout << "# (comma-separated values within parenthesis)\n";
    std::cout << "s_interval = # DEFAULT = (1.0, 6.0) in log10 space\n";
    std::cout << "h_interval = # DEFAULT = (0.0, 0.5)\n";
    std::cout << "# to fix a parameter, provide identical values\n\n";
    
    std::cout << "num_threads = # integer, DEFAULT = use all available cores\n";
    std::cout << "num_pilot = # integer, size of pilot run to improve prior dist.\n";
    std::cout << "num_sims = # integer, number of simulations\n";
    std::cout << "num_meiosis = # integer, DEFAULT = 150\n";
    std::cout << "num_accepted = # integer, number of accepted pilot simulations\n\n";
    
    std::cout << "**********************************************************\n\n\n";
    
    std::cout << "For more information, please email gvbarroso@gmail.com\n\n";
    return(0);
  }  
  
  /////////////////////////////////////////////////////////////////////////////////////
  bpp::BppApplication tides(argc, argv, "TIDES");
  
  tides.startTimer();
  std::map<std::string, std::string> params = tides.getParams();

  std::string frmap = bpp::BppApplicationTools::getAFilePath(
"female_rmap_path", params, "none");

  std::string mrmap = bpp::BppApplicationTools::getAFilePath(
"male_rmap_path", params, "none");
  
  RecombinationMap femaleRecMap(frmap, "f");
  RecombinationMap maleRecMap(mrmap, "m");
  
  std::vector<double> selInter = bpp::BppApplicationTools::getVectorParameter<double>
("s_interval", params, ',', "(1.0,6.0)", "", true, 4);

  std::vector<double> domInter =  bpp::BppApplicationTools::getVectorParameter<double>
("h_interval", params, ',', "(0.0,0.6)", "", true, 4);
  ///////////////////////////////////////////////////////////////////////////////////////
  

  ///////////////////////////////////////////////////////////////////////////////////////
  size_t numThreads = bpp::BppApplicationTools::getParameter<size_t>("num_threads",
params, std::thread::hardware_concurrency(), "", true, 4);
  
  size_t numMeiosis = bpp::BppApplicationTools::getParameter<size_t>("num_meiosis",
params, 150, "", true, 4); 

  // updates global variable
  NUM_AVAIL_THREADS = numThreads;
  std::vector<std::shared_ptr<Child>> kids(0); // all kids
  Simulator abc(selInter, domInter);
  
  {
    /*
     * start of fake scope
     */

    ///////////////////////////////////////////////////////////////////////////////////
    std::string compression = bpp::BppApplicationTools::getStringParameter(
"seq_compression", params, "none"); 

    std::vector<std::string> seqFileNames(0);
    seqFileNames.push_back(bpp::BppApplicationTools::getStringParameter("mothers_file",
params, "")); 
    seqFileNames.push_back(bpp::BppApplicationTools::getStringParameter("fathers_file",
params, "")); 
    
    bool simKids = bpp::BppApplicationTools::getBooleanParameter("sim_kids", params,
false);
    
    if(!simKids)
      seqFileNames.push_back(bpp::BppApplicationTools::getStringParameter("children_file",
params, "")); 
    
    size_t numFiles = seqFileNames.size(); 
        
    // type [moms | dads | kids] -> individuals -> Haploid pairs
    std::vector<std::vector<std::pair<HaploidGenome, HaploidGenome>>> tmp(numFiles,
std::vector<std::pair<HaploidGenome, HaploidGenome>>(0));
    
    // the invdividual ids from the VCF header files, to assemble trios
    std::vector<std::vector<std::string>> indvIds(numFiles, std::vector<std::string>(0));
    
    std::cout << "Processing input files..."; std::cout.flush();
    
    auto parseFiles = [&] (size_t idx) 
    {
      Vcf vcfReader(seqFileNames[idx], compression);
      tmp[idx] = vcfReader.parseSeqs(femaleRecMap, maleRecMap); 
      indvIds[idx] = vcfReader.getIndvNames(); // from VCF header file
    };
  
    if(numThreads >= numFiles)
    {
      std::thread* threadVector = new std::thread[numFiles];
  
      for(size_t i = 0; i < numFiles; ++i)
        threadVector[i] = std::thread(parseFiles, i);
  
      for(size_t i = 0; i < numFiles; ++i)
        threadVector[i].join();
      
      delete [] threadVector;
    }
    else 
      for(size_t i = 0; i < numFiles; ++i)
        parseFiles(i);
    
    
    // Number of moms and dads may differ => 2 loops needed here:
    std::vector<std::shared_ptr<RealDiploid>> moms(0);
    moms.reserve(tmp[0].size());
    for(size_t i = 0; i < tmp[0].size(); ++i)
      moms.emplace_back(std::make_shared<RealDiploid>(RealDiploid(tmp[0][i], indvIds[0]
[i], "f")));  
    
    std::vector<std::shared_ptr<RealDiploid>> dads(0);
    dads.reserve(tmp[1].size());
    for(size_t i = 0; i < tmp[1].size(); ++i)
      dads.emplace_back(std::make_shared<RealDiploid>(RealDiploid(tmp[1][i], indvIds[1]
[i], "m")));
    
    std::cout << " done." << std::endl << std::endl;  
    ///////////////////////////////////////////////////////////////////////////////////
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////
    if(simKids) // if we are simulating the 'observed' data from fathers and mothers
    {
      double s = bpp::BppApplicationTools::getDoubleParameter("sim_s", params, 1.); 
      double h = bpp::BppApplicationTools::getDoubleParameter("sim_h", params, 0.5);
      size_t numKidsPerCouple = bpp::BppApplicationTools::getParameter<size_t>(
"kids_per_couple", params, 10); 
      size_t seed = bpp::BppApplicationTools::getParameter<size_t>("seed", params, 42);
      
      // set seed to get the same results when fitting multiple models to the same data
      // important for model comparison / selection!
      std::mt19937 gen;
      gen.seed(seed); 
      
      std::shuffle(std::begin(moms), std::end(moms), gen);
      std::shuffle(std::begin(dads), std::end(dads), gen);
      
      std::vector<double> fitMomi(moms.size());
      std::vector<double> fitDadi(dads.size());
      std::vector<double> fitKids(kids.size());
      
      std::uniform_real_distribution<double> runif(0., 1.); // for viability selection
      
      bool dfe = bpp::BppApplicationTools::getBooleanParameter("dfe", params, false, "",
true, 4); 
      
      // HashTables with selection coefficients (1 per chr)
      SnpInfo selCoeffs; // for dfe == true only, but declared here for conevience
      
      // if SNPs have different s values fitness computation is tricky
      if(dfe) 
      {
        Tools dfePb;
        std::cout << "Engineering kids with DFE from bedgraph files.\n" << std::endl;

        // bedgraphs with duplicate SNP coordinates (@ start, end) and s as 4th col.
        GenomicRanges sbMoms("dfe_mothers.bedgraph");
        GenomicRanges sbDads("dfe_fathers.bedgraph");
       
        GenomicRanges sbParents(std::move(Utils::merge_bedgraphs(sbMoms, sbDads)));
        
        std::cout << "Building reference SNP panel..." << std::endl;
        
        // builds reference SNPs for parental generation
        Haplotype refSnps = moms[0]->getHaplotypePair().first.getGenome(); // init
        for(size_t i = 0; i < moms.size(); ++i)
        {
          std::string chr = "";
        
          // chr iterators
          auto stop = std::end(moms[i]->getHaplotypePair().first.getGenome());

          auto m1 = std::begin(moms[i]->getHaplotypePair().first.getGenome());
          auto m2 = std::begin(moms[i]->getHaplotypePair().second.getGenome());
        
          auto f1 = std::begin(dads[i]->getHaplotypePair().first.getGenome());
          auto f2 = std::begin(dads[i]->getHaplotypePair().second.getGenome());
        
          for(; m1 != stop; ++m1, ++m2, ++f1, ++f2)
          {
            chr = m1->first; // shared by all iterators
          
            refSnps.at(chr).insert(std::end(refSnps.at(chr)),
                                   std::begin(m1->second),
                                   std::end(m1->second));
           
            refSnps.at(chr).insert(std::end(refSnps.at(chr)), 
                                   std::begin(m2->second),
                                   std::end(m2->second));
          
            refSnps.at(chr).insert(std::end(refSnps.at(chr)),
                                   std::begin(f1->second),
                                   std::end(f1->second));
          
            refSnps.at(chr).insert(std::end(refSnps.at(chr)), 
                                   std::begin(f2->second),
                                   std::end(f2->second));
          }
          
          double prog = static_cast<double>(i + 1) / static_cast<double>(moms.size());
          dfePb.print_pb(prog);
        }
        
        // deletes duplicates
        for(auto it = std::begin(refSnps); it != std::end(refSnps); ++it)
        {
          std::sort(std::begin(it->second), std::end(it->second));
          it->second.erase(std::unique(std::begin(it->second), std::end(it->second)),
                           std::end(it->second));
        }
        
        dfePb.reset_pb();
        std::cout << std::endl << std::endl;
        
        std::cout << "Assigning selection coefficients to SNPs...";
        std::cout.flush();
        
        // init s (to 0) for every SNP in dataset
        for(auto itChr = std::begin(refSnps); itChr != std::end(refSnps); ++itChr)
        {
          std::string chr = itChr->first;  
          
          size_t tableSize = itChr->second.size() + itChr->second.size() / 2;
          
          HashTable table(tableSize, itChr->second.back());
          table.init(itChr->second); 
          
          selCoeffs.emplace(std::make_pair(std::move(chr), std::move(table))); 
        }
        
        // updates s on hash tables
        for(auto itC = std::begin(refSnps); itC != std::end(refSnps); ++itC)
        {
          for(auto itSnp = std::begin(itC->second); itSnp != std::end(itC->second); ++itSnp)
          {
            double focalS = sbParents.fetchValue(itC->first, *itSnp);
            selCoeffs.at(itC->first).increment(*itSnp, focalS);
          }
        }
        std::cout << "done." << std::endl << std::endl;
        
        std::cout << "Computing fitness for females..." << std::endl;
        
        for(size_t i = 0; i < moms.size(); ++i)
        {
          double fit = 1.;
          // left haplotype
          for(auto itC = std::begin(moms[i]->getHaplotypePair().first.getGenome());
itC != std::end(moms[i]->getHaplotypePair().first.getGenome()); ++itC)
          {
            std::string chr = itC->first;
            for(auto itS = std::begin(moms[i]->getHaplotypePair().first.getGenome().at(chr));
itS != std::end(moms[i]->getHaplotypePair().first.getGenome().at(chr)); ++itS)
            {
              bool homo = std::binary_search(std::begin(
moms[i]->getHaplotypePair().second.getGenome().at(chr)), std::end(
moms[i]->getHaplotypePair().second.getGenome().at(chr)), *itS);
              // if homozygous
              if(homo)
                fit *= (1. + selCoeffs.at(chr).getValue(*itS));
              // if heterozygous
              else
                fit *= (1. + selCoeffs.at(chr).getValue(*itS) * h);
            }
          }
          // right haplotype
          for(auto itC = std::begin(moms[i]->getHaplotypePair().second.getGenome());
itC != std::end(moms[i]->getHaplotypePair().second.getGenome()); ++itC)
          {
            std::string chr = itC->first;
            for(auto itS = std::begin(moms[i]->getHaplotypePair().second.getGenome().at(chr));
itS != std::end(moms[i]->getHaplotypePair().second.getGenome().at(chr)); ++itS)
            {
              bool homo = std::binary_search(std::begin(
moms[i]->getHaplotypePair().first.getGenome().at(chr)), std::end(
moms[i]->getHaplotypePair().first.getGenome().at(chr)), *itS);
              // second set of chromosomes, we only look for heterozygous pos.
              if(!homo)
                fit *= (1. + selCoeffs.at(chr).getValue(*itS) * h);
            }
          }
          fitMomi[i] = fit;
          
          double prog = static_cast<double>(i + 1) / static_cast<double>(moms.size());
          dfePb.print_pb(prog);
        }
        dfePb.reset_pb();
        std::cout << std::endl << std::endl;
        
        std::cout << "Computing fitness for males..." << std::endl;
        
        for(size_t i = 0; i < dads.size(); ++i)
        {
          double fit = 1.;
          // left haplotype
          for(auto itC = std::begin(dads[i]->getHaplotypePair().first.getGenome());
itC != std::end(dads[i]->getHaplotypePair().first.getGenome()); ++itC)
          {
            std::string chr = itC->first;
            for(auto itS = std::begin(dads[i]->getHaplotypePair().first.getGenome().at(chr));
itS != std::end(dads[i]->getHaplotypePair().first.getGenome().at(chr)); ++itS)
            {
              bool homo = std::binary_search(std::begin(
dads[i]->getHaplotypePair().second.getGenome().at(chr)), std::end(
dads[i]->getHaplotypePair().second.getGenome().at(chr)), *itS);
              // if homozygous
              if(homo)
                fit *= (1. + selCoeffs.at(chr).getValue(*itS));
              // if heterozygous
              else
                fit *= (1. + selCoeffs.at(chr).getValue(*itS) * h);
            }
          }
          // right haplotype
          for(auto itC = std::begin(dads[i]->getHaplotypePair().second.getGenome());
itC != std::end(dads[i]->getHaplotypePair().second.getGenome()); ++itC)
          {
            std::string chr = itC->first;
            for(auto itS = std::begin(dads[i]->getHaplotypePair().second.getGenome().at(chr));
itS != std::end(dads[i]->getHaplotypePair().second.getGenome().at(chr)); ++itS)
            {
              bool homo = std::binary_search(std::begin(
dads[i]->getHaplotypePair().first.getGenome().at(chr)), std::end(
dads[i]->getHaplotypePair().first.getGenome().at(chr)), *itS);
              // second set of chromosomes, we only look for heterozygous pos.
              if(!homo)
                fit *= (1. + selCoeffs.at(chr).getValue(*itS) * h);
            }
          }
          fitDadi[i] = fit;
          
          double prog = static_cast<double>(i + 1) / static_cast<double>(dads.size());
          dfePb.print_pb(prog);
        }
      }
      else
      {
        std::cout << "Engineering kids with s = " << s << "; h = " << h << "..."; 
        std::cout.flush();
      
        // fitness moms
        for(size_t i = 0; i < moms.size(); ++i)
        { 
          double fitHomo = std::exp(-(moms[i]->getNumHomoSites() * s));
          double fitHet = std::exp(-(moms[i]->getNumHetSites() * h * s));
        
          fitMomi[i] = fitHomo * fitHet;
        }
      
        // fitness dads
        for(size_t i = 0; i < dads.size(); ++i)
        { 
          double fitHomo = std::exp(-(dads[i]->getNumHomoSites() * s));
          double fitHet = std::exp(-(dads[i]->getNumHetSites() * h * s));
        
          fitDadi[i] = fitHomo * fitHet;
        }
      }
      
      // viability selection in moms
      auto itM = std::begin(moms);
      auto itF = std::begin(fitMomi);
      while(itM != std::end(moms))
      {
        double r = runif(gen);
      
        if(*itF < r)
        {
          itM = moms.erase(itM);
          itF = fitMomi.erase(itF);
        }
        else
        {
          ++itM;
          ++itF;
        }
      }
    
      // viability selection in dads
      auto itD = std::begin(dads);
      auto itG = std::begin(fitDadi);
      while(itD != std::end(dads))
      {
        double r = runif(gen);
      
        if(*itG < r)
        {
          itD = dads.erase(itD);
          itG = fitDadi.erase(itG);
        }
        else
        {
          ++itD;
          ++itG;
        }
      }
      
      size_t numCouples = std::min(moms.size(), dads.size());
      if(numCouples < 1)
        throw bpp::Exception("TIDES::No adults survided viability selection!");
      
      kids.reserve(numCouples * numKidsPerCouple);
    
      // matches males and females at random and generates kids
      for(size_t i = 0; i < numCouples; ++i) // non-overlapping couples
      {
        for(size_t j = 0; j < numKidsPerCouple; ++j)
        {
          // for reproducibility (when simulating data and doing model selection)
          moms[i]->setSeed(seed + j); 
          dads[i]->setSeed(seed + j);
          
          auto z = std::make_pair(HaploidGenome(moms[i]->generateGamete(femaleRecMap)),
                                  HaploidGenome(dads[i]->generateGamete(maleRecMap)));
         
          std::string name = "k" + bpp::TextTools::toString(j + i * numKidsPerCouple);
          kids.emplace_back(std::make_shared<Child>(Child(z, moms[i], dads[i], name)));
        }
        moms[i]->resetSeed(); 
        dads[i]->resetSeed(); 
      }
    
      fitKids.resize(kids.size());
      if(dfe)
      {
        std::cout << std::endl << std::endl << "Computing fitness for kids..." << std::endl;
        Tools dfePb2;
        
        for(size_t i = 0; i < kids.size(); ++i)
        {
          double fit = 1.;
          for(auto itC = std::begin(kids[i]->getHaplotypePair().first.getGenome());
itC != std::end(kids[i]->getHaplotypePair().first.getGenome()); ++itC)
          {
            std::string chr = itC->first;
            for(auto itS = std::begin(kids[i]->getHaplotypePair().first.getGenome().at(chr));
itS != std::end(kids[i]->getHaplotypePair().first.getGenome().at(chr)); ++itS)
            {
              bool homo = std::binary_search(std::begin(
kids[i]->getHaplotypePair().second.getGenome().at(chr)), std::end(
kids[i]->getHaplotypePair().second.getGenome().at(chr)), *itS);
              
              // if homozygous
              if(homo)
                fit *= (1. + selCoeffs.at(chr).getValue(*itS));
              // if heterozygous
              else
                fit *= (1. + selCoeffs.at(chr).getValue(*itS) * h);
            }
          }
          for(auto itC = std::begin(kids[i]->getHaplotypePair().second.getGenome());
itC != std::end(kids[i]->getHaplotypePair().second.getGenome()); ++itC)
          {
            std::string chr = itC->first;
            for(auto itS = std::begin(kids[i]->getHaplotypePair().second.getGenome().at(chr));
itS != std::end(kids[i]->getHaplotypePair().second.getGenome().at(chr)); ++itS)
            {
              bool homo = std::binary_search(std::begin(
kids[i]->getHaplotypePair().first.getGenome().at(chr)), std::end(
kids[i]->getHaplotypePair().first.getGenome().at(chr)), *itS);
              // second set of chromosomes, we only look for heterozygous pos.
              if(!homo)
                fit *= (1. + selCoeffs.at(chr).getValue(*itS) * h);
            }
          }
          fitKids[i] = fit;
          
          double prog = static_cast<double>(i + 1) / static_cast<double>(kids.size());
          dfePb2.print_pb(prog);
        }
      }
      else
      {
        for(size_t i = 0; i < kids.size(); ++i)
        { 
          double fitHomo = std::exp(-(kids[i]->getNumHomoSites() * s));
          double fitHet = std::exp(-(kids[i]->getNumHetSites() * h * s));
        
          fitKids[i] = fitHomo * fitHet;
        }
      }
      
      // viability selection in kids
      auto itK = std::begin(kids);
      auto itJ = std::begin(fitKids);
      while(itK != std::end(kids))
      {
        double r = runif(gen);
       
        if(*itJ < r)
        {
          itK = kids.erase(itK);
          itJ = fitKids.erase(itJ);
        }
        else
        {
          ++itK;
          ++itJ;
        }
      }
      std::cout << " done." << std::endl << std::endl;
      
      // handles potential down-sampling (e.g. to test power)
      size_t ds = bpp::BppApplicationTools::getParameter<size_t>("downsample", params,
kids.size(), "", true, 4);

      if(ds > 0 && ds < kids.size())
      {
        std::cout << "Downsampling from " << kids.size() << " to " << ds << " children."
<< std::endl << std::endl;  
        std::shuffle(std::begin(kids), std::end(kids), gen);
        kids.erase(std::begin(kids) + ds, std::end(kids));
      }
            
      std::cout << "Generating " << numMeiosis << " theoretical siblings for each of " <<
kids.size() << " children..." << std::endl; 

      Tools pb;
      for(size_t i = 0; i < kids.size(); ++i)
      {
        kids[i]->genPotentialSiblings(numMeiosis, femaleRecMap, maleRecMap);
        double prog = static_cast<double>(i + 1) / static_cast<double>(kids.size());
        pb.print_pb(prog);
      }
      std::cout << std::endl << std::endl;
    }
    else // if !simKids, i.e., if we are reading children's genomes from VCF
    {
      // reads trio_ids.txt with the correspondence between kids and parents
      std::string triosFileName = bpp::BppApplicationTools::getAFilePath("trio_ids_file",
params, "none");
      std::ifstream trioIDsFile;
      trioIDsFile.open(triosFileName);

      // constructs children with pointers to their respective mothers and fathers
      std::vector<std::vector<std::string>> trioIDs(0, std::vector<std::string>(3));
      if(trioIDsFile.is_open())
      {
        std::vector<std::string> splitLine(0);  
        std::string focalTrioLine = "";
      
        getline(trioIDsFile, focalTrioLine); // skips header
      
        while(getline(trioIDsFile, focalTrioLine))
        {
          boost::split(splitLine, focalTrioLine, [=](char c) {return c == '\t';});
          trioIDs.push_back(splitLine);
        }
      }
      else
        throw bpp::Exception("Could not open " + triosFileName);
    
      size_t numKids = trioIDs.size(); // sibs and half-sibs are included
    
      std::cout << "Assembling trios..." << std::endl;
    
      Tools pb;
      std::string name, momName, dadName = "";
      kids.reserve(numKids);
      for(size_t i = 0; i < numKids; ++i)
      {
        name = trioIDs[i][0];
        momName = trioIDs[i][1];
        dadName = trioIDs[i][2];
      
        if(name != indvIds[2][i])
          throw bpp::Exception("Children ID's in VCF don't match trio_ids file!");
      
        // searches for mom and dad with matching names to forward parental ptrs to kids 
        auto itM = std::find_if(std::begin(moms), std::end(moms),
                                [&](std::shared_ptr<RealDiploid> m) 
                                { return m->getName() == momName; });
      
        auto itD = std::find_if(std::begin(dads), std::end(dads),
                                [&](std::shared_ptr<RealDiploid> d) 
                                { return d->getName() == dadName; });
            
        kids.emplace_back(std::make_shared<Child>(Child(tmp[2][i], *itM, *itD, name)));
      
        double prog = static_cast<double>(i + 1) / static_cast<double>(numKids);
        pb.print_pb(prog);
      }
      std::cout << " done." << std::endl << std::endl;
      
      
      std::cout << "Generating " << numMeiosis << " theoretical siblings for each of " <<
kids.size() << " children..." << std::endl; 

      pb.reset_pb();
      for(size_t i = 0; i < kids.size(); ++i)
      {
        kids[i]->genPotentialSiblings(numMeiosis, femaleRecMap, maleRecMap);
        double prog = static_cast<double>(i + 1) / static_cast<double>(kids.size());
        pb.print_pb(prog);
      }
      std::cout << std::endl << std::endl;
      
      
      // handling de novo mutations
      double mu = bpp::BppApplicationTools::getDoubleParameter("mu", params, 0.);
      double seqLen = bpp::BppApplicationTools::getDoubleParameter("seq_len", params, 0.);
      double lambda = 2. * mu * seqLen; // factor of two because of diploidy
    
      if(lambda > 0.) // if a valid mutation rate is provided in the options file
      {
        std::cout << "Adding de novo mutations with rate " << lambda << " per zygote...";
        std::cout.flush();
      
        for(size_t i = 0; i < numKids; ++i)
          kids[i]->addDeNovoMutations(lambda); // [~Pois(lambda)]
        
        std::cout << " done." << std::endl << std::endl;  
      }
      else // else we will remove de novo mutations from children
      {
        std::cout << "Building reference SNP panel..." << std::endl; std::cout.flush();
        // not parallelised because appendTrioSnps modifies refSnps
        // init (Haplotype typedef in HaploidGenome.h)
        Haplotype refSnps = kids[0]->getMother()->getHaplotypePair().first.getGenome(); // init
          
        pb.reset_pb();
        for(size_t i = 0; i < numKids; ++i)
        {
          auto mm = std::begin(kids[i]->getMother()->getHaplotypePair().first.getGenome());
          auto mmStop = std::end(kids[i]->getMother()->getHaplotypePair().first.getGenome());
   
          auto mp = std::begin(kids[i]->getMother()->getHaplotypePair().second.getGenome());
          auto pm = std::begin(kids[i]->getFather()->getHaplotypePair().first.getGenome());
          auto pp = std::begin(kids[i]->getFather()->getHaplotypePair().second.getGenome());
          auto egg = std::begin(kids[i]->getHaplotypePair().first.getGenome());
          auto sperm = std::begin(kids[i]->getHaplotypePair().second.getGenome());

          std::string chr = "";
          // since they all have the same chr's, only need to check 
          // the stop criterium for one of the iterators
          for(; mm != mmStop; ++mm, ++mp, ++pm, ++pp, ++egg, ++sperm)
          {
            chr = mm->first; // shared by all iterators

            refSnps.at(chr).insert(std::end(refSnps.at(chr)),
                                   std::begin(mm->second),
                                   std::end(mm->second));
            refSnps.at(chr).insert(std::end(refSnps.at(chr)),
                                   std::begin(mp->second),
                                   std::end(mp->second));
            refSnps.at(chr).insert(std::end(refSnps.at(chr)),
                                   std::begin(pm->second),
                                   std::end(pm->second)); 
            refSnps.at(chr).insert(std::end(refSnps.at(chr)),
                                   std::begin(pp->second),
                                   std::end(pp->second)); 
            refSnps.at(chr).insert(std::end(refSnps.at(chr)),
                                   std::begin(egg->second),
                                   std::end(egg->second));
            refSnps.at(chr).insert(std::end(refSnps.at(chr)),
                                   std::begin(sperm->second),
                                   std::end(sperm->second)); 
          }

          // deletes duplicates
          for(auto it = std::begin(refSnps); it != std::end(refSnps); ++it)
          {
            std::sort(std::begin(it->second), std::end(it->second));
            it->second.erase(std::unique(std::begin(it->second), std::end(it->second)),
                             std::end(it->second));
          }
          
          double prog = static_cast<double>(i + 1) / static_cast<double>(numKids);
          pb.print_pb(prog);
        }
        
        // moms
        std::vector<std::shared_ptr<RealDiploid>> momPtrs(0); 
        for(auto itKid = std::begin(kids); itKid != std::end(kids); ++itKid)
        {            
          // only includes unique mothers (since trios can be overlapping)
          bool isUnique = std::none_of(std::begin(momPtrs), std::end(momPtrs),
                          [&](std::shared_ptr<RealDiploid> p) 
                          { return p->getName() == (*itKid)->getMother()->getName(); });
      
          if(isUnique)
            momPtrs.push_back(std::make_shared<RealDiploid>(*((*itKid)->getMother())));
        }
        
        std::vector<std::shared_ptr<Diploid>> uniqueMoms(momPtrs.size());
        for(size_t i = 0; i < momPtrs.size(); i++)
        {
          std::shared_ptr<Diploid> m = std::dynamic_pointer_cast<Diploid>(momPtrs[i]);
          uniqueMoms[i] = m;
        }
        SnpInfo momsFreqs(std::move(Utils::fetch_snp_freqs(refSnps, uniqueMoms)));
    
    
        // dads
        std::vector<std::shared_ptr<RealDiploid>> dadPtrs(0);
        for(auto itKid = std::begin(kids); itKid != std::end(kids); ++itKid)
        {            
          // only includes unique fathers
          bool isUnique = std::none_of(std::begin(dadPtrs), std::end(dadPtrs),
                          [&](std::shared_ptr<RealDiploid> p) 
                          { return p->getName() == (*itKid)->getFather()->getName(); });
      
          if(isUnique)
            dadPtrs.push_back(std::make_shared<RealDiploid>(*((*itKid)->getFather())));
        }
        
        std::vector<std::shared_ptr<Diploid>> uniqueDads(dadPtrs.size());
        for(size_t i = 0; i < dadPtrs.size(); i++)
        {
          std::shared_ptr<Diploid> m = std::dynamic_pointer_cast<Diploid>(dadPtrs[i]);
          uniqueDads[i] = m;
        }
        SnpInfo dadsFreqs(std::move(Utils::fetch_snp_freqs(refSnps, uniqueDads))); 
          
        
        SnpInfo parentalFreqs; // sum of momsFreqs and dadsFreqs
        auto itChrMom = std::begin(momsFreqs);
        auto itChrDad = std::begin(dadsFreqs);
    
        // momsFreqs and dadsFreqs have same chr's and lengths (of course)
        for(; itChrMom != std::end(momsFreqs); ++itChrMom, ++itChrDad)
        {
          HashTable parentsTbl(std::move(itChrMom->second));
      
          for(size_t i = 0; i < parentsTbl.getTableSize(); ++i) 
          {
            HashNode* entry = parentsTbl.getTable()[i];
        
            while(entry != NULL)
            {
              size_t momKey = entry->getKey();
              size_t momVal = entry->getValue();
              size_t dadVal = itChrDad->second.getValue(momKey);

              entry->setValue(momVal + dadVal); // sums frequency
          
              // moves to next node in same bucket
              entry = entry->getNext(); 
            }
          }
          
          parentalFreqs.emplace(std::make_pair(itChrMom->first, std::move(parentsTbl)));
        }
        
        std::cout << "\n\nRemoving de novo SNPs from children..." << std::endl;
        // this method automatically updates counts of sites for individuals
        Utils::remove_denovo_snps(kids, parentalFreqs);
        std::cout << std::endl;
        
        std::cout << "Generating " << numMeiosis << " theoretical siblings for each of "
<< numKids << " children..." << std::endl; 

        pb.reset_pb();
        for(size_t i = 0; i < numKids; ++i)
        {
          kids[i]->genPotentialSiblings(numMeiosis, femaleRecMap, maleRecMap);
          double prog = static_cast<double>(i + 1) / static_cast<double>(numKids);
          pb.print_pb(prog);
        }
        std::cout << std::endl << std::endl;
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////
    
    /*
     * end of fake scope
     */
  } 
  //////////////////////////////////////////////////////////////////////////////////
  
  
   
  //////////////////////////////////////////////////////////////////////////////////
   size_t numSims = bpp::BppApplicationTools::getParameter<size_t>("num_sims", params,
18000000, "", true, 0);
    
  size_t numPilot = bpp::BppApplicationTools::getParameter<size_t>("num_pilot", params,
3000000, "", true, 0);
  
  size_t numAccepted = bpp::BppApplicationTools::getParameter<size_t>("num_accepted",
params, numPilot / 1000, "", true, 0);
  
  abc.computeObsStats(kids);
  
  std::cout << "Running pilot simulation of size " << numPilot << "..." << std::endl;
  abc.setPilot(true);
  abc.init(numPilot);
  abc.pilot(kids, numPilot, numAccepted); 

  std::cout << "Evaluating " << numSims << " parameter draws..." << std::endl;
  abc.setPilot(false);
  abc.simulate(kids, numSims); 
    
  std::cout << std::endl << std::endl;
  std::cout << "Writing (simulated) data to files..."; std::cout.flush();
  abc.writeSimParamsToFiles(numSims);
  abc.writeStatsToFiles(numSims);
  
  std::cout << " done." << std::endl << std::endl;
  /////////////////////////////////////////////////////////////////////////////////////


  tides.done();
  
  return (0);
}
