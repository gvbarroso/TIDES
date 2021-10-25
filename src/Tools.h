/*
 * Authors: Gustavo V. Barroso 
 * Created: 30/11/2019
 * Last modified: 13/07/2020
 *
 */


#ifndef _TOOLS_
#define _TOOLS_

#include <utility>
#include <chrono>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>

#include "BppExceptions.h"
#include "BppTextTools.h"
#include "GenomicRanges.h"
#include "HaploidGenome.h"

class Tools {

private:
  // timer
  std::chrono::time_point<std::chrono::high_resolution_clock> startTimePoint_;
  std::chrono::time_point<std::chrono::high_resolution_clock> endTimePoint_;

  //progress bar
  std::string firstPartOfpBar;
  std::string lastPartOfpBar;
  std::string pBarFiller;
  std::string pBarUpdater;
  
  size_t amountOfFiller_;
  size_t pBarLength_;
  size_t currUpdateVal_; 
  double neededProgress_; 

public:
  Tools():
  startTimePoint_(),
  endTimePoint_(),
  firstPartOfpBar("["),
  lastPartOfpBar("]"),
  pBarFiller("|"),
  pBarUpdater("/-\\|"),
  amountOfFiller_(-1),
  pBarLength_(50),
  currUpdateVal_(0),
  neededProgress_(1.)
  { }

  void print_pb(double newProgress);

  void reset_pb()
  {
    currUpdateVal_ = 0;
    amountOfFiller_ = 0;
  }

  void stop_timer(double mult, const std::string& task, const std::string& unit);

  void start_timer()
  {
    startTimePoint_ = std::chrono::high_resolution_clock::now();
  }
  
};

#endif
