/*
 * Authors: Gustavo V. Barroso 
 * Created: 30/11/2019
 * Last modified: 28/04/2019
 *
 */


#include "Tools.h"

void Tools::print_pb(double newProgress) 
{
  // updates
  amountOfFiller_ = static_cast<size_t>((newProgress / neededProgress_) *
static_cast<double>(pBarLength_));

  currUpdateVal_ %= pBarUpdater.length();
  std::cout << "\r" << firstPartOfpBar; 
    
  for(size_t a = 0; a < amountOfFiller_; a++)
    std::cout << pBarFiller;

  std::cout << pBarUpdater[currUpdateVal_];

  size_t percentage = static_cast<size_t>(100 * (newProgress / neededProgress_));
  if(percentage > 100)
    percentage = 100;
  
  std::cout << lastPartOfpBar << " (" << percentage << "%)" << std::flush;

  currUpdateVal_ += 1.;
}

void Tools::stop_timer(double mult, const std::string& task, const std::string& unit) 
{
  endTimePoint_ = std::chrono::high_resolution_clock::now();

  auto start = std::chrono::time_point_cast<std::chrono::microseconds>
(startTimePoint_).time_since_epoch().count();

  auto end = std::chrono::time_point_cast<std::chrono::microseconds>
(endTimePoint_).time_since_epoch().count();

  auto duration = end - start;
  double conv = duration / mult;
    
  std::cout << "Duration of " << task << " = " << conv << " (" << unit << ")\n";
}

