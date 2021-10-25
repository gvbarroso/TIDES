/*
  This file is a copy of BppApplication.cpp in the bpp-core libraries available in
  (https://github.com/BioPP/bpp-core/blob/master/src/Bpp/App/BppApplication.h)
  
  It's authors are Guillaume Deuchst, Julien Dutheil, Sylvain Gaillard and Francois Gindraud
  
  Copied on May 25 2020 for convenience reasons
*/

#include "BppApplication.h"
#include "BppAttributesTools.h"
#include "BppApplicationTools.h"

// From the STL:
#include <iostream>

using namespace bpp;

BppApplication::BppApplication(int argc, char* argv[], const std::string& name): appName_(name), params_(), timerStarted_(false)
{
  std::cout << "Parsing options:" << std::endl;  
  params_ = AttributesTools::parseOptions(argc, argv);
  BppApplicationTools::warningLevel = BppApplicationTools::getIntParameter("--warning", params_, 0, "", true, 3);
  bool noint = BppApplicationTools::getBooleanParameter("--noninteractive", params_, false, "", true, 3);
  BppApplicationTools::interactive = !noint;
}

void BppApplication::startTimer()
{
  BppApplicationTools::startTimer();
  timerStarted_ = true;
}

void BppApplication::done()
{
  std::cout << appName_ << "'s done. Bye." << std::endl;
  if (timerStarted_)
    BppApplicationTools::displayTime("Total execution time:");
}
