/*
  This file is a copy of BppApplication.h in the bpp-core libraries available in
  (https://github.com/BioPP/bpp-core/blob/master/src/Bpp/App/BppApplication.h)
  
  It's authors are Guillaume Deuchst, Julien Dutheil, Sylvain Gaillard and Francois Gindraud
  
  Copied on May 25 2020 for convenience reasons
*/

#ifndef _BPPAPPLICATION_H_
#define _BPPAPPLICATION_H_

#include "BppExceptions.h"

// From the STL:
#include <string>
#include <map>

namespace bpp
{

  class BppApplication
  {
    private:
      std::string appName_;
      mutable std::map<std::string, std::string> params_;
      bool timerStarted_;

    public:
      BppApplication(int argc, char* argv[], const std::string& name);

    public:
      void startTimer();
      void done();

      std::map<std::string, std::string>& getParams() { return params_; }

      const std::string& getParam(const std::string& name) const
      {
        if (params_.find(name) == params_.end()) throw Exception("BppApplication::getParam(). Parameter '" + name + "' not found.");
        return params_[name];
      }
      
      std::string& getParam(const std::string& name) { return params_[name]; }

  };

} //end of namespace bpp;

#endif // _BPPAPPLICATION_H_
