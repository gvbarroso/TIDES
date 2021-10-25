/*
  This file is a copy of Exception.cpp in the bpp-core libraries available in
  (https://github.com/BioPP/bpp-core/blob/master/src/Bpp/Exceptions.cpp)
  
  It's authors are Guillaume Deuchst, Julien Dutheil, Sylvain Gaillard and Francois Gindraud
  
  Copied on May 25 2020 for convenience reasons
*/

#include <string>
#include <utility>

#include "BppExceptions.h"

namespace bpp
{
  Exception::Exception(std::string text)
    : message_(std::move(text))
  {
  }
  const char* Exception::what() const noexcept { return message_.c_str(); }
  const std::string& Exception::message() const noexcept { return message_; }

  IOException::IOException(std::string text)
    : Exception(std::move(text))
  {
  }

  NullPointerException::NullPointerException(std::string text)
    : Exception(std::move(text))
  {
  }

  ZeroDivisionException::ZeroDivisionException(std::string text)
    : Exception(std::move(text))
  {
  }

  BadIntegerException::BadIntegerException(std::string text, int badInt)
    : Exception(std::move(text) + '(' + std::to_string(badInt) + ')')
    , badInt_(badInt)
  {
  }
  int BadIntegerException::getBadInteger() const { return badInt_; }

  BadNumberException::BadNumberException(std::string text, double badNumber)
    : Exception(std::move(text) + '(' + std::to_string(badNumber) + ')')
    , badNumber_(badNumber)
  {
  }
  double BadNumberException::getBadNumber() const { return badNumber_; }

  NumberFormatException::NumberFormatException(std::string text, std::string badNumber)
    : Exception(std::move(text) + '(' + std::move(badNumber) + ')')
    , badNumber_(badNumber)
  {
  }
  const std::string& NumberFormatException::getBadNumber() const { return badNumber_; }

  IndexOutOfBoundsException::IndexOutOfBoundsException(std::string text, std::size_t badInt, std::size_t lowerBound,
                                                       std::size_t upperBound)
    : Exception("out of [" + std::to_string(lowerBound) + ", " + std::to_string(upperBound) + "])" + std::move(text))
    , badIndex_(badInt)
    , bounds_{{lowerBound, upperBound}}
  {
  }
  const std::array<std::size_t, 2>& IndexOutOfBoundsException::getBounds() const { return bounds_; }
  std::size_t IndexOutOfBoundsException::getBadIndex() const { return badIndex_; }

  BadSizeException::BadSizeException(std::string text, std::size_t badSize, std::size_t correctSize)
    : Exception("Incorrect size " + std::to_string(badSize) + ", expected " + std::to_string(correctSize) + ". " +
                std::move(text))
    , badSize_(badSize)
    , correctSize_(correctSize)
  {
  }
  std::size_t BadSizeException::getBadSize() const { return badSize_; }
  std::size_t BadSizeException::getCorrectSize() const { return correctSize_; }

  OutOfRangeException::OutOfRangeException(std::string text, double badValue, double lowerBound, double upperBound)
    : Exception(std::to_string(badValue) + " out of [" + std::to_string(lowerBound) + ", " +
                std::to_string(upperBound) + "])" + std::move(text))
    , badValue_(badValue)
    , bounds_{{lowerBound, upperBound}}
  {
  }
  double OutOfRangeException::getBadValue() const { return badValue_; }
  double OutOfRangeException::getLowerBound() const { return bounds_[0]; }
  double OutOfRangeException::getUpperBound() const { return bounds_[1]; }

  NotImplementedException::NotImplementedException(std::string text)
    : Exception(std::move(text))
  {
  }
} // namespace bpp
