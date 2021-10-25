/*
  This file is a copy of NestedStringTokenizer.h in the bpp-core libraries available in
  (https://github.com/BioPP/bpp-core/blob/master/src/Bpp/Text/NestedStringTokenizer.h)
  
  It's authors are Guillaume Deuchst, Julien Dutheil, Sylvain Gaillard and Francois Gindraud
  
  Copied on May 25 2020 for convenience reasons
*/

#ifndef _NESTEDSTRINGTOKENIZER_H_
#define _NESTEDSTRINGTOKENIZER_H_

//From the STL:
#include <deque>
#include <string>

#include "BppStringTokenizer.h"
#include "BppExceptions.h"

namespace bpp
{

  /**
   * @brief An improved tokenizer for strings.
   *
   * Splits a string according to a given (set of) delimiter(s).
   * Delimiters in certains blocks ({}, [], etc) are ignored.
   */
  class NestedStringTokenizer:
    public StringTokenizer
  {
  public:
		
    /**
     * @brief Build a new StringTokenizer from a string.
     *
     * @param s          The string to parse.
     * @param open       Opening block.
     * @param end        Ending block.
     * @param delimiters Chars that must be considered as delimiters.
     * @param solid      If true, delimiters is considered as a single bloc delimiter.
     */
    NestedStringTokenizer(const std::string& s, const std::string& open, const std::string& end, const std::string& delimiters = " \t\n\f\r", bool solid = false);
		
    virtual ~NestedStringTokenizer() {}
	
  public:
		
    /**
     * @brief Get the next available token.
     * If no token is availbale, throw an Exception.
     *
     * @return The next token if there is one.
     */
    const std::string& nextToken();


    /**
     * @brief This function is not supported for nested tokenizers.
     *
     * @return An empty string.
     */
    std::string unparseRemainingTokens() const { return ""; }
  };

} //end of namespace bpp;

#endif	//_NESTEDSTRINGTOKENIZER_H_
