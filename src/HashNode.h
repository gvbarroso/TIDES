/*
 * Authors: Gustavo V. Barroso 
 * Created: 02/03/2020
 * Last modified: 11/07/2020
 *
 */


#ifndef _HASHNODE_
#define _HASHNODE_

#include <vector>
#include <string>
#include <algorithm>
#include <utility>

// hashing with chaining
class HashNode {
    
private:
    size_t key_; // SNP coordinate
    double value_; // SNP frequency
    // next bucket with the same prehash
    HashNode* next_;
    
public:
  HashNode(size_t key, double value):
  key_(key),
  value_(value),
  next_(NULL) 
  { }
  
  HashNode(const HashNode& other):
  key_(0),
  value_(0),
  next_(NULL) 
  {
    key_ = other.key_;
    value_ = other.value_;
    
    if(other.next_ != NULL)
      next_ = new HashNode(*other.next_); // deep copy
  }

  size_t getKey() const
  {
    return key_;
  }

  double getValue() const
  {
    return value_;
  }
    
  void setValue(double value)
  {
    HashNode::value_ = value;
  }

  HashNode* getNext() const
  {
    return next_;
  }

  void setNext(HashNode* next)
  {
    HashNode::next_ = next;
  }
  
};

#endif
