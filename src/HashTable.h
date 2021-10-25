/*
 * Authors: Gustavo V. Barroso 
 * Created: 02/03/2020
 * Last modified: 28/07/2020
 *
 */

 
#ifndef _HASHTABLE_
#define _HASHTABLE_

#include <random>
#include <iostream>

#include "HashNode.h"


class HashTable {
    
private:
  HashNode** table_;
    
  size_t tableSize_;
  size_t largePrime_; 
  
  size_t a_;
  size_t b_;
    
public:
  HashTable(size_t tableSize, size_t universeSize): 
  table_(),
  tableSize_(tableSize),
  largePrime_(0),
  a_(0),
  b_(0)
  {
    // constructs zero initialized hash table of size
    table_ = new HashNode* [tableSize_]();
    largePrime_ = findLargePrime_(universeSize);
    
    std::mt19937 gen;  
    std::random_device rd;
    gen.seed(rd());  
    
    std::uniform_int_distribution<size_t> dist(0, largePrime_ - 1); 
  
    a_ = dist(gen);
    b_ = dist(gen);    
  }
  
  HashTable(): 
  table_(nullptr),
  tableSize_(0),
  largePrime_(0),
  a_(0),
  b_(0)
  { } 
  
  HashTable(HashTable&& other):
  table_(nullptr),
  tableSize_(0),
  largePrime_(0),
  a_(0),
  b_(0)
  {
    table_ = other.table_;
    tableSize_ = other.tableSize_;
    largePrime_ = other.largePrime_;
    a_ = other.a_;
    b_ = other.b_;
    
    other.table_ = nullptr;
    other.tableSize_ = 0;
    other.largePrime_ = 0;
    other.a_ = 0;
    other.b_ = 0;
  }
  
  HashTable(const HashTable& other):
  table_(nullptr),
  tableSize_(0),
  largePrime_(0),
  a_(0),
  b_(0)
  {
    std::cout << "copy constructor of HashTable: this will likely blow up!" << std::endl;
    // deep copy 
    table_ = new HashNode* [other.tableSize_]; 
    for(size_t i = 0; i < other.tableSize_; ++i)
    { 
      if(other.table_[i] != NULL)
        table_[i] = new HashNode(*(other.table_[i]));
    }
    
    tableSize_ = other.tableSize_;
    largePrime_ = other.largePrime_;
    a_ = other.a_;
    b_ = other.b_;
  }
  
  HashTable& operator=(const HashTable& other)
  {
    std::cout << "assignment ( = ) of HashTable: this will likely blow up!" << std::endl;
    if(this == &other)
      return *this;
    else
    {
      delete [] table_;
      // deep copy 
      table_ = new HashNode* [other.tableSize_]; 
      for(size_t i = 0; i < other.tableSize_; ++i)
      { 
        if(other.table_[i] != NULL)
          table_[i] = new HashNode(*(other.table_[i]));
      }
      
      tableSize_ = other.tableSize_;
      largePrime_ = other.largePrime_;
      a_ = other.a_;
      b_ = other.b_;
    }
    
    return *this;
  }
  
  HashTable& operator=(HashTable& other)
  {
    if(this == &other)
      return *this;
    else
    {
      delete [] table_;
      // deep copy 
      table_ = new HashNode* [other.tableSize_]; 
      for(size_t i = 0; i < other.tableSize_; ++i)
      { 
        if(other.table_[i] != NULL)
          table_[i] = new HashNode(*(other.table_[i]));
      }
    
      tableSize_ = other.tableSize_;
      largePrime_ = other.largePrime_;
      a_ = other.a_;
      b_ = other.b_;
    }
    
    return *this;
  }
  
  ~HashTable() 
  {
    // destroys all buckets one by one
    for(size_t i = 0; i < tableSize_; ++i)
    {
      HashNode* entry = table_[i];
     
      while(entry != NULL) 
      {
        HashNode* prev = entry;
        entry = entry->getNext();
      
        delete prev;
      }
      
      table_[i] = NULL;
    }

    delete [] table_;
  }

public:
  void reset(); // re-set all values (HashNode::value_) to 0
  
  void setAllValues(double val); // set all values to val

  // initializes all values (HashNode::value_) to 0
  void init(const std::vector<size_t>& chrSnps); 
  
  void init(const std::vector<size_t>& chrSnps, double val); // init all values to val
  
  void setValue(size_t key, double val); 
  
  void increment(size_t key, double val);
  
  void insert(size_t key, double val);

  void remove(size_t key);
  
  double getValue(size_t key);
  
  size_t getTableSize()
  {
    return tableSize_;
  }
  
  HashNode**& getTable()
  {
    return table_;
  }
  
private:
  size_t preHash_(size_t key);
  
  size_t findLargePrime_(size_t tableSize);

};

#endif
