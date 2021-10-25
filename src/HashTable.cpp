/*
 * Authors: Gustavo V. Barroso 
 * Created: 02/03/2020
 * Last modified: 28/07/2020
 *
 * Modified from:
 * https://medium.com/@aozturk/simple-hash-map-hash-table-implementation-in-c-931965904250
 */


#include "HashTable.h"


void HashTable::reset()
{
  for(size_t i = 0; i < tableSize_; ++i)
  {
    HashNode* entry = table_[i];
    
    while(entry != NULL) 
    {
      entry->setValue(0.);
      entry = entry->getNext();
    }
  }
}

void HashTable::setAllValues(double val) 
{
  for(size_t i = 0; i < tableSize_; ++i)
  {
    HashNode* entry = table_[i];
    
    while(entry != NULL) 
    {
      entry->setValue(val);
      entry = entry->getNext();
    }
  }
}

void HashTable::init(const std::vector<size_t>& chrSnps)
{
  for(auto itSnp = std::begin(chrSnps); itSnp != std::end(chrSnps); ++itSnp)
    insert(*itSnp, 0.);
}

void HashTable::init(const std::vector<size_t>& chrSnps, double val)
{
  for(auto itSnp = std::begin(chrSnps); itSnp != std::end(chrSnps); ++itSnp)
    insert(*itSnp, val);
}

void HashTable::increment(size_t key, double val) 
{
  size_t hashValue = preHash_(key);
  
  HashNode* entry = table_[hashValue];

  while(entry->getKey() != key)
    entry = entry->getNext();

  double newVal = entry->getValue() + val;
  entry->setValue(newVal);
}

void HashTable::setValue(size_t key, double val) 
{
  size_t hashValue = preHash_(key);
  
  HashNode* entry = table_[hashValue];

  while(entry->getKey() != key)
    entry = entry->getNext();

  entry->setValue(val);
}

void HashTable::insert(size_t key, double value) 
{
  size_t hashValue = preHash_(key);
  
  HashNode* prev = NULL;
  HashNode* entry = table_[hashValue];

  while(entry != NULL && entry->getKey() != key)
  {
    prev = entry;
    entry = entry->getNext();
  }

  if(entry == NULL)
  {
    entry = new HashNode(key, value);
        
    if(prev == NULL) // insert as first bucket 
      table_[hashValue] = entry;
    else
      prev->setNext(entry);
  }
  else // just update the value
    entry->setValue(value);
}

void HashTable::remove(size_t key) 
{
  size_t hashValue = preHash_(key);
  
  HashNode* prev = NULL;
  HashNode* entry = table_[hashValue];

  while(entry != NULL && entry->getKey() != key) 
  {
    prev = entry;
    entry = entry->getNext();
  }

  if(entry == NULL) // key not found
    return;
  else
  {
    if(prev == NULL) // remove first bucket of the list
      table_[hashValue] = entry->getNext();
    else
      prev->setNext(entry->getNext());
      
    delete entry;
  }
}

double HashTable::getValue(size_t key)
{
  size_t hashValue = preHash_(key);
  
  HashNode* entry = table_[hashValue];
  
  while(entry->getKey() != key)
    entry = entry->getNext();

  return entry->getValue();
}

size_t HashTable::findLargePrime_(size_t number)
{
  auto isPrime = [=](size_t n)  
  {  
    // corner cases  
    if(n <= 1)
      return false;  
    if(n <= 3)
      return true;  
    
    // this is checked so that we can skip   
    // middle five numbers in below loop  
    if(n % 2 == 0 || n % 3 == 0)
      return false;  
    
    for(size_t i = 5; i * i <= n; i = i + 6)  
      if(n % i == 0 || n % (i + 2) == 0)  
        return false;  
    
    return true;  
  }; 

  size_t prime = number; 
  bool found = false; 
  
  while(!found) 
  { 
    ++prime; 
  
    if(isPrime(prime)) 
      found = true; 
  }
  
  return prime; 
}

size_t HashTable::preHash_(size_t key)
{
  return (((a_ * key) + b_) % largePrime_) % tableSize_;
}
