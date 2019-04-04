#include "MUQ/Utilities/MultiIndices/MultiIndex.h"

#include <stdexcept>

#include <iostream>

using namespace muq::Utilities;


MultiIndex::MultiIndex(unsigned lengthIn) : length(lengthIn),
                                            maxValue(0),
                                            totalOrder(0)
{}


MultiIndex::MultiIndex(Eigen::RowVectorXi const& indIn) : MultiIndex(indIn.size())
{
  maxValue = 0;
  totalOrder = 0;

  for(int i=0; i<indIn.size(); ++i){
    if( indIn[i]>0 ){
      nzInds[i] = indIn[i];
      maxValue = std::max<int>(maxValue, indIn[i]);
      totalOrder += indIn[i];
    }
  }
}

MultiIndex::MultiIndex(std::initializer_list<unsigned> const& indIn) : MultiIndex(indIn.size())
{
  maxValue = 0;
  totalOrder = 0;

  unsigned i = 0;
  for(auto it = indIn.begin(); it != indIn.end(); ++it){
    if( *it > 0 ){
      nzInds[i] = *it;

      maxValue = std::max<int>(maxValue, *it);
      totalOrder += *it;

      i++;
    }
  }
}


Eigen::RowVectorXi MultiIndex::GetVector() const
{
  Eigen::RowVectorXi output = Eigen::RowVectorXi::Zero(length);

  for(auto it = nzInds.begin(); it!=nzInds.end(); ++it)
    output(it->first) = it->second;

  return output;
}

bool MultiIndex::SetValue(unsigned ind, unsigned val)
{
  if(ind>length){
    throw std::out_of_range("Tried to set the value of index " + std::to_string(ind) + " on an multiindex with only " + std::to_string(length) + " components.");
  }else{

    auto it = nzInds.find(ind);
    if(it != nzInds.end()){
      it->second = val;
    }else{
      nzInds[ind] = val;
    }


    // Update the total and maximum order values after updating multi index
    totalOrder = 0;
    maxValue = 0;

    for (const auto& value : nzInds){
      totalOrder += value.second;
      maxValue = std::max(maxValue, value.second);
    }

    return it!=nzInds.end();
  }
}

unsigned int MultiIndex::NumNz() const
{
    unsigned int numNz = 0;
    for(auto& part : nzInds)
      numNz += int(part.second > 0);

    return numNz;
}


unsigned MultiIndex::MultiIndex::GetValue(unsigned ind) const
{
  if(ind>length){
    throw std::out_of_range("Tried to access index " + std::to_string(ind) + " of a multiindex with only " + std::to_string(length) + " components.");
  }else{
    auto searchIter = nzInds.find(ind);
    if(searchIter != nzInds.end()){
      return searchIter->second;
    }else{
      return 0;
    }
  }
}


void MultiIndex::SetLength(unsigned newLength)
{
  if(newLength > length){
    length = newLength;
  }else{

    auto it = nzInds.begin();
    while(it!= nzInds.end()){
      if (it->first >= newLength) {
	      it = nzInds.erase(it);
      } else {
	      it++;
      }
    }

    // Update the stored summaries
    length = newLength;
    maxValue = 0;
    totalOrder = 0;
    for(auto it = nzInds.begin(); it!=nzInds.end(); ++it){
      maxValue = std::max(maxValue, it->second);
      totalOrder += it->second;
    }

  }
}

bool MultiIndex::operator!=(const MultiIndex &b){

  if( (b.length != length) || (b.maxValue != maxValue) || (b.totalOrder != totalOrder))
    return true;

  if(b.nzInds.size() != nzInds.size())
    return true;

  // Loop through the nonzero indices
  auto bit = b.nzInds.begin();
  auto it = nzInds.begin();
  for(; it!=nzInds.end(); ++it){
    if(it->first != bit->first)
      return true;
    if(it->second != bit->second)
      return true;
  }

  return false;
}

bool MultiIndex::operator==(const MultiIndex &b){
  return !( *this != b);
}

bool MultiIndex::operator<(const MultiIndex &b){

  if(totalOrder<b.totalOrder){
    return true;
  }else if(totalOrder>b.totalOrder){
      return false;
  }else if(maxValue<b.maxValue){
    return true;
  }else if(maxValue>b.maxValue){
    return false;
  }else{

    for(int i=0; i<std::min<unsigned>(length, b.length); ++i){
      if(GetValue(i)<b.GetValue(i)){
        return true;
      }else if(GetValue(i)>b.GetValue(i)){
        return false;
      }
    }

    // it should never get to this point unless the multiindices are equal
    return false;
  }

}
