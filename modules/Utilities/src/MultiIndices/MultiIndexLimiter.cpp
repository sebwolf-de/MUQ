
#include <boost/serialization/export.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "MUQ/Utilities/MultiIndices/MultiIndexLimiter.h"
//#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"

/*
bool muq::Utilities::GeneralLimiter::IsFeasible(std::shared_ptr<MultiIndex> multi) const{
  return limitingSet->IsActive(multi);
};
*/

bool muq::Utilities::DimensionLimiter::IsFeasible(std::shared_ptr<MultiIndex> multi) const{
  
  for(auto pair = multi->GetNzBegin(); pair!=multi->GetNzEnd(); ++pair){
    if(((pair->first<lowerDim)||(pair->first>=lowerDim+length))&&(pair->second!=0))
      return false;
  }
  return true;
};

bool muq::Utilities::MaxOrderLimiter::IsFeasible(std::shared_ptr<MultiIndex> multi) const{
  
  if(maxOrders.size()==0){
    return (multi->Max() <= maxOrder);
  }else{
    assert(multi->GetLength()<=maxOrders.size());
    
    if(multi->Max() <= vectorMin)
      return true;
    
    for(auto nzIter = multi->GetNzBegin(); nzIter!=multi->GetNzEnd(); ++nzIter){
      if(nzIter->second>maxOrders(nzIter->first)){
        return false;
      }
    }
    return true;
  }
};


//BOOST_CLASS_EXPORT(muq::Utilities::DimensionLimiter)
//BOOST_CLASS_EXPORT(muq::Utilities::TotalOrderLimiter)
//BOOST_CLASS_EXPORT(muq::Utilities::MaxOrderLimiter)
//BOOST_CLASS_EXPORT(muq::Utilities::NoLimiter)
//BOOST_CLASS_EXPORT(muq::Utilities::AndLimiter)
//BOOST_CLASS_EXPORT(muq::Utilities::OrLimiter)
//BOOST_CLASS_EXPORT(muq::Utilities::XorLimiter)
