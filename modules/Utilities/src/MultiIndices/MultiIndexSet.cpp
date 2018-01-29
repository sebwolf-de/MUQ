#include "MUQ/Utilities/multiIndex/MultiIndexSet.h"

#include <algorithm>

using namespace muq::Utilities;
using namespace std;

std::shared_ptr<MultiIndexSet> muq::Utilities::operator+=( std::shared_ptr<MultiIndexSet> x, const std::shared_ptr<MultiIndexSet>& y)
{
  (*x)+=(*y);
  return x;
}

std::shared_ptr<MultiIndexSet> muq::Utilities::operator+=( std::shared_ptr<MultiIndexSet> x, const std::shared_ptr<MultiIndex>& y)
{
  (*x)+=y;
  return x;
}

MultiIndexSet::MultiIndexSet(const unsigned dimIn,
                             shared_ptr<MultiIndexLimiter> limiterIn,
                             std::shared_ptr<MultiIndexPool> poolIn ) : maxOrders(Eigen::VectorXu::Zero(dimIn)),
                                                                        dim(dimIn),
                                                                        setID(poolIn->GetNextId()),
                                                                        pool(poolIn),
                                                                        limiter(limiterIn)
{
  assert(pool!=nullptr);
};

void MultiIndexSet::SetLimiter(std::shared_ptr<MultiIndexLimiter> limiterIn){
  
  // first, make sure no active terms in the set currently break the new limiter.  If a term is inactive, remove all edges tied to it
  for(int localInd=0; localInd<local2global.size(); ++localInd)
  {
    if(IsActive(localInd)){
      assert(limiterIn->IsFeasible(pool->at(local2global.at(localInd))));
    }else{
      
      if(!limiterIn->IsFeasible(pool->at(local2global.at(localInd)))){
        for(int inNode : inEdges[localInd])
          outEdges[inNode].erase(localInd);
        inEdges[localInd].clear();
      }
    }
  }
  
  // copy the limiter
  limiter = limiterIn;
}

std::shared_ptr<MultiIndexSet> MultiIndexSet::CloneExisting(std::shared_ptr<MultiIndexSet> original)
{
  auto output = make_shared<MultiIndexSet>(original->dim, original->limiter, original->pool);
  output->active2local  = original->active2local;
  output->local2active  = original->local2active;
  output->local2global = original->local2global;
  output->outEdges    = original->outEdges;
  output->inEdges     = original->inEdges;
  output->multi2local = original->multi2local;
  output->maxOrders   = original->maxOrders;
  //output->numInds     = original->numInds;
  
  // now, make sure the pool knows that all these are active
  for(int i=0; i<output->active2local.size(); ++i)
    output->pool->SetActive(output->local2global[output->active2local[i]],output->setID,i);
  
  return output;
}

Eigen::MatrixXu MultiIndexSet::GetAllMultiIndices() const
{
  
  Eigen::MatrixXu output(active2local.size(), dim);
  for(int i=0; i<active2local.size(); ++i){
    auto multi = pool->at(local2global[active2local[i]]);
    int multiDim = multi->GetDimension();
    output.row(i).head(multiDim) = multi->GetMulti();
  }
  
  return output;
}

Eigen::RowVectorXu MultiIndexSet::IndexToMulti(unsigned int const i) const{
  assert(i<active2local.size());
  return pool->at(local2global[active2local[i]])->GetMulti();
}

int MultiIndexSet::MultiToIndex(const std::shared_ptr<MultiIndex> input) const{
 
  auto localIter = multi2local.find(input);
  
  if(localIter!=multi2local.end()){
    return local2active[localIter->second];
  }else{
    return -1;
  }
}

int MultiIndexSet::MultiToIndex(Eigen::RowVectorXu const& multiIn) const{
  auto multiPtr = make_shared<MultiIndex>(multiIn);
  return MultiToIndex(multiPtr);
}


//void MultiIndexSet::UpdateMaxOrders(){
//  for(auto ind : active2local){
//    for(auto pair : allMultis[ind]->nzInds)
//      maxOrders(pair.first) = std::max<unsigned>(maxOrders(pair.first),pair.second);
//  }
//}

//void MultiIndexSet::UpdateMultiMap(){
//  
//  // clear the map
//  multi2local = std::map<std::shared_ptr<MultiIndex>, unsigned int, MultiPtrComp>();
//  
//  // loop through all the multis and update the map with negative indices
//  for(int i=0; i<allMultis.size(); ++i)
//    multi2local[allMultis[i]] = i;   
//}

int MultiIndexSet::AddActive(std::shared_ptr<MultiIndex>                                        newNode,
                             std::map<std::shared_ptr<MultiIndex>, unsigned int, MultiPtrComp>::iterator iter)
{
  if(iter!=multi2local.end()){
  
    int localInd = iter->second;
    
    // the node already exists, all we have to do is make it active
    Activate(localInd);
      
    return local2active[localInd];
    
  }else if(limiter->IsFeasible(newNode)){
  
    // the multiindex is not active or inactive, add it to the list
    int newGlobalInd = pool->AddMulti(newNode);
    
    local2global.push_back(newGlobalInd);
    int newLocalInd = local2global.size()-1;
    
    inEdges.push_back(std::set<int>());
    outEdges.push_back(std::set<int>());
    
    // store it as an active node
    local2active.push_back(-1);
    multi2local[newNode] = newLocalInd;
    
    Activate(newLocalInd);
    
    return local2active[newLocalInd];
    
  }else{
    // the node was not feasible, so return a negative index
    return -1;
  }

}

int MultiIndexSet::AddActive(std::shared_ptr<MultiIndex> newNode)
{
  // first, check to see if the index is in the list of nodes  
  auto iter = multi2local.find(newNode);
  return AddActive(newNode,iter);
}

int MultiIndexSet::AddInactive(std::shared_ptr<MultiIndex>                                        newNode,
                                    std::map<std::shared_ptr<MultiIndex>, unsigned int, MultiPtrComp>::iterator iter){

  if(iter==multi2local.end()){
    
    if(limiter->IsFeasible(newNode)){
    
      // the multiindex is not stored, add it to the list
      int newGlobalInd = pool->AddMulti(newNode);
      
      local2global.push_back(newGlobalInd);
      int newLocalInd = local2global.size()-1;
      
      inEdges.push_back(std::set<int>());
      outEdges.push_back(std::set<int>());
    
      // store it as an inactive node
      local2active.push_back(-1);
      multi2local[newNode] = newLocalInd;
    
      AddForwardNeighbors(newLocalInd,false);
      AddBackwardNeighbors(newLocalInd,false);
    
      return newLocalInd;
      
    }else{
      // the node was not feasible, so return a negative index
      return -1;
    }
    
  }else{
    return iter->second;
  }
  
}

int MultiIndexSet::AddInactive(std::shared_ptr<MultiIndex> newNode){
  
  auto iter = multi2local.find(newNode);
  return AddInactive(newNode,iter);
  
}

bool MultiIndexSet::IsActive(std::shared_ptr<MultiIndex> multiIndex) const
{
  auto iter = multi2local.find(multiIndex);
  if(iter!=multi2local.end()){
    return IsActive(iter->second);
  }else{
    return false;
  }
}

bool MultiIndexSet::IsActive(unsigned int localIndex) const
{
  if(local2active[localIndex]!=-1){
    return true;
  }else{
    return false;
  }
}

bool MultiIndexSet::IsAdmissible(unsigned int localIndex) const
{
  auto multi = pool->GetMulti(local2global.at(localIndex));
  if(!limiter->IsFeasible(multi))
    return false;
  
  if(IsActive(localIndex))
    return true;

  // count the number of input edges that are coming from active indices
  int numAdmiss = 0;
  for(int inNode : inEdges.at(localIndex)){
    if(IsActive(inNode))
      numAdmiss++;
  }
  
  if(numAdmiss==multi->nzInds.size()){
    return true;
  }else{
    return false;
  }
}

bool MultiIndexSet::IsAdmissible(std::shared_ptr<MultiIndex> multiIndex) const
{
  auto iter = multi2local.find(multiIndex);
  if(iter==multi2local.end()){
    return false;
  }else{
    return IsAdmissible(iter->second);
  }
}

bool MultiIndexSet::IsAdmissible(Eigen::RowVectorXu const& multiIndex) const
{
  return IsAdmissible(make_shared<MultiIndex>(multiIndex));
}

bool MultiIndexSet::IsExpandable(unsigned int activeIndex) const
{
  assert(activeIndex<active2local.size());

  // an index is expandable when at least one forward neighbor is admissible but not active (i.e. outedge)
  
  // loop through the outgoing edges for this node
  for(int nextInd : outEdges[active2local[activeIndex]]){
    if(!IsActive(nextInd)&&IsAdmissible(nextInd))
      return true;
  }
  return false;
}

void MultiIndexSet::Activate(int localIndex)
{
  assert(localIndex<local2global.size());
  
  // the index is already in the global set, if the value is non-negative, it is also active and we don't need to do anything
  if(local2active[localIndex]<0)
  {
    auto newNode = pool->GetMulti(local2global[localIndex]);
    
    // now add the index to the active set
    active2local.push_back(localIndex);
    int newActiveInd = active2local.size()-1;
    local2active[localIndex] = newActiveInd;
    pool->SetActive(local2global[localIndex],setID,newActiveInd);
    
    // update the maximum order
    for(auto pair : newNode->nzInds)
      maxOrders(pair.first) = std::max<unsigned>(maxOrders(pair.first),pair.second);
    
    AddForwardNeighbors(localIndex,true);
    AddBackwardNeighbors(localIndex,true);
  }
}

void MultiIndexSet::Activate(std::shared_ptr<MultiIndex> multiIndex)
{
  auto iter = multi2local.find(multiIndex);
  
  assert(iter!=multi2local.end());
  assert(IsAdmissible(iter->second));
  
  Activate(iter->second);
}

void MultiIndexSet::AddForwardNeighbors(unsigned int localIndex, bool addInactive)
{
  Eigen::RowVectorXu base = pool->at(local2global[localIndex])->GetMulti();
  shared_ptr<MultiIndex> newNode;
  for(unsigned int i=0; i<base.size(); ++i)
  {
    base(i)++;
    
    newNode = make_shared<MultiIndex>(base);
    if(limiter->IsFeasible(newNode)){
      auto iter = multi2local.find(newNode);
    
      if(iter!=multi2local.end()){
        inEdges.at(iter->second).insert(localIndex);
        outEdges.at(localIndex).insert(iter->second);
      }else if(addInactive){
        AddInactive(newNode,iter);
      }
    }
    base(i)--;
  }
}


Eigen::MatrixXu  MultiIndexSet::GetAdmissibleForwardNeighbors(unsigned int activeIndex)
{
  assert(activeIndex<active2local.size());
  unsigned int localInd = active2local.at(activeIndex);
  
  vector<shared_ptr<MultiIndex>> output;
  for( auto neighbor : outEdges[localInd])
  {
    if(IsAdmissible(neighbor))
      output.push_back(pool->GetMulti(local2global[neighbor]));
  }
  
  if(output.size()>0){
  
    Eigen::MatrixXu outputMat(output.size(),dim);
    for(int row =0; row<output.size();++row)
      outputMat.row(row).head(output[row]->GetDimension()) = output[row]->GetMulti();
    
    return outputMat;
    
  }else{
    return Eigen::MatrixXu();
  }
}

void MultiIndexSet::AddBackwardNeighbors(unsigned int localIndex, bool addInactive)
{
  Eigen::RowVectorXu base = pool->at(local2global[localIndex])->GetMulti();
  shared_ptr<MultiIndex> newNode;
  for(unsigned int i=0; i<base.size(); ++i)
  {
    if(base(i)>0){
    
      base(i)--;
      newNode = make_shared<MultiIndex>(base);
      if(limiter->IsFeasible(newNode)){
        auto iter = multi2local.find(newNode);
    
        if(iter!=multi2local.end()){
          outEdges.at(iter->second).insert(localIndex);
          inEdges.at(localIndex).insert(iter->second);
        }else if(addInactive){
          AddInactive(newNode,iter);
        }
      }
      base(i)++;
    }
  }
}

Eigen::VectorXu MultiIndexSet::Expand(unsigned int activeIndex)
{
  assert(activeIndex<active2local.size());

  vector<unsigned int> newIndices;
  int localIndex = active2local[activeIndex];
  
  // loop through the forward neighbors of this index
  std::set<int> tempSet = outEdges.at(localIndex);
  for(int neighbor : tempSet)
  {
    if(IsAdmissible(neighbor)&&(!IsActive(neighbor))){
      Activate(neighbor);
      newIndices.push_back(local2active[neighbor]);
    }
  }
  
  // return the vector of newly activated indices
  return Eigen::Map<Eigen::VectorXu>(&newIndices[0],newIndices.size());
}

Eigen::VectorXu MultiIndexSet::ForciblyExpand(unsigned int const activeIndex)
{
  assert(activeIndex<active2local.size());
  
  vector<unsigned int> newIndices;
  int localIndex = active2local.at(activeIndex);
  
  // loop through the forward neighbors of this index
  std::set<int> tempSet = outEdges.at(localIndex);
  for(int neighbor : tempSet)
    ForciblyActivate(neighbor,newIndices);
  
  // return the vector of newly activated indices
  return Eigen::Map<Eigen::VectorXu>(&newIndices[0],newIndices.size());
  
}

void MultiIndexSet::ForciblyActivate(int localIndex, std::vector<unsigned int> &newIndices){


  if(!IsActive(localIndex)){
  
    // make the node active and add inactive neighbors if necessary, this also updates the edges and enables the loop below
    Activate(localIndex); 
    newIndices.push_back(local2active.at(localIndex));
    
    // now, fill in all of the previous neighbors
    std::set<int> tempSet = inEdges.at(localIndex);
    for(int ind : tempSet)
      ForciblyActivate(ind,newIndices);
      
  }
}

Eigen::VectorXu MultiIndexSet::ForciblyActivate(std::shared_ptr<MultiIndex> multiIndex){

  assert(limiter->IsFeasible(multiIndex));
  
  auto iter = multi2local.find(multiIndex);
  vector<unsigned int> newIndices;
  
  // if we found the multiindex and it is active, there is nothing to do
  if(iter!=multi2local.end()){
    ForciblyActivate(iter->second,newIndices);
  }else{
    // Add the new index as an active node
    int newLocalInd = AddInactive(multiIndex,iter);
    ForciblyActivate(newLocalInd,newIndices);
  }
  
  return Eigen::Map<Eigen::VectorXu>(&newIndices[0],newIndices.size());
}

MultiIndexSet& MultiIndexSet::operator+=(const MultiIndexSet& rhs)
{
  Union(rhs);
  return *this;
}

int MultiIndexSet::Union(const MultiIndexSet& rhs)
{
  int oldTerms = size();
  
  for(int i = 0; i < rhs.local2active.size(); ++i) {
    
    auto newMulti = rhs.pool->at(rhs.local2global.at(i));
    if(limiter->IsFeasible(newMulti)){
      if(rhs.local2active[i]<0){
        AddInactive(newMulti);
      }else{
        AddActive(newMulti);
      }
    }
  }
  
  return size() - oldTerms;
}

MultiIndexSet& MultiIndexSet::operator+=(std::shared_ptr<MultiIndex> rhs)
{
  AddActive(rhs);
  return *this;
}

MultiIndexSet& MultiIndexSet::operator+=(Eigen::RowVectorXu const& multiIndex){
  AddActive(make_shared<MultiIndex>(multiIndex));
  return *this;
}

