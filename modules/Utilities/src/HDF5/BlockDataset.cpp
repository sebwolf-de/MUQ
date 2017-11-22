#include "MUQ/Utilities/HDF5/BlockDataset.h"

using namespace muq::Utilities;


BlockDataset& BlockDataset::operator=(boost::any const& val) {

    AnyWriterMapType& map = *GetAnyWriterMap();
    map[val.type()](val, *this);
    
    return *this;
}

std::shared_ptr<BlockDataset::AnyWriterMapType> BlockDataset::GetAnyWriterMap() {
    
  static std::shared_ptr<AnyWriterMapType> map;

  if( !map )
    map = std::make_shared<AnyWriterMapType>();

  return map;
}

REGISTER_HDF5BLOCK_ANYTYPE(double)
REGISTER_HDF5BLOCK_ANYTYPE(float)
REGISTER_HDF5BLOCK_ANYTYPE(int)
REGISTER_HDF5BLOCK_ANYTYPE(unsigned)
