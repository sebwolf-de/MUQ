#include "MUQ/Utilities/SaveSchedulerBase.h"

using namespace muq;
using namespace Utilities;

bool SaveSchedulerBase::ShouldSave(int step) {
  if (step % thinIncr == 0)
    return true
  else
    return false
}
