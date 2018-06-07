#ifndef SAVESCHEDULERBASE_H_
#define SAVESCHEDULERBASE_H_

namespace muq {
  namespace Utilities {

    class SaveSchedulerBase {

    public:
      SaveSchedulerBase() {};

      bool ShouldSave(int step);

    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
