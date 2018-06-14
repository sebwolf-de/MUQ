#ifndef SAVESCHEDULERBASE_H_
#define SAVESCHEDULERBASE_H_

namespace muq {
  namespace SamplingAlgorithms {

    class SaveSchedulerBase {

    public:
      SaveSchedulerBase() {};

      virtual bool ShouldSave(int step) = 0;

    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
