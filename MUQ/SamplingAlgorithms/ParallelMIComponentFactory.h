#ifndef PARALLELMICOMPONENTFACTORY_H_
#define PARALLELMICOMPONENTFACTORY_H_

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif

#include "spdlog/spdlog.h"
#include "MUQ/SamplingAlgorithms/ParallelAbstractSamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {

    class ParallelMIComponentFactory : public MIComponentFactory {

    public:

      ParallelMIComponentFactory (std::shared_ptr<parcer::Communicator> comm, std::shared_ptr<parcer::Communicator> global_comm, std::shared_ptr<MIComponentFactory> componentFactory)
       : comm(comm), global_comm(global_comm), componentFactory(componentFactory)
      {

        if (comm->GetRank() != 0) {
          //Dune::Timer timer_idle;
          //Dune::Timer timer_full;
          while (true) {
            spdlog::trace("Parallel factory rank {} waiting...", comm->GetRank());
            //timer_idle.start();
            ControlFlag command = comm->Recv<ControlFlag>(0, WorkgroupTag);
            //timer_idle.stop();
            spdlog::trace("Parallel factory rank {} received command", comm->GetRank(), command);
            if (command == ControlFlag::FINALIZE) {
              samplingProblems.clear(); // Tear down models synchronously
              comm->Barrier();
              spdlog::trace("Parallel factory rank {} passed finalize barrier", comm->GetRank());
              break;
            }
            if (command == ControlFlag::INIT_PROBLEM) {
              auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(0, WorkgroupTag));
              int id = comm->Recv<int>(0, WorkgroupTag);
              spdlog::trace("Parallel factory rank {} building model index {}", comm->GetRank(), *index);
              samplingProblems[id] = componentFactory->SamplingProblem(index);//std::make_shared<MySamplingProblem>(index, comm, id, measurements);
            }
            else if (command == ControlFlag::LOGDENSITY) {
              int id = comm->Recv<int>(0, WorkgroupTag);
              auto state = std::make_shared<SamplingState>(comm->Recv<Eigen::VectorXd>(0, WorkgroupTag));
              samplingProblems[id]->LogDensity(0.0, state, AbstractSamplingProblem::SampleType::Proposed); // TODO: Correct arguments!
            }
            else if (command == ControlFlag::TEST) {
              int id = comm->Recv<int>(0, WorkgroupTag);
              auto state = std::make_shared<SamplingState>(comm->Recv<Eigen::VectorXd>(0, WorkgroupTag));

              double density;
              //std::dynamic_pointer_cast<MySamplingProblem>(samplingProblems[id])->test(state, "test", density);
              std::cerr << "Not implemented!!!" << std::endl;
            }
            else if (command == ControlFlag::QOI) {
              int id = comm->Recv<int>(0, WorkgroupTag);

              samplingProblems[id]->QOI();
            } else {
              std::cerr << "Unexpected command!" << std::endl;
              exit(43);
            }
          }
          //std::cout << "Worker " << global_comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;
        }
      }

      void finalize() {
        if (comm->GetRank() != 0)
          return;
        spdlog::trace("Parallel factory sending finalize");
        for (int dest = 1; dest < comm->GetSize(); dest++)
          comm->Send(ControlFlag::FINALIZE, dest, WorkgroupTag);
        samplingProblems.clear(); // Tear down models synchronously
        comm->Barrier();
        spdlog::trace("Parallel factory finalized");
      }

      virtual std::shared_ptr<MCMCProposal> Proposal (std::shared_ptr<MultiIndex> const& index, std::shared_ptr<AbstractSamplingProblem> const& samplingProblem) override {
        return componentFactory->Proposal(index, samplingProblem);
      }

      virtual std::shared_ptr<MultiIndex> FinestIndex() override {
        return componentFactory->FinestIndex();
      }

      virtual std::shared_ptr<MCMCProposal> CoarseProposal (std::shared_ptr<MultiIndex> const& index,
                                                            std::shared_ptr<AbstractSamplingProblem> const& coarseProblem,
                                                            std::shared_ptr<SingleChainMCMC> const& coarseChain) override {
        return componentFactory->CoarseProposal(index, coarseProblem, coarseChain);
      }

      virtual std::shared_ptr<AbstractSamplingProblem> SamplingProblem (std::shared_ptr<MultiIndex> const& index) override {
        //int idcnt = index->GetValue(0);
        idcnt++;
        if (comm->GetRank() == 0) {
          spdlog::debug("Rank {} requesting model {} from parallel factory", comm->GetRank(), *index);
          for (int dest = 1; dest < comm->GetSize(); dest++) {
            comm->Send(ControlFlag::INIT_PROBLEM, dest, WorkgroupTag);
            comm->Send(*index, dest, WorkgroupTag);
            comm->Send(idcnt, dest, WorkgroupTag);
          }
        }
        samplingProblems[idcnt] = std::make_shared<ParallelAbstractSamplingProblem>(comm, idcnt, componentFactory->SamplingProblem(index));
        return samplingProblems[idcnt];
      }

      virtual std::shared_ptr<MIInterpolation> Interpolation (std::shared_ptr<MultiIndex> const& index) override {
        return componentFactory->Interpolation(index);
      }

      virtual Eigen::VectorXd StartingPoint (std::shared_ptr<MultiIndex> const& index) override {
        return componentFactory->StartingPoint(index);
      }

    private:
      int idcnt = 0;
      std::shared_ptr<parcer::Communicator> comm;
      std::shared_ptr<parcer::Communicator> global_comm;
      std::shared_ptr<MIComponentFactory> componentFactory;

      std::map<int, std::shared_ptr<AbstractSamplingProblem>> samplingProblems;

    };

  }
}

#endif

#endif
