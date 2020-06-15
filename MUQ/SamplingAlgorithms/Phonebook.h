#ifndef PHONEBOOK_H_
#define PHONEBOOK_H_

#if MUQ_HAS_MPI


#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string>
#include <map>
#include <mpi.h>
#if MUQ_HAS_OTF2
#include <otf2/otf2.h>
#include <otf2/OTF2_MPI_Collectives.h>

class OTF2Tracer {


  static OTF2_TimeStamp
  get_time( void )
  {
    double t = MPI_Wtime() * 1e9;
    return ( uint64_t )t;
  }
  static OTF2_FlushType
  pre_flush( void*            userData,
             OTF2_FileType    fileType,
             OTF2_LocationRef location,
             void*            callerData,
             bool             final )
  {
    return OTF2_FLUSH;
  }
  static OTF2_TimeStamp
  post_flush( void*            userData,
              OTF2_FileType    fileType,
              OTF2_LocationRef location )
  {
    return get_time();
  }
  constexpr static OTF2_FlushCallbacks flush_callbacks =
  {
    .otf2_pre_flush  = pre_flush,
    .otf2_post_flush = post_flush
  };

public:

  OTF2Tracer(std::string archive_path, std::string archive_name) {
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    archive = OTF2_Archive_Open( archive_path.c_str(),
                                 archive_name.c_str(),
                                 OTF2_FILEMODE_WRITE,
                                 1024 * 1024 /* event chunk size */,
                                 4 * 1024 * 1024 /* def chunk size */,
                                 OTF2_SUBSTRATE_POSIX,
                                 OTF2_COMPRESSION_NONE );
    OTF2_Archive_SetFlushCallbacks( archive, &flush_callbacks, NULL );
    OTF2_MPI_Archive_SetCollectiveCallbacks( archive,
                                             MPI_COMM_WORLD,
                                             MPI_COMM_NULL );
    OTF2_Archive_OpenEvtFiles( archive );

    evt_writer = OTF2_Archive_GetEvtWriter( archive, rank );

    epoch_start = get_time();


  }



  void enterRegion(int region) {
    ensureRegionName(region);
    OTF2_EvtWriter_Enter( evt_writer,
                          NULL,
                          get_time(),
                          region );
  }

  void leaveRegion(int region) {
    ensureRegionName(region);
    OTF2_EvtWriter_Leave( evt_writer,
                          NULL,
                          get_time(),
                          region );
  }

private:
  void ensureRegionName(int region) {
    if (regionNames.count(region) == 0)
      regionNames[region] = "Unnamed region " + std::to_string(region);
  }

public:

  void setRegionName(int region, std::string name) {
    regionNames[region] = name;
  }

  /*~OTF2Tracer() {
   *
}*/



  void write() {
    uint64_t epoch_end = get_time();
    OTF2_Archive_CloseEvtWriter( archive, evt_writer );
    OTF2_Archive_CloseEvtFiles( archive );
    OTF2_Archive_OpenDefFiles( archive );
    OTF2_DefWriter* def_writer = OTF2_Archive_GetDefWriter( archive,
                                                            rank );
    OTF2_Archive_CloseDefWriter( archive, def_writer );
    OTF2_Archive_CloseDefFiles( archive );
    uint64_t global_epoch_start;
    MPI_Reduce( &epoch_start,
                &global_epoch_start,
                1, OTF2_MPI_UINT64_T, MPI_MIN,
                0, MPI_COMM_WORLD );
    uint64_t global_epoch_end;
    MPI_Reduce( &epoch_end,
                &global_epoch_end,
                1, OTF2_MPI_UINT64_T, MPI_MAX,
                0, MPI_COMM_WORLD );
    if ( 0 == rank )
    {
      OTF2_GlobalDefWriter* global_def_writer = OTF2_Archive_GetGlobalDefWriter( archive );
      OTF2_GlobalDefWriter_WriteClockProperties( global_def_writer,
                                                 1000000000,
                                                 global_epoch_start,
                                                 global_epoch_end - global_epoch_start + 1 );
      OTF2_GlobalDefWriter_WriteString( global_def_writer, 0, "" );
      OTF2_GlobalDefWriter_WriteString( global_def_writer, 1, "Master Thread" );
      OTF2_GlobalDefWriter_WriteString( global_def_writer, 2, "MPI_Barrier" );
      OTF2_GlobalDefWriter_WriteString( global_def_writer, 3, "PMPI_Barrier" );
      OTF2_GlobalDefWriter_WriteString( global_def_writer, 4, "barrier" );
      OTF2_GlobalDefWriter_WriteString( global_def_writer, 5, "MyHost" );
      OTF2_GlobalDefWriter_WriteString( global_def_writer, 6, "node" );
      OTF2_GlobalDefWriter_WriteString( global_def_writer, 7, "MPI" );
      OTF2_GlobalDefWriter_WriteString( global_def_writer, 8, "MPI_COMM_WORLD" );
      int num_strings = 9;

      for ( const auto &regionNamePair : regionNames ) {

        OTF2_GlobalDefWriter_WriteString( global_def_writer, num_strings, regionNamePair.second.c_str() );
        OTF2_GlobalDefWriter_WriteRegion( global_def_writer,
                                          regionNamePair.first /* id */,
                                          num_strings /* region name  */,
                                          num_strings /* alternative name */,
                                          num_strings /* description */,
                                          OTF2_REGION_ROLE_CODE,
                                          OTF2_PARADIGM_MPI,
                                          OTF2_REGION_FLAG_NONE,
                                          7 /* source file */,
                                          0 /* begin lno */,
                                          0 /* end lno */ );
        num_strings++;
      }
      OTF2_GlobalDefWriter_WriteSystemTreeNode( global_def_writer,
                                                0 /* id */,
                                                5 /* name */,
                                                6 /* class */,
                                                OTF2_UNDEFINED_SYSTEM_TREE_NODE /* parent */ );

      for ( int r = 0; r < size; r++ )
      {
        char process_name[ 32 ];
        sprintf( process_name, "MPI Rank %d", r );
        OTF2_GlobalDefWriter_WriteString( global_def_writer,
                                          num_strings + r,
                                          process_name );
        OTF2_GlobalDefWriter_WriteLocationGroup( global_def_writer,
                                                 r /* id */,
                                                 num_strings + r /* name */,
                                                 OTF2_LOCATION_GROUP_TYPE_PROCESS,
                                                 0 /* system tree */ );
        OTF2_GlobalDefWriter_WriteLocation( global_def_writer,
                                            r /* id */,
                                            1 /* name */,
                                            OTF2_LOCATION_TYPE_CPU_THREAD,
                                            4 /* # events */,
                                            r /* location group */ );
      }
      uint64_t comm_locations[ size ];
      for ( int r = 0; r < size; r++ )
      {
        comm_locations[ r ] = r;
      }
      OTF2_GlobalDefWriter_WriteGroup( global_def_writer,
                                       0 /* id */,
                                       7 /* name */,
                                       OTF2_GROUP_TYPE_COMM_LOCATIONS,
                                       OTF2_PARADIGM_MPI,
                                       OTF2_GROUP_FLAG_NONE,
                                       size,
                                       comm_locations );
      OTF2_GlobalDefWriter_WriteGroup( global_def_writer,
                                       1 /* id */,
                                       0 /* name */,
                                       OTF2_GROUP_TYPE_COMM_GROUP,
                                       OTF2_PARADIGM_MPI,
                                       OTF2_GROUP_FLAG_NONE,
                                       size,
                                       comm_locations );
      OTF2_GlobalDefWriter_WriteComm( global_def_writer,
                                      0 /* id */,
                                      8 /* name */,
                                      1 /* group */,
                                      OTF2_UNDEFINED_COMM /* parent */ );
      OTF2_Archive_CloseGlobalDefWriter( archive,
                                         global_def_writer );
    }

    OTF2_Archive_Close( archive );
  }

private:
  int size;
  int rank;
  uint64_t epoch_start;

  OTF2_Archive* archive;
  OTF2_EvtWriter* evt_writer;

  std::map<int, std::string> regionNames;
};
#else
class OTF2Tracer {
public:

  OTF2Tracer(std::string archive_path, std::string archive_name) {}

  void enterRegion(int region) {}

  void leaveRegion(int region) {}

  void setRegionName(int region, std::string name) {}
  void write() {}
};
#endif





enum TracerRegions {
  Setup,
  Finalize,
  PhonebookBusy,
  BurnIn,
  Sampling,
  CollectorBusy,
  RetrievingProposal
};


static OTF2Tracer* tracer;







#include "spdlog/fmt/ostr.h"
#include "spdlog/spdlog.h"
#include <parcer/Communicator.h>
#include "MUQ/SamplingAlgorithms/ParallelFlags.h"
#include "MUQ/Utilities/Cereal/MultiIndexSerializer.h"
#include <deque>
#include <chrono>

namespace muq {
  namespace SamplingAlgorithms {

    class PhonebookServer {
    public:

      PhonebookServer(std::shared_ptr<parcer::Communicator> comm)
        : comm(comm)
      {
      }

      void Run() {
        //Dune::Timer timer_idle;
        //Dune::Timer timer_full;


        while (true) {
          MPI_Status status;
          //timer_idle.start();


          if (scheduler_active) {
            for ( auto &indexWorkerListPair : phonebook ) {
              std::shared_ptr<MultiIndex> index = indexWorkerListPair.first;
              WorkerList& workerList = indexWorkerListPair.second;
              bool work_run_out = workerList.NumWorkersReady() == workerList.NumWorkers() && getNumQueuedTasks(index, workerList) == 0.0;

              if (workerList.NormalizedRegisteredReadyCounter() > .5 || work_run_out || workerList.recheck) {

                if (workerList.NormalizedRegisteredReadyCounter() > .5) {
                  workerList.ResetTimer();
                }

                if (workerList.NumWorkers() <= 1) // Only reschedule if we have more than 1 worker on this model
                  continue;

                if (workerList.NumWorkersReady() == 0) {
                  workerList.recheck = true;
                  continue;
                }
                workerList.recheck = false;

                spdlog::debug("Timer triggered for {}, idle fraction {}, {} workers on that model, work run out {}", *indexWorkerListPair.first, workerList.GetIdleFraction(), workerList.NumWorkers(), work_run_out);

                //if (index->GetValue(0) == 1) // FIXME: Fix rescheduling model 1
                //  continue;

                double my_load = getLoadFactor(index, workerList);
                double largest_others_load = .0;
                std::shared_ptr<MultiIndex> most_loaded_index = nullptr;
                for ( auto &indexWorkerListPair : phonebook ) {
                  std::shared_ptr<MultiIndex> index = indexWorkerListPair.first;
                  WorkerList& workerList = indexWorkerListPair.second;
                  double others_load = getLoadFactor(index, workerList);
                  spdlog::debug("Load on model {}: {}", *index, others_load);
                  if (others_load > largest_others_load) {
                    largest_others_load = others_load;
                    most_loaded_index = index;
                  }
                }
                assert (most_loaded_index != nullptr);

                if (*index == *most_loaded_index)
                  continue;

                if (work_run_out || my_load + .5 / (double)workerList.NumWorkers() + .01 < largest_others_load) {

                  spdlog::debug("Reassigning from model {} to model {}, {} ready here", *index, *most_loaded_index, workerList.NumWorkersReady());

                  const int RootNode = 0; // Send to root

                  int rescheduleRank = workerList.GetWorkersReady()[0];
                  UnRegister(index, rescheduleRank); // Already unregister this rank so it won't be offered for sampling anymore while being rescheduled!

                  comm->Send(ControlFlag::SCHEDULING_NEEDED, RootNode, ControlTag);
                  comm->Send(*index, RootNode, ControlTag);
                  comm->Send(rescheduleRank, RootNode, ControlTag);
                  comm->Send(*most_loaded_index, RootNode, ControlTag);

                  //workerList.ResetTimer();

                  scheduler_active = false;
                  break; // Phonebook needs to be ready for further communication after rescheduling a process, so don't reschedule another one right now
                }


                /*if (workerList.GetIdleFraction() > .3 || work_run_out) {
                  spdlog::debug("HIGH IDLE FRACTION IN {}", *indexWorkerListPair.first);

                  // Find busiest model
                  std::shared_ptr<MultiIndex> busiestModelIndex = nullptr;
                  double busiest_idle_fraction = 1.0;
                  for ( auto &indexWorkerListPair : phonebook ) {
                    WorkerList& workerList = indexWorkerListPair.second;

                    spdlog::debug("Model {} has idle fraction {}", *indexWorkerListPair.first, workerList.GetIdleFraction());

                    // TODO: Add #pending high priority requests / #processes to model load
                    int high_priority_in_queue = 0;
                    for (auto request_iter = requests.begin(); request_iter != requests.end(); request_iter++) {
                      std::shared_ptr<MultiIndex> index = std::get<0>(*request_iter);
                      bool high_priority = std::get<2>(*request_iter);
                      if (*index == *indexWorkerListPair.first && high_priority)
                        high_priority_in_queue++;
                    }
                    double current_idle_fraction = workerList.GetIdleFraction();
                    current_idle_fraction -= (double)high_priority_in_queue / workerList.NumWorkers();


                    if (current_idle_fraction < busiest_idle_fraction) {
                      busiest_idle_fraction = current_idle_fraction;
                      busiestModelIndex = indexWorkerListPair.first;
                    }
                  }
                  assert (busiestModelIndex != nullptr);


                  const int RootNode = 0; // Send to root

                  comm->Send(ControlFlag::SCHEDULING_NEEDED, RootNode, ControlTag);
                  comm->Send(*indexWorkerListPair.first, RootNode, ControlTag);
                  comm->Send(*busiestModelIndex, RootNode, ControlTag);

                  scheduler_active = false;
                  break; // Phonebook needs to be ready for further communication after rescheduling a process, so don't reschedule another one right now
                }*/



              }

            }
          }

          ControlFlag command = comm->Recv<ControlFlag>(MPI_ANY_SOURCE, ControlTag, &status);
          //timer_idle.stop();

          for ( auto &indexWorkerListPair : phonebook ) {
            WorkerList& workerList = indexWorkerListPair.second;
            workerList.tick();
          }



          tracer->enterRegion(TracerRegions::PhonebookBusy);

          if (command == ControlFlag::GET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));
            bool high_priority = comm->Recv<bool>(status.MPI_SOURCE, ControlTag);
            if (high_priority)
              requests.push_front(std::tie(index, status.MPI_SOURCE, high_priority));
            else
              requests.push_back(std::tie(index, status.MPI_SOURCE, high_priority));
          } else if (command == ControlFlag::SCHEDULING_DONE) {
            scheduler_active = true;
          } else if (command == ControlFlag::SCHEDULING_STOP) {
            scheduler_active = false;
          } else if (command == ControlFlag::SET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);

            if (!phonebook.count(index)) {
              phonebook[index] = WorkerList();
            }
            phonebook[index].AddWorker(rank);
            //comm->Send(ControlFlag::HANDSHAKE, status.MPI_SOURCE, ControlTag);
            spdlog::trace("Phonebook entry for {} set", *index);
          } else if (command == ControlFlag::UNSET_WORKGROUP) {
            auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);

            UnRegister(index, rank);
          } else if (command == ControlFlag::GET_WORKGROUPS) {
            auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));

            if (!phonebook.count(index)) {
              std::cerr << "getting workers for nonexistent model!" << std::endl;
            }
            spdlog::trace("Getting workers from phonebook map");
            std::vector<int> sendvec = phonebook[index].GetWorkers();
            spdlog::trace("Sending {} workgroups", sendvec.size());
            comm->Send(sendvec, status.MPI_SOURCE, ControlTag);
          } else if (command == ControlFlag::GET_LARGEST_INDEX) {
            if (phonebook.empty()) {
              comm->Send(-1, status.MPI_SOURCE, ControlTag);
              spdlog::trace("Sent empty largest index");
            } else {
              comm->Send(*phonebook.rbegin()->first, status.MPI_SOURCE, ControlTag);
              spdlog::trace("Sent largest index {} unset", *phonebook.rbegin()->first);
            }
          } else if (command == ControlFlag::SET_WORKER_READY) {

            auto index = std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(status.MPI_SOURCE, ControlTag));
            int rank = comm->Recv<int>(status.MPI_SOURCE, ControlTag, &status);
            if (!phonebook.count(index)) {
              std::cerr << "setting ready for nonexistent model!" << std::endl;
              continue;
            }
            phonebook[index].SetWorkerReady(rank);
          } else if (command == ControlFlag::QUIT) {
            tracer->leaveRegion(TracerRegions::PhonebookBusy);
            spdlog::trace("Rank {} quit", comm->GetRank());
            break;
          }

          for (auto request_iter = requests.begin(); request_iter != requests.end();) {
            std::shared_ptr<MultiIndex> index = std::get<0>(*request_iter);
            int sender = std::get<1>(*request_iter);
            //bool high_priority = std::get<2>(*request_iter);
            /*if (!high_priority) {
              request_iter++;
              continue;
            }*/

            if (!phonebook.count(index)) {
              std::cerr << "checking request for nonexistent model!" << std::endl;
            }
            int workerRank = phonebook[index].NextWorker();
            if (workerRank == -1) {
              request_iter++;
            } else {
              comm->Send(workerRank, sender, ControlTag);
              request_iter = requests.erase(request_iter);
            }
          }



          /*for (auto request_iter = requests.begin(); request_iter != requests.end();) {
            std::shared_ptr<MultiIndex> index = std::get<0>(*request_iter);
            int sender = std::get<1>(*request_iter);
            bool high_priority = std::get<2>(*request_iter);
            if (high_priority) {
              request_iter++;
              continue;
            }

            if (!phonebook.count(index)) {
              std::cerr << "checking request for nonexistent model!" << std::endl;
            }
            int workerRank = phonebook[index].NextWorker();
            if (workerRank == -1) {
              request_iter++;
            } else {
              comm->Send(workerRank, sender, ControlTag);
              request_iter = requests.erase(request_iter);
            }
          }*/

          tracer->leaveRegion(TracerRegions::PhonebookBusy);
        }
        //std::cout << "Phonebook " << comm->GetRank() << " idle time:\t" << timer_idle.elapsed() << " of:\t" << timer_full.elapsed() << std::endl;

      }

    private:

      void UnRegister(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        /*if (!phonebook.count(modelIndex)) {
          std::cerr << "unsetting nonexistent entry!" << std::endl;
        }*/
        phonebook[modelIndex].RemoveWorker(rank);
        if (phonebook[modelIndex].NumWorkers() == 0) {
          spdlog::debug("Phonebook erasing modelIndex", *modelIndex);
          phonebook.erase(modelIndex);
        }
        //comm->Send(ControlFlag::HANDSHAKE, status.MPI_SOURCE, ControlTag);
        spdlog::trace("Phonebook entry for {} unset", *modelIndex);
      }

      class WorkerList {
      public:

        WorkerList() {
          ResetTimer();
        }

        int NextWorker() {
          //assert(workers.size() > 0);
          if (workers.size() == 0)
            return -1;

          if (workers_ready.size() == 0)
            return -1;

          //tick();
          int worker = workers_ready.front();
          workers_ready.pop_front();
          return worker;
        }
        void AddWorker(int worker) {
          //tick();
          workers.push_back(worker);
        }
        void RemoveWorker(int worker) {
          //tick();
          workers_ready.erase(std::remove(workers_ready.begin(), workers_ready.end(), worker), workers_ready.end());
          workers.erase(std::remove(workers.begin(), workers.end(), worker), workers.end());
          spdlog::trace("Removed worker {}, {} remaining", worker, NumWorkers());
        }
        int NumWorkers() const {
          return workers.size();
        }
        int NumWorkersReady() const {
          return workers_ready.size();
        }
        void SetWorkerReady(int worker) {
          //tick();
          registeredReadyCounter++;

          workers_ready.push_back(worker);
        }
        std::vector<int> GetWorkers() {
          return workers;
        }
        std::deque<int> GetWorkersReady() {
          return workers_ready;
        }

        double GetIdleFraction() {
          if (total_time == std::chrono::nanoseconds::zero())
            return .0;
          return ((double)idle_time.count()) / (double)(total_time.count());
        }
        double NormalizedRegisteredReadyCounter() {
          if (workers.size() == 0)
            return -1;
          return (double)registeredReadyCounter / (double)workers.size();
        }
        void ResetTimer() {
          registeredReadyCounter = 0;
          idle_time = std::chrono::nanoseconds::zero();
          total_time = std::chrono::nanoseconds::zero();
          begin_tick = std::chrono::high_resolution_clock::now();
          last_tick = begin_tick;
        }
        void tick() {
          std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();

          std::chrono::nanoseconds measurement_period = now - last_tick;
          idle_time += measurement_period * workers_ready.size();
          total_time += measurement_period * workers.size();

          last_tick = now;
        }
      public:
        bool recheck = false;

      private:
        int registeredReadyCounter;
        std::chrono::nanoseconds idle_time;
        std::chrono::nanoseconds total_time;
        std::chrono::high_resolution_clock::time_point begin_tick, last_tick;

        std::vector<int> workers;
        std::deque<int> workers_ready;
      };


      double getNumQueuedTasks (std::shared_ptr<MultiIndex> index, WorkerList& worker_list) {
        double in_queue = 0;
        for (auto request_iter = requests.begin(); request_iter != requests.end(); request_iter++) {
          std::shared_ptr<MultiIndex> index_iter = std::get<0>(*request_iter);
          bool high_priority = std::get<2>(*request_iter);
          if (*index_iter == *index)
            in_queue += high_priority ? 1.0 : .1;
        }
        return in_queue;
      }

      double getLoadFactor (std::shared_ptr<MultiIndex> index, WorkerList& worker_list) {
        double load_factor = 1.0 - worker_list.GetIdleFraction();
        load_factor += getNumQueuedTasks(index, worker_list) / worker_list.NumWorkers();
        return load_factor;
      }



      std::deque< std::tuple<std::shared_ptr<MultiIndex>, int, bool> > requests;


      std::map<std::shared_ptr<MultiIndex>, WorkerList, MultiPtrComp> phonebook;
      std::shared_ptr<parcer::Communicator> comm;
      bool scheduler_active = true;
    };


    class PhonebookClient {

    public:

      PhonebookClient(std::shared_ptr<parcer::Communicator> comm, int phonebookRank)
        : comm(comm), phonebookRank(phonebookRank)
      {
      }

      int Query(std::shared_ptr<MultiIndex> remoteIndex, bool high_priority) {
        comm->Send(ControlFlag::GET_WORKGROUP, phonebookRank, ControlTag);
        comm->Send<MultiIndex>(*remoteIndex, phonebookRank, ControlTag);
        comm->Send(high_priority, phonebookRank, ControlTag);
        return comm->Recv<int>(phonebookRank, ControlTag);
      }

      std::vector<int> GetWorkgroups(std::shared_ptr<MultiIndex> modelIndex) {
        spdlog::debug("GetWorkgroups call for model {}", *modelIndex);
        comm->Send(ControlFlag::GET_WORKGROUPS, phonebookRank, ControlTag);
        comm->Send(*modelIndex, phonebookRank, ControlTag);
        spdlog::debug("GetWorkgroups call for model {}, retrieving", *modelIndex);
        std::vector<int> ret = comm->Recv<std::vector<int>>(phonebookRank, ControlTag);
        spdlog::debug("GetWorkgroups call for model {}, returning", *modelIndex);
        return ret;
      }

      std::shared_ptr<MultiIndex> LargestIndex() {
        comm->Send(ControlFlag::GET_LARGEST_INDEX, phonebookRank, ControlTag);
        return std::make_shared<MultiIndex>(comm->Recv<MultiIndex>(phonebookRank, ControlTag));
      }

      void Register(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        comm->Send(ControlFlag::SET_WORKGROUP, phonebookRank, ControlTag);
        comm->Send(*modelIndex, phonebookRank, ControlTag);
        comm->Ssend(rank, phonebookRank, ControlTag);
      }

      void UnRegister(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        comm->Send(ControlFlag::UNSET_WORKGROUP, phonebookRank, ControlTag);
        comm->Send(*modelIndex, phonebookRank, ControlTag);
        comm->Ssend(rank, phonebookRank, ControlTag);
        //if (comm->Recv<ControlFlag>(phonebookRank, ControlTag) != ControlFlag::HANDSHAKE)
        //  std::cerr << "Failed handshake in UnRegister()!" << std::endl;
      }

      void SetWorkerReady(std::shared_ptr<MultiIndex> modelIndex, int rank) {
        comm->Send(ControlFlag::SET_WORKER_READY, phonebookRank, ControlTag);
        comm->Send(*modelIndex, phonebookRank, ControlTag);
        comm->Send(rank, phonebookRank, ControlTag);
      }

      void SchedulingDone() {
        comm->Send(ControlFlag::SCHEDULING_DONE, phonebookRank, ControlTag);
      }
      void SchedulingStop() {
        comm->Ssend(ControlFlag::SCHEDULING_STOP, phonebookRank, ControlTag);
      }

    private:
      std::shared_ptr<parcer::Communicator> comm;
      int phonebookRank;
    };
  }
}

#endif

#endif
