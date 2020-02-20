#ifndef Prodtest_h
#define Prodtest_h

#include <iostream>
#include <string>
#include <memory>
#include <cuda_runtime.h>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAContextState.h"
#include "HeterogeneousCore/CUDAServices/interface/CUDAService.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_noncached_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_unique_ptr.h"

#include "Prod_run.h"

class Prodtest: public edm::stream::EDProducer<edm::ExternalWork> 
{
 public:
  explicit Prodtest(const edm::ParameterSet& ps);
  ~Prodtest() override;

  virtual void acquire(edm::Event const&, edm::EventSetup const&, edm::WaitingTaskWithArenaHolder) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;

 private:
  int nhits_ = 0;
  edm::EDGetTokenT<HGCeeUncalibratedRecHitCollection> token_;
  const std::string collection_name_ = "HeterogeneousHGCeeUncalibratedRecHits";
  edm::Handle<HGCeeUncalibratedRecHitCollection> handle_ee_;

  //geometry
  void set_geometry_(const edm::EventSetup&);
  std::unique_ptr<hgcal::RecHitTools> tools_;

  //memory
  void host(const int&, float*&, cudautils::host::unique_ptr<float[]>&);
  void device(const int&, float*&, float*&, cudautils::device::unique_ptr<float[]>&);
  cudautils::host::unique_ptr<float[]> h_mem_;
  cudautils::device::unique_ptr<float[]> d_mem_;

  float *h_ = nullptr;
  float *d1_ = nullptr;
  float *d2_ = nullptr;
};

#endif //Prodtest_h
