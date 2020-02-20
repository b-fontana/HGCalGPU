#include "Prod_test.h"

Prodtest::Prodtest(const edm::ParameterSet& ps):
  token_(consumes<HGCUncalibratedRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEUncalibRecHitsTok")))
{
  nhits_ = 5000000;
  tools_.reset(new hgcal::RecHitTools());

  h_ = new float();
  d1_ = new float();
  d2_ = new float();

  //_allocate memory for hits on the host
  host(nhits_, h_, h_mem_);
  //_allocate memory for hits on the device
  device(nhits_, d1_, d2_, d_mem_);

  produces<HGCUncalibratedRecHitCollection>(collection_name_);
}

Prodtest::~Prodtest()
{
}

void Prodtest::host(const int& nhits, float *&h, cudautils::host::unique_ptr<float[]>& mem)
{
  mem = cudautils::make_host_unique<float[]>(nhits_ * sizeof(float), 0);
  h = mem.get();
}

void Prodtest::device(const int& nhits, float *&d1, float *&d2, cudautils::device::unique_ptr<float[]>& mem)
{
  mem = cudautils::make_device_unique<float[]>(nhits_ * sizeof(float) * 2, 0);
  d1 = mem.get();
  d2 = d1 + nhits_;
}

void Prodtest::acquire(edm::Event const& event, edm::EventSetup const& setup, edm::WaitingTaskWithArenaHolder w) 
{
  tools_->getEventSetup(setup);
  set_geometry_(setup);
  event.getByToken(token_, handle_ee_);
  //const auto &hits_ee = *handle_ee_;
  Prod_run prod_run(nhits_, h_, d1_, d2_);
  prod_run.run_kernels();
  std::cout << "end acquire" << std::endl;
}

void Prodtest::produce(edm::Event& event, const edm::EventSetup& setup)
{
  std::cout << "start produce" << std::endl;
  auto rechits = std::make_unique<HGCUncalibratedRecHitCollection>( *handle_ee_ );
  event.put(std::move(rechits), collection_name_);
  std::cout << "end produce" << std::endl;
}

void Prodtest::set_geometry_(const edm::EventSetup& setup)
{
  tools_->getEventSetup(setup);
  //rechitMaker_->set(es);
  std::string handle_str;
  handle_str = "HGCalEESensitive";
  edm::ESHandle<HGCalGeometry> handle;
  setup.get<IdealGeometryRecord>().get(handle_str, handle);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Prodtest);
