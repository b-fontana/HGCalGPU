#ifndef _TYPES_CUH_
#define _TYPES_CUH_

#include <vector>

class HGCConstantVectorData {
 public:
  std::vector<double> fCPerMIP;
  std::vector<double> cee;
  std::vector<double> noise_fC;
  std::vector<float> rcorr;
  std::vector<float> weights;
  std::vector<size_t> sizes = {0, 0, 0, 0, 0};
};

class HGCeeUncalibratedRecHitConstantData {
 public:
  double hgcEE_keV2DIGI_;
  double hgceeUncalib2GeV_;
  double *hgcEE_fCPerMIP_;
  double *hgcEE_cee_;
  double *hgcEE_noise_fC_;
  float *rcorr_;
  float *weights_;
  uint32_t rangeMatch_;
  uint32_t rangeMask_;
  bool hgcEE_isSiFE_;
  size_t nbytes;
};

class HGChefUncalibratedRecHitConstantData {
 public:
  double hgcHEF_keV2DIGI_;
  double hgchefUncalib2GeV_;
  double *hgcHEF_fCPerMIP_;
  double *hgcHEF_cee_;
  double *hgcHEF_noise_fC_;
  float *rcorr_;
  float *weights_;
  uint32_t rangeMatch_;
  uint32_t rangeMask_;
  bool hgcHEF_isSiFE_;
  size_t nbytes;
};

class HGChebUncalibratedRecHitConstantData {
 public:
  double hgcHEB_keV2DIGI_;
  double hgchebUncalib2GeV_;
  double hgcHEB_noise_MIP_;
  double *hgcHEB_fCPerMIP_;
  double *hgcHEB_cee_;
  float *rcorr_;
  float *weights_;
  uint32_t rangeMatch_;
  uint32_t rangeMask_;
  bool hgcHEB_isSiFE_;
  size_t nbytes;
};

class HGCUncalibratedRecHitSoA {
public:
  float *amplitude;
  float *pedestal;
  float *jitter;
  float *chi2;
  float *OOTamplitude;
  float *OOTchi2;
  uint32_t *flags;
  uint32_t *aux;
  uint32_t *id;
  size_t nbytes;
};

class HGCRecHitSoA {
 public:
  float *energy;
  float *time;
  uint32_t *id;
  uint32_t *flags;
  uint32_t *flagBits;
  size_t nbytes;
};

#endif /* _TYPES_H_ */
