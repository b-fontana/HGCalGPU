#ifndef _TYPES_CUH_
#define _TYPES_CUH_

#include <vector>

class HGCConstantVectorData {
 public:
  std::vector<double> fCPerMIP;
  std::vector<double> cce;
  std::vector<double> noise_fC;
  std::vector<double> rcorr;
  std::vector<double> weights;
};

class HGCeeUncalibratedRecHitConstantData {
 public:
  double hgcEE_keV2DIGI_;
  double hgceeUncalib2GeV_;
  double *hgcEE_fCPerMIP_;  
  double *hgcEE_cce_;
  double *hgcEE_noise_fC_;
  double *rcorr_;
  float *weights_;
  uint32_t rangeMatch_;
  uint32_t rangeMask_;
  bool hgcEE_isSiFE_;
  size_t ndoubles;
  size_t nbytes;
  size_t ndelem;
  size_t nfelem;
  size_t nuelem;
  size_t nbelem;
  size_t s_hgcEE_fCPerMIP_;
  size_t s_hgcEE_cce_;
  size_t s_hgcEE_noise_fC_;
  size_t s_rcorr_;
  size_t s_weights_;
};

class HGChefUncalibratedRecHitConstantData {
 public:
  double hgcHEF_keV2DIGI_;
  double hgchefUncalib2GeV_;
  double *hgcHEF_fCPerMIP_;
  double *hgcHEF_cce_;
  double *hgcHEF_noise_fC_;
  double *rcorr_;
  float *weights_;
  uint32_t rangeMatch_;
  uint32_t rangeMask_;
  bool hgcHEF_isSiFE_;
  size_t nbytes;
  size_t ndelem;
  size_t nfelem;
  size_t nuelem;
  size_t nbelem;
  size_t s_hgcHEF_fCPerMIP_;
  size_t s_hgcHEF_cce_;
  size_t s_hgcHEF_noise_fC_;
  size_t s_rcorr_;
  size_t s_weights_;
};

class HGChebUncalibratedRecHitConstantData {
 public:
  double hgcHEB_keV2DIGI_;
  double hgchebUncalib2GeV_;
  double hgcHEB_noise_MIP_;
  double *hgcHEB_fCPerMIP_;
  double *hgcHEB_cce_;
  double *rcorr_;
  float *weights_;
  uint32_t rangeMatch_;
  uint32_t rangeMask_;
  bool hgcHEB_isSiFE_;
  size_t nbytes;
  size_t ndelem;
  size_t nfelem;
  size_t nuelem;
  size_t nbelem;
  size_t s_hgcHEF_fCPerMIP_;
  size_t s_hgcHEF_cce_;
  size_t s_hgcHEF_noise_fC_;
  size_t s_rcorr_;
  size_t s_rweights_;
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
