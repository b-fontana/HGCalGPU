///
/// \class l1t::Stage2Layer1FirmwareFactory
///
///
/// \author: Jim Brooke
///

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TCalorimeter/interface/Stage2PreProcessorFirmware.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer1FirmwareFactory.h"

#include "CondFormats/L1TObjects/interface/CaloParams.h"

using namespace std;
using namespace edm;

l1t::Stage2Layer1FirmwareFactory::ReturnType
l1t::Stage2Layer1FirmwareFactory::create(const FirmwareVersion & fwv, CaloParams* params) {

  ReturnType p;
  unsigned v = fwv.firmwareVersion();
  
  switch (v){
  case 1:
    p = ReturnType(new Stage2PreProcessorFirmwareImp1(fwv, params));
      break;
  default:
    // Invalid Firmware, log an error:
    LogError("l1t|caloStage2") << "Invalid firmware version requested: " << v << "\n";
    break;
  }
  
  return p;

}
