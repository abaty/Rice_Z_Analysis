#ifndef VertexCompositeNtuple_h
#define VertexCompositeNtuple_h

// Header file for ROOT classes
#include <TROOT.h>
#include <TChain.h>
#include <TInterpreter.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TFile.h>
#include <TVector3.h>

// Header file for c++ classes
#include <iostream>
#include <string>
#include <vector>
#include <map>

// Header file for the classes stored in the TChain


typedef std::vector< std::vector<UChar_t> > UCharVecVec;
const UInt_t NCAND = 500;
const UInt_t NGEN  = 10;


class VertexCompositeNtuple {

public :

  VertexCompositeNtuple();
  virtual ~VertexCompositeNtuple();
  virtual Bool_t       GetTree         (const std::vector< std::string >&, const std::string& treeName="dimucontana");
  virtual Bool_t       GetTree         (const std::string&, const std::string& treeName="dimucontana");
  virtual Int_t        GetEntry        (Long64_t);
  virtual Long64_t     GetEntries      (void) const { return (fChain_ ? fChain_->GetEntries() : -1); }
  virtual Long64_t     GetTreeEntries  (void) const { return ((fChain_ && fChain_->GetTree()) ? fChain_->GetTree()->GetEntriesFast() : -1); }
  virtual Int_t        GetTreeNumber   (void) const { return fCurrent_; }
  virtual void         Clear           (void);

  // EVENT INFO GETTERS
  UInt_t    RunNb()                       { SetBranch("RunNb");                       return RunNb_;                      }
  UInt_t    LSNb()                        { SetBranch("LSNb");                        return LSNb_;                       }
  UInt_t    EventNb()                     { SetBranch("EventNb");                     return EventNb_;                    }
  Short_t   nPV()                         { SetBranch("nPV");                         return nPV_;                        }
  Float_t   bestvtxX()                    { SetBranch("bestvtxX");                    return bestvtxX_;                   }
  Float_t   bestvtxY()                    { SetBranch("bestvtxY");                    return bestvtxY_;                   }
  Float_t   bestvtxZ()                    { SetBranch("bestvtxZ");                    return bestvtxZ_;                   }
  Short_t   centrality()                  { SetBranch("centrality");                  return centrality_;                 }
  Int_t     Npixel()                      { SetBranch("Npixel");                      return Npixel_;                     }
  Float_t   HFsumETPlus()                 { SetBranch("HFsumETPlus");                 return HFsumETPlus_;                }
  Float_t   HFsumETMinus()                { SetBranch("HFsumETMinus");                return HFsumETMinus_;               }
  Float_t   ZDCPlus()                     { SetBranch("ZDCPlus");                     return ZDCPlus_;                    }
  Float_t   ZDCMinus()                    { SetBranch("ZDCMinus");                    return ZDCMinus_;                   }
  Int_t     Ntrkoffline()                 { SetBranch("Ntrkoffline");                 return Ntrkoffline_;                }
  Short_t*  trigPrescale()                { SetBranch("trigPrescale");                return trigPrescale_;               }
  Bool_t*   trigHLT()                     { SetBranch("trigHLT");                     return trigHLT_;                    }
  Bool_t*   evtSel()                      { SetBranch("evtSel");                      return evtSel_;                     }

  // EVENT PLANE GETTERS
  Float_t   ephfpSumW()                   { SetBranch("ephfpSumW");                   return ephfpSumW_;                  }
  Float_t*  ephfpAngle()                  { SetBranch("ephfpAngle");                  return ephfpAngle_;                 }
  Float_t*  ephfpQ()                      { SetBranch("ephfpQ");                      return ephfpQ_;                     }
  Float_t   ephfmSumW()                   { SetBranch("ephfmSumW");                   return ephfmSumW_;                  }
  Float_t*  ephfmAngle()                  { SetBranch("ephfmAngle");                  return ephfmAngle_;                 }
  Float_t*  ephfmQ()                      { SetBranch("ephfmQ");                      return ephfmQ_;                     }
  Float_t   eptrackmidSumW()                   { SetBranch("eptrackmidSumW");                   return eptrackmidSumW_;                  }
  Float_t*  eptrackmidAngle()                  { SetBranch("eptrackmidAngle");                  return eptrackmidAngle_;                 }
  Float_t*  eptrackmidQ()                      { SetBranch("eptrackmidQ");                      return eptrackmidQ_;                     }

  // CANDIDATE INFO GETTERS
  UInt_t    candSize()                    { SetBranch("candSize");                    return candSize_;                   }
  Float_t*  pT()                          { SetBranch("pT");                          return pT_;                         }
  Float_t*  eta()                         { SetBranch("eta");                         return eta_;                        }
  Float_t*  y()                           { SetBranch("y");                           return y_;                          }
  Float_t*  phi()                         { SetBranch("phi");                         return phi_;                        }
  Float_t*  mass()                        { SetBranch("mass");                        return mass_;                       }
  Float_t*  flavor()                      { SetBranch("flavor");                      return flavor_;                     }
  Float_t*  VtxProb()                     { SetBranch("VtxProb");                     return VtxProb_;                    }
  Float_t*  V3DCosPointingAngle()         { SetBranch("3DCosPointingAngle");          return V3DCosPointingAngle_;        }
  Float_t*  V3DPointingAngle()            { SetBranch("3DPointingAngle");             return V3DPointingAngle_;           }
  Float_t*  V2DCosPointingAngle()         { SetBranch("2DCosPointingAngle");          return V2DCosPointingAngle_;        }
  Float_t*  V2DPointingAngle()            { SetBranch("2DPointingAngle");             return V2DPointingAngle_;           }
  Float_t*  V3DDecayLengthSignificance()  { SetBranch("3DDecayLengthSignificance");   return V3DDecayLengthSignificance_; }
  Float_t*  V3DDecayLength()              { SetBranch("3DDecayLength");               return V3DDecayLength_;             }
  Float_t*  V2DDecayLengthSignificance()  { SetBranch("2DDecayLengthSignificance");   return V2DDecayLengthSignificance_; }
  Float_t*  V2DDecayLength()              { SetBranch("2DDecayLength");               return V2DDecayLength_;             }
  Float_t*  zDCASignificanceDaugther1()   { SetBranch("zDCASignificanceDaugther1");   return zDCASignificanceDaugther1_;  }
  Float_t*  xyDCASignificanceDaugther1()  { SetBranch("xyDCASignificanceDaugther1");  return xyDCASignificanceDaugther1_; }
  Bool_t*   HighPurityDaugther1()         { SetBranch("HighPurityDaugther1");         return HighPurityDaugther1_;        }
  Float_t*  NHitD1()                      { SetBranch("NHitD1");                      return NHitD1_;                     }
  Float_t*  pTD1()                        { SetBranch("pTD1");                        return pTD1_;                       }
  Float_t*  pTerrD1()                     { SetBranch("pTerrD1");                     return pTerrD1_;                    }
  Float_t*  EtaD1()                       { SetBranch("EtaD1");                       return EtaD1_;                      }
  Float_t*  PhiD1()                       { SetBranch("PhiD1");                       return PhiD1_;                      }
  Short_t*  chargeD1()                    { SetBranch("chargeD1");                    return chargeD1_;                   }
  Float_t*  dedxHarmonic2D1()             { SetBranch("dedxHarmonic2D1");             return dedxHarmonic2D1_;            }
  Float_t*  zDCASignificanceDaugther2()   { SetBranch("zDCASignificanceDaugther2");   return zDCASignificanceDaugther2_;  }
  Float_t*  xyDCASignificanceDaugther2()  { SetBranch("xyDCASignificanceDaugther2");  return xyDCASignificanceDaugther2_; }
  Bool_t*   HighPurityDaugther2()         { SetBranch("HighPurityDaugther2");         return HighPurityDaugther2_;        }
  Float_t*  NHitD2()                      { SetBranch("NHitD2");                      return NHitD2_;                     }
  Float_t*  pTD2()                        { SetBranch("pTD2");                        return pTD2_;                       }
  Float_t*  pTerrD2()                     { SetBranch("pTerrD2");                     return pTerrD2_;                    }
  Float_t*  EtaD2()                       { SetBranch("EtaD2");                       return EtaD2_;                      }
  Float_t*  PhiD2()                       { SetBranch("PhiD2");                       return PhiD2_;                      }
  Short_t*  chargeD2()                    { SetBranch("chargeD2");                    return chargeD2_;                   }
  Float_t*  dedxHarmonic2D2()             { SetBranch("dedxHarmonic2D2");             return dedxHarmonic2D2_;            }
  Bool_t*   isSwap()                      { SetBranch("isSwap");                      return isSwap_;                     }
  Int_t*    idmom_reco()                  { SetBranch("idmom_reco");                  return idmom_reco_;                 }
  Bool_t*   matchGEN()                    { SetBranch("matchGEN");                    return matchGEN_;                   }
  Int_t*    PIDD1()                       { SetBranch("PIDD1");                       return PIDD1_;                      }
  Int_t*    PIDD2()                       { SetBranch("PIDD2");                       return PIDD2_;                      }

  // MUON INFO GETTERS
  Bool_t*   OneStMuon1()                  { SetBranch("OneStMuon1");                  return OneStMuon1_;                 }
  Bool_t*   PFMuon1()                     { SetBranch("PFMuon1");                     return PFMuon1_;                    }
  Bool_t*   GlbMuon1()                    { SetBranch("GlbMuon1");                    return GlbMuon1_;                   }
  Bool_t*   trkMuon1()                    { SetBranch("trkMuon1");                    return trkMuon1_;                   }
  Bool_t*   tightMuon1()                  { SetBranch("tightMuon1");                  return tightMuon1_;                 }
  Bool_t*   softMuon1()                   { SetBranch("softMuon1");                   return softMuon1_;                  }
  Bool_t*   hybridMuon1()                 { SetBranch("hybridMuon1");                 return hybridMuon1_;                }
  Bool_t*   HPMuon1()                     { SetBranch("HPMuon1");                     return HPMuon1_;                    }
  UCharVecVec trigMuon1()                 { SetBranch("trigMuon1");                   return GET(trigMuon1_);             }
  Short_t*  nMatchedStationD1()           { SetBranch("nMatchedStationD1");           return nMatchedStationD1_;          }
  Short_t*  nTrackerLayerD1()             { SetBranch("nTrackerLayerD1");             return nTrackerLayerD1_;            }
  Short_t*  nPixelLayerD1()               { SetBranch("nPixelLayerD1");               return nPixelLayerD1_;              }
  Short_t*  nPixelHitD1()                 { SetBranch("nPixelHitD1");                 return nPixelHitD1_;                }
  Short_t*  nMuonHitD1()                  { SetBranch("nMuonHitD1");                  return nMuonHitD1_;                 }
  Float_t*  GlbTrkChiD1()                 { SetBranch("GlbTrkChiD1");                 return GlbTrkChiD1_;                }
  Float_t*  muondXYD1()                   { SetBranch("muondXYD1");                   return muondXYD1_;                  }
  Float_t*  muondZD1()                    { SetBranch("muondZD1");                    return muondZD1_;                   }
  Float_t*  dXYD1()                       { SetBranch("dXYD1");                       return dXYD1_;                      }
  Float_t*  dZD1()                        { SetBranch("dZD1");                        return dZD1_;                       }
  Short_t*  nMatchedChamberD1()           { SetBranch("nMatchedChamberD1");           return nMatchedChamberD1_;          }
  Float_t*  EnergyDepositionD1()          { SetBranch("EnergyDepositionD1");          return EnergyDepositionD1_;         }
  Float_t*  dx1_seg()                     { SetBranch("dx1_seg");                     return dx1_seg_;                    }
  Float_t*  dy1_seg()                     { SetBranch("dy1_seg");                     return dy1_seg_;                    }
  Float_t*  dxSig1_seg()                  { SetBranch("dxSig1_seg");                  return dxSig1_seg_;                 }
  Float_t*  dySig1_seg()                  { SetBranch("dySig1_seg");                  return dySig1_seg_;                 }
  Float_t*  ddxdz1_seg()                  { SetBranch("ddxdz1_seg");                  return ddxdz1_seg_;                 }
  Float_t*  ddydz1_seg()                  { SetBranch("ddydz1_seg");                  return ddydz1_seg_;                 }
  Float_t*  ddxdzSig1_seg()               { SetBranch("ddxdzSig1_seg");               return ddxdzSig1_seg_;              }
  Float_t*  ddydzSig1_seg()               { SetBranch("ddydzSig1_seg");               return ddydzSig1_seg_;              }
  Bool_t*   OneStMuon2()                  { SetBranch("OneStMuon2");                  return OneStMuon2_;                 }
  Bool_t*   PFMuon2()                     { SetBranch("PFMuon2");                     return PFMuon2_;                    }
  Bool_t*   GlbMuon2()                    { SetBranch("GlbMuon2");                    return GlbMuon2_;                   }
  Bool_t*   trkMuon2()                    { SetBranch("trkMuon2");                    return trkMuon2_;                   }
  Bool_t*   tightMuon2()                  { SetBranch("tightMuon2");                  return tightMuon2_;                 }
  Bool_t*   softMuon2()                   { SetBranch("softMuon2");                   return softMuon2_;                  }
  Bool_t*   hybridMuon2()                 { SetBranch("hybridMuon2");                 return hybridMuon2_;                }
  Bool_t*   HPMuon2()                     { SetBranch("HPMuon2");                     return HPMuon2_;                    }
  UCharVecVec trigMuon2()                 { SetBranch("trigMuon2");                   return GET(trigMuon2_);             }
  Short_t*  nMatchedStationD2()           { SetBranch("nMatchedStationD2");           return nMatchedStationD2_;          }
  Short_t*  nTrackerLayerD2()             { SetBranch("nTrackerLayerD2");             return nTrackerLayerD2_;            }
  Short_t*  nPixelLayerD2()               { SetBranch("nPixelLayerD2");               return nPixelLayerD2_;              }
  Short_t*  nPixelHitD2()                 { SetBranch("nPixelHitD2");                 return nPixelHitD2_;                }
  Short_t*  nMuonHitD2()                  { SetBranch("nMuonHitD2");                  return nMuonHitD2_;                 }
  Float_t*  GlbTrkChiD2()                 { SetBranch("GlbTrkChiD2");                 return GlbTrkChiD2_;                }
  Float_t*  muondXYD2()                   { SetBranch("muondXYD2");                   return muondXYD2_;                  }
  Float_t*  muondZD2()                    { SetBranch("muondZD2");                    return muondZD2_;                   }
  Float_t*  dXYD2()                       { SetBranch("dXYD2");                       return dXYD2_;                      }
  Float_t*  dZD2()                        { SetBranch("dZD2");                        return dZD2_;                       }
  Short_t*  nMatchedChamberD2()           { SetBranch("nMatchedChamberD2");           return nMatchedChamberD2_;          }
  Float_t*  EnergyDepositionD2()          { SetBranch("EnergyDepositionD2");          return EnergyDepositionD2_;         }
  Float_t*  dx2_seg()                     { SetBranch("dx2_seg");                     return dx2_seg_;                    }
  Float_t*  dy2_seg()                     { SetBranch("dy2_seg");                     return dy2_seg_;                    }
  Float_t*  dxSig2_seg()                  { SetBranch("dxSig2_seg");                  return dxSig2_seg_;                 }
  Float_t*  dySig2_seg()                  { SetBranch("dySig2_seg");                  return dySig2_seg_;                 }
  Float_t*  ddxdz2_seg()                  { SetBranch("ddxdz2_seg");                  return ddxdz2_seg_;                 }
  Float_t*  ddydz2_seg()                  { SetBranch("ddydz2_seg");                  return ddydz2_seg_;                 }
  Float_t*  ddxdzSig2_seg()               { SetBranch("ddxdzSig2_seg");               return ddxdzSig2_seg_;              }
  Float_t*  ddydzSig2_seg()               { SetBranch("ddydzSig2_seg");               return ddydzSig2_seg_;              }

  // GEN INFO GETTERS
  Float_t   weight_gen()                  { SetBranch("weight_gen");                  return weight_gen_;                 }
  UInt_t    candSize_gen()                { SetBranch("candSize_gen");                return candSize_gen_;               }
  Float_t*  pT_gen()                      { SetBranch("pT_gen");                      return pT_gen_;                     }
  Float_t*  eta_gen()                     { SetBranch("eta_gen");                     return eta_gen_;                    }
  Float_t*  y_gen()                       { SetBranch("y_gen");                       return y_gen_;                      }
  Short_t*  status_gen()                  { SetBranch("status_gen");                  return status_gen_;                 }
  Int_t*    PID_gen()                     { SetBranch("PID_gen");                     return PID_gen_;                    }
  Short_t*    DecayID_gen()                 { SetBranch("DecayID_gen");                 return DecayID_gen_;                    }
  Int_t*    MotherID_gen()                { SetBranch("MotherID_gen");                return MotherID_gen_;               }
  Short_t*  RecIdx_gen()                  { SetBranch("RecIdx_gen");                  return RecIdx_gen_;                 }
  Int_t*    PIDD1_gen()                   { SetBranch("PIDD1_gen");                   return PIDD1_gen_;                  }
  Short_t*  chargeD1_gen()                { SetBranch("chargeD1_gen");                return chargeD1_gen_;               }
  Float_t*  pTD1_gen()                    { SetBranch("pTD1_gen");                    return pTD1_gen_;                   }
  Float_t*  EtaD1_gen()                   { SetBranch("EtaD1_gen");                   return EtaD1_gen_;                  }
  Float_t*  PhiD1_gen()                   { SetBranch("PhiD1_gen");                   return PhiD1_gen_;                  }
  Int_t*    PIDD2_gen()                   { SetBranch("PIDD2_gen");                   return PIDD2_gen_;                  }
  Short_t*  chargeD2_gen()                { SetBranch("chargeD2_gen");                return chargeD2_gen_;               }
  Float_t*  pTD2_gen()                    { SetBranch("pTD2_gen");                    return pTD2_gen_;                   }
  Float_t*  EtaD2_gen()                   { SetBranch("EtaD2_gen");                   return EtaD2_gen_;                  }
  Float_t*  PhiD2_gen()                   { SetBranch("PhiD2_gen");                   return PhiD2_gen_;                  }

  // EXTRA GETTERS
  Bool_t    tightMuon1  (const UInt_t& iC, const std::string& type="");
  Bool_t    tightMuon2  (const UInt_t& iC, const std::string& type="");
  Bool_t    hybridMuon1 (const UInt_t& iC, const std::string& type="");
  Bool_t    hybridMuon2 (const UInt_t& iC, const std::string& type="");
  Bool_t    softMuon1   (const UInt_t& iC, const std::string& type="");
  Bool_t    softMuon2   (const UInt_t& iC, const std::string& type="");
  Bool_t    tightCand   (const UInt_t& iC, const std::string& type="") { return (tightMuon1(iC, type) && tightMuon2(iC, type));   }
  Bool_t    hybridCand  (const UInt_t& iC, const std::string& type="") { return (hybridMuon1(iC, type) && hybridMuon2(iC, type)); }
  Bool_t    softCand    (const UInt_t& iC, const std::string& type="") { return (softMuon1(iC, type) && softMuon2(iC, type));     }
  Bool_t    trigCand    (const UInt_t& iT, const UInt_t& iC, const bool& OR=false) { if (trigMuon1().size()<=iT) { std::cout << "[ERROR] Trigger index1: "<<iT<<">"<<trigMuon1().size() << std::endl; return false; }; return (OR ? (trigMuon1()[iT][iC] || trigMuon2()[iT][iC]) : (trigMuon1()[iT][iC] && trigMuon2()[iT][iC])); }
  Double_t  phiAsym     (const UInt_t& iC);


 private:

  virtual Long64_t  LoadTree        (Long64_t);
  virtual char      GetBranchStatus (const std::string&);
  virtual void      SetBranch       (const std::string&);
  virtual void      InitTree        (void);
  virtual Int_t     LoadEntry       (void) { return fChain_->GetEntry(entry_); }
  virtual void      GenerateDictionaries (void);

  template <typename T> T GET(T* x) { return ( (x) ? *x : T() ); }
  
  std::map<std::string, TChain*> fChainM_;
  TChain*   fChain_; // DONT USE SMART POINTERS
  Int_t     fCurrent_=-1;
  Long64_t  entry_;
  
  static const UInt_t NEP   = 3;
  static const UInt_t NTRG  = 15;
  static const UInt_t NSEL  = 10;

  // EVENT INFO VARIABLES
  UInt_t            RunNb_=0;
  UInt_t            LSNb_=0;
  UInt_t            EventNb_=0;
  Short_t           nPV_=-1.;
  Float_t           bestvtxX_=-99.;
  Float_t           bestvtxY_=-99.;
  Float_t           bestvtxZ_=-99.;
  Short_t           centrality_=-1;
  Int_t             Npixel_=-1;
  Float_t           HFsumETPlus_=-1.;
  Float_t           HFsumETMinus_=-1.;
  Float_t           ZDCPlus_=-1.;
  Float_t           ZDCMinus_=-1.;
  Int_t             Ntrkoffline_=-1;
  Short_t           trigPrescale_[NTRG]={0};
  Bool_t            trigHLT_[NTRG]={0};
  Bool_t            evtSel_[NSEL]={0};

  // EVENT PLANE VARIABLES
  Float_t           ephfpSumW_=-99.;
  Float_t           ephfpAngle_[NEP]={0};
  Float_t           ephfpQ_[NEP]={0};
  Float_t           ephfmSumW_=-99.;
  Float_t           ephfmAngle_[NEP]={0};
  Float_t           ephfmQ_[NEP]={0};
  Float_t           eptrackmidSumW_=-99.;
  Float_t           eptrackmidAngle_[NEP]={0};
  Float_t           eptrackmidQ_[NEP]={0};

  // CANDIDATE INFO VARIABLES
  UInt_t            candSize_=0;
  Float_t           pT_[NCAND]={0};   //[candSize]
  Float_t           eta_[NCAND]={0};   //[candSize]
  Float_t           y_[NCAND]={0};   //[candSize]
  Float_t           phi_[NCAND]={0};   //[candSize]
  Float_t           mass_[NCAND]={0};   //[candSize]
  Float_t           flavor_[NCAND]={0};   //[candSize]
  Float_t           VtxProb_[NCAND]={0};   //[candSize]
  Float_t           V3DCosPointingAngle_[NCAND]={0};   //[candSize]
  Float_t           V3DPointingAngle_[NCAND]={0};   //[candSize]
  Float_t           V2DCosPointingAngle_[NCAND]={0};   //[candSize]
  Float_t           V2DPointingAngle_[NCAND]={0};   //[candSize]
  Float_t           V3DDecayLengthSignificance_[NCAND]={0};   //[candSize]
  Float_t           V3DDecayLength_[NCAND]={0};   //[candSize]
  Float_t           V2DDecayLengthSignificance_[NCAND]={0};   //[candSize]
  Float_t           V2DDecayLength_[NCAND]={0};   //[candSize]
  Float_t           zDCASignificanceDaugther1_[NCAND]={0};   //[candSize]
  Float_t           xyDCASignificanceDaugther1_[NCAND]={0};   //[candSize]
  Bool_t            HighPurityDaugther1_[NCAND]={0};   //[candSize]
  Float_t           NHitD1_[NCAND]={0};   //[candSize]
  Float_t           pTD1_[NCAND]={0};   //[candSize]
  Float_t           pTerrD1_[NCAND]={0};   //[candSize]
  Float_t           EtaD1_[NCAND]={0};   //[candSize]
  Float_t           PhiD1_[NCAND]={0};   //[candSize]
  Short_t           chargeD1_[NCAND]={0};   //[candSize]
  Float_t           dedxHarmonic2D1_[NCAND]={0};   //[candSize]
  Float_t           zDCASignificanceDaugther2_[NCAND]={0};   //[candSize]
  Float_t           xyDCASignificanceDaugther2_[NCAND]={0};   //[candSize]
  Bool_t            HighPurityDaugther2_[NCAND]={0};   //[candSize]
  Float_t           NHitD2_[NCAND]={0};   //[candSize]
  Float_t           pTD2_[NCAND]={0};   //[candSize]
  Float_t           pTerrD2_[NCAND]={0};   //[candSize]
  Float_t           EtaD2_[NCAND]={0};   //[candSize]
  Float_t           PhiD2_[NCAND]={0};   //[candSize]
  Short_t           chargeD2_[NCAND]={0};   //[candSize]
  Float_t           dedxHarmonic2D2_[NCAND]={0};   //[candSize]
  Bool_t            isSwap_[NCAND]={0};   //[candSize]
  Int_t             idmom_reco_[NCAND]={0};   //[candSize]
  Bool_t            matchGEN_[NCAND]={0};   //[candSize]
  Int_t             PIDD1_[NCAND]={0};   //[candSize]
  Int_t             PIDD2_[NCAND]={0};   //[candSize]

  // MUON INFO VARIABLES
  Bool_t            OneStMuon1_[NCAND]={0};   //[candSize]
  Bool_t            PFMuon1_[NCAND]={0};   //[candSize]
  Bool_t            GlbMuon1_[NCAND]={0};   //[candSize]
  Bool_t            trkMuon1_[NCAND]={0};   //[candSize]
  Bool_t            tightMuon1_[NCAND]={0};   //[candSize]
  Bool_t            softMuon1_[NCAND]={0};   //[candSize]
  Bool_t            hybridMuon1_[NCAND]={0};   //[candSize]
  Bool_t            HPMuon1_[NCAND]={0};   //[candSize]
  UCharVecVec*      trigMuon1_=0;
  Short_t           nMatchedStationD1_[NCAND]={0};   //[candSize]
  Short_t           nTrackerLayerD1_[NCAND]={0};   //[candSize]
  Short_t           nPixelLayerD1_[NCAND]={0};   //[candSize]
  Short_t           nPixelHitD1_[NCAND]={0};   //[candSize]
  Short_t           nMuonHitD1_[NCAND]={0};   //[candSize]
  Float_t           GlbTrkChiD1_[NCAND]={0};   //[candSize]
  Float_t           muondXYD1_[NCAND]={0};   //[candSize]
  Float_t           muondZD1_[NCAND]={0};   //[candSize]
  Float_t           dXYD1_[NCAND]={0};   //[candSize]
  Float_t           dZD1_[NCAND]={0};   //[candSize]
  Short_t           nMatchedChamberD1_[NCAND]={0};   //[candSize]
  Float_t           EnergyDepositionD1_[NCAND]={0};   //[candSize]
  Float_t           dx1_seg_[NCAND]={0};   //[candSize]
  Float_t           dy1_seg_[NCAND]={0};   //[candSize]
  Float_t           dxSig1_seg_[NCAND]={0};   //[candSize]
  Float_t           dySig1_seg_[NCAND]={0};   //[candSize]
  Float_t           ddxdz1_seg_[NCAND]={0};   //[candSize]
  Float_t           ddydz1_seg_[NCAND]={0};   //[candSize]
  Float_t           ddxdzSig1_seg_[NCAND]={0};   //[candSize]
  Float_t           ddydzSig1_seg_[NCAND]={0};   //[candSize]
  Bool_t            OneStMuon2_[NCAND]={0};   //[candSize]
  Bool_t            PFMuon2_[NCAND]={0};   //[candSize]
  Bool_t            GlbMuon2_[NCAND]={0};   //[candSize]
  Bool_t            trkMuon2_[NCAND]={0};   //[candSize]
  Bool_t            tightMuon2_[NCAND]={0};   //[candSize]
  Bool_t            softMuon2_[NCAND]={0};   //[candSize]
  Bool_t            hybridMuon2_[NCAND]={0};   //[candSize]
  Bool_t            HPMuon2_[NCAND]={0};   //[candSize]
  UCharVecVec*      trigMuon2_=0;
  Short_t           nMatchedStationD2_[NCAND]={0};   //[candSize]
  Short_t           nTrackerLayerD2_[NCAND]={0};   //[candSize]
  Short_t           nPixelLayerD2_[NCAND]={0};   //[candSize]
  Short_t           nPixelHitD2_[NCAND]={0};   //[candSize]
  Short_t           nMuonHitD2_[NCAND]={0};   //[candSize]
  Float_t           GlbTrkChiD2_[NCAND]={0};   //[candSize]
  Float_t           muondXYD2_[NCAND]={0};   //[candSize]
  Float_t           muondZD2_[NCAND]={0};   //[candSize]
  Float_t           dXYD2_[NCAND]={0};   //[candSize]
  Float_t           dZD2_[NCAND]={0};   //[candSize]
  Short_t           nMatchedChamberD2_[NCAND]={0};   //[candSize]
  Float_t           EnergyDepositionD2_[NCAND]={0};   //[candSize]
  Float_t           dx2_seg_[NCAND]={0};   //[candSize]
  Float_t           dy2_seg_[NCAND]={0};   //[candSize]
  Float_t           dxSig2_seg_[NCAND]={0};   //[candSize]
  Float_t           dySig2_seg_[NCAND]={0};   //[candSize]
  Float_t           ddxdz2_seg_[NCAND]={0};   //[candSize]
  Float_t           ddydz2_seg_[NCAND]={0};   //[candSize]
  Float_t           ddxdzSig2_seg_[NCAND]={0};   //[candSize]
  Float_t           ddydzSig2_seg_[NCAND]={0};   //[candSize]

  // GEN INFO VARIABLES
  Float_t           weight_gen_=-99.;
  UInt_t            candSize_gen_=0;
  Float_t           pT_gen_[NGEN]={0};   //[candSize_gen]
  Float_t           eta_gen_[NGEN]={0};   //[candSize_gen]
  Float_t           y_gen_[NGEN]={0};   //[candSize_gen]
  Short_t           status_gen_[NGEN]={0};   //[candSize_gen]
  Int_t             PID_gen_[NGEN]={0};   //[candSize_gen]
  Short_t             DecayID_gen_[NGEN]={0};   //[candSize_gen]
  Int_t             MotherID_gen_[NGEN]={0};   //[candSize_gen]
  Short_t           RecIdx_gen_[NGEN]={0};   //[candSize_gen]
  Int_t             PIDD1_gen_[NGEN]={0};   //[candSize_gen]
  Short_t           chargeD1_gen_[NGEN]={0};   //[candSize_gen]
  Float_t           pTD1_gen_[NGEN]={0};   //[candSize_gen]
  Float_t           EtaD1_gen_[NGEN]={0};   //[candSize_gen]
  Float_t           PhiD1_gen_[NGEN]={0};   //[candSize_gen]
  Int_t             PIDD2_gen_[NGEN]={0};   //[candSize_gen]
  Short_t           chargeD2_gen_[NGEN]={0};   //[candSize_gen]
  Float_t           pTD2_gen_[NGEN]={0};   //[candSize_gen]
  Float_t           EtaD2_gen_[NGEN]={0};   //[candSize_gen]
  Float_t           PhiD2_gen_[NGEN]={0};   //[candSize_gen]

  // EVENT INFO BRANCHES
  TBranch          *b_RunNb;   //!
  TBranch          *b_LSNb;   //!
  TBranch          *b_EventNb;   //!
  TBranch          *b_nPV;   //!
  TBranch          *b_bestvtxX;   //!
  TBranch          *b_bestvtxY;   //!
  TBranch          *b_bestvtxZ;   //!
  TBranch          *b_centrality;   //!
  TBranch          *b_Npixel;   //!
  TBranch          *b_HFsumETPlus;   //!
  TBranch          *b_HFsumETMinus;   //!
  TBranch          *b_ZDCPlus;   //!
  TBranch          *b_ZDCMinus;   //!
  TBranch          *b_Ntrkoffline;   //!
  TBranch          *b_trigPrescale;   //!
  TBranch          *b_trigHLT;   //!
  TBranch          *b_evtSel;   //!

  // EVENT PLANE BRANCHES
  TBranch          *b_ephfpAngle;   //!
  TBranch          *b_ephfpQ;   //!
  TBranch          *b_ephfpSumW;   //!
  TBranch          *b_ephfmAngle;   //!
  TBranch          *b_ephfmQ;   //!
  TBranch          *b_ephfmSumW;   //!
  TBranch          *b_eptrackmidAngle;   //!
  TBranch          *b_eptrackmidQ;   //!
  TBranch          *b_eptrackmidSumW;   //!

  // CANDIDATE INFO BRANCHES
  TBranch          *b_candSize;   //!
  TBranch          *b_pT;   //!
  TBranch          *b_eta;   //!
  TBranch          *b_y;   //!
  TBranch          *b_phi;   //!
  TBranch          *b_mass;   //!
  TBranch          *b_flavor;   //!
  TBranch          *b_VtxProb;   //!
  TBranch          *b_3DCosPointingAngle;   //!
  TBranch          *b_3DPointingAngle;   //!
  TBranch          *b_2DCosPointingAngle;   //!
  TBranch          *b_2DPointingAngle;   //!
  TBranch          *b_3DDecayLengthSignificance;   //!
  TBranch          *b_3DDecayLength;   //!
  TBranch          *b_2DDecayLengthSignificance;   //!
  TBranch          *b_2DDecayLength;   //!
  TBranch          *b_zDCASignificanceDaugther1;   //!
  TBranch          *b_xyDCASignificanceDaugther1;   //!
  TBranch          *b_HighPurityDaugther1;   //!
  TBranch          *b_NHitD1;   //!
  TBranch          *b_pTD1;   //!
  TBranch          *b_pTerrD1;   //!
  TBranch          *b_EtaD1;   //!
  TBranch          *b_PhiD1;   //!
  TBranch          *b_chargeD1;   //!
  TBranch          *b_dedxHarmonic2D1;   //!
  TBranch          *b_zDCASignificanceDaugther2;   //!
  TBranch          *b_xyDCASignificanceDaugther2;   //!
  TBranch          *b_HighPurityDaugther2;   //!
  TBranch          *b_NHitD2;   //!
  TBranch          *b_pTD2;   //!
  TBranch          *b_pTerrD2;   //!
  TBranch          *b_EtaD2;   //!
  TBranch          *b_PhiD2;   //!
  TBranch          *b_chargeD2;   //!
  TBranch          *b_dedxHarmonic2D2;   //!
  TBranch          *b_isSwap;   //!
  TBranch          *b_idmom_reco;   //!
  TBranch          *b_matchGEN;   //!
  TBranch          *b_PIDD1;   //!
  TBranch          *b_PIDD2;   //!

  // MUON INFO BRANCHES
  TBranch          *b_OneStMuon1;   //!
  TBranch          *b_PFMuon1;   //!
  TBranch          *b_GlbMuon1;   //!
  TBranch          *b_trkMuon1;   //!
  TBranch          *b_tightMuon1;   //!
  TBranch          *b_softMuon1;   //!
  TBranch          *b_hybridMuon1;   //!
  TBranch          *b_HPMuon1;   //!
  TBranch          *b_trigMuon1;   //!
  TBranch          *b_nMatchedStationD1;   //!
  TBranch          *b_nTrackerLayerD1;   //!
  TBranch          *b_nPixelLayerD1;   //!
  TBranch          *b_nPixelHitD1;   //!
  TBranch          *b_nMuonHitD1;   //!
  TBranch          *b_GlbTrkChiD1;   //!
  TBranch          *b_muondXYD1;   //!
  TBranch          *b_muondZD1;   //!
  TBranch          *b_dXYD1;   //!
  TBranch          *b_dZD1;   //!
  TBranch          *b_nMatchedChamberD1;   //!
  TBranch          *b_EnergyDepositionD1;   //!
  TBranch          *b_dx1_seg;   //!
  TBranch          *b_dy1_seg;   //!
  TBranch          *b_dxSig1_seg;   //!
  TBranch          *b_dySig1_seg;   //!
  TBranch          *b_ddxdz1_seg;   //!
  TBranch          *b_ddydz1_seg;   //!
  TBranch          *b_ddxdzSig1_seg;   //!
  TBranch          *b_ddydzSig1_seg;   //!
  TBranch          *b_OneStMuon2;   //!
  TBranch          *b_PFMuon2;   //!
  TBranch          *b_GlbMuon2;   //!
  TBranch          *b_trkMuon2;   //!
  TBranch          *b_tightMuon2;   //!
  TBranch          *b_softMuon2;   //!
  TBranch          *b_hybridMuon2;   //!
  TBranch          *b_HPMuon2;   //!
  TBranch          *b_trigMuon2;   //!
  TBranch          *b_nMatchedStationD2;   //!
  TBranch          *b_nTrackerLayerD2;   //!
  TBranch          *b_nPixelLayerD2;   //!
  TBranch          *b_nPixelHitD2;   //!
  TBranch          *b_nMuonHitD2;   //!
  TBranch          *b_GlbTrkChiD2;   //!
  TBranch          *b_muondXYD2;   //!
  TBranch          *b_muondZD2;   //!
  TBranch          *b_dXYD2;   //!
  TBranch          *b_dZD2;   //!
  TBranch          *b_nMatchedChamberD2;   //!
  TBranch          *b_EnergyDepositionD2;   //!
  TBranch          *b_dx2_seg;   //!
  TBranch          *b_dy2_seg;   //!
  TBranch          *b_dxSig2_seg;   //!
  TBranch          *b_dySig2_seg;   //!
  TBranch          *b_ddxdz2_seg;   //!
  TBranch          *b_ddydz2_seg;   //!
  TBranch          *b_ddxdzSig2_seg;   //!
  TBranch          *b_ddydzSig2_seg;   //!
  
  // GEN INFO BRANCHES
  TBranch          *b_weight_gen;   //!
  TBranch          *b_candSize_gen;   //!
  TBranch          *b_pT_gen;   //!
  TBranch          *b_eta_gen;   //!
  TBranch          *b_y_gen;   //!
  TBranch          *b_status_gen;   //!
  TBranch          *b_PID_gen;   //!
  TBranch          *b_DecayID_gen;   //!
  TBranch          *b_MotherID_gen;   //!
  TBranch          *b_RecIdx_gen;   //!
  TBranch          *b_PIDD1_gen;   //!
  TBranch          *b_chargeD1_gen;   //!
  TBranch          *b_pTD1_gen;   //!
  TBranch          *b_EtaD1_gen;   //!
  TBranch          *b_PhiD1_gen;   //!
  TBranch          *b_PIDD2_gen;   //!
  TBranch          *b_chargeD2_gen;   //!
  TBranch          *b_pTD2_gen;   //!
  TBranch          *b_EtaD2_gen;   //!
  TBranch          *b_PhiD2_gen;   //!
};

VertexCompositeNtuple::VertexCompositeNtuple() : fChain_(0)
{
};

VertexCompositeNtuple::~VertexCompositeNtuple()
{
  if (fChain_) { const auto& f = fChain_->GetCurrentFile(); if (f) { f->Close(); delete f; }; fChain_->Reset(); }
  for (auto& c : fChainM_) { if (c.second) { c.second->Reset(); } }
};

Bool_t VertexCompositeNtuple::GetTree(const std::string& fileName, const std::string& treeName)
{
  const auto& fileNames = std::vector<std::string>({fileName});
  return GetTree(fileNames, treeName);
};

Bool_t VertexCompositeNtuple::GetTree(const std::vector< std::string >& inFileName, const std::string& treeName)
{
  // Check the File Names
  auto fileName = inFileName;
  for (auto& f : fileName) { if (f.rfind("/store/", 0)==0) { f = "root://cms-xrd-global.cern.ch/" + f; } }
  // Open the input files
  const auto& f = TFile::Open(fileName[0].c_str(), "READ");
  if (!f || !f->IsOpen() || f->IsZombie()) { std::cout << "[ERROR] Failed to open file: " << fileName[0] << std::endl; return false; }
  // Extract the input TChains
  std::cout << "[INFO] Extracting tree: " << treeName.c_str() << std::endl;
  fChainM_.clear();
  TDirectory * dir;
  if(treeName.size()==0){
    if(f->GetListOfKeys()->Contains("VertexCompositeNtuple")){
      fChainM_["VertexCompositeNtuple"] = new TChain("VertexCompositeNtuple" , "VertexCompositeNtuple");
    } 
    if (fChainM_.size()==0) { std::cout << "[ERROR] fChain VertexCompositeNtuple was not created, some input files are   missing" << std::endl; return false; }
  }
  else{
    if (fileName[0].rfind("root://", 0)==0) { dir = (TDirectory*)f->Get(treeName.c_str()); }
    else { dir = (TDirectory*)f->Get((fileName[0]+":/"+treeName).c_str()); }
    if (!dir) { std::cout << "[ERROR] Failed to open directory: " << treeName << std::endl; return false; }
    if (dir->GetListOfKeys()->Contains("VertexCompositeNtuple")) { fChainM_["VertexCompositeNtuple"] = new TChain((treeName+"/VertexCompositeNtuple").c_str() , "VertexCompositeNtuple"); }
    if (fChainM_.size()==0) { std::cout << "[ERROR] fChain VertexCompositeNtuple was not created, some input files are   missing" << std::endl; return false; }
  } 

  // Add the files in the TChain
  for (auto& c : fChainM_) {
    if(treeName.size()==0){
      for (auto& f : fileName) { c.second->Add(Form("%s/%s", f.c_str(), c.first.c_str())); }; c.second->GetEntries();
    }
    else{
      for (auto& f : fileName) { c.second->Add(Form("%s/%s/%s", f.c_str(), treeName.c_str(), c.first.c_str())); }; c.second->GetEntries();
    }
  }
  for (auto& c : fChainM_) { if (!c.second) { std::cout << "[ERROR] fChain " << c.first << " was not created, some input files are missing" << std::endl; return false; } }
  // Initialize the input TChains (set their branches)
  InitTree();
  // Add Friend TChains
  if (fChain_) { delete fChain_; }
  fChain_ = (TChain*)fChainM_.begin()->second->Clone(treeName.c_str());
  for (auto& c : fChainM_) {
    if (c.second != fChain_) {
      c.second->SetMakeClass(1); // For the proper setup.
      fChain_->AddFriend(c.second, c.first.c_str(), kTRUE); // Add the Friend TChain
    }
  }
  if (!fChain_) return false;
  // Set All Branches to Status 0
  fChain_->SetBranchStatus("*",0);
  //
  return true;
};

Int_t VertexCompositeNtuple::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  entry_ = entry;
  if (LoadTree(entry_) < 0) return -1;
  Clear();
  const auto& status = LoadEntry();
  // Check contents of entry
  if (candSize_ >= NCAND) { std::cout << "[ERROR] Reconstructed candidate size ("<<candSize_<<") is larger than "<<NCAND << std::endl; return -9; }
  if (candSize_gen_ >= NGEN) { std::cout << "[ERROR] Generated candidate size ("<<candSize_gen_<<") is larger than "<<NGEN << std::endl; return -9; }
  return status;
};

Long64_t VertexCompositeNtuple::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain_) return -5;
  const auto& centry = fChain_->LoadTree(entry);
  if (fChain_->GetTreeNumber() != fCurrent_) { fCurrent_ = fChain_->GetTreeNumber(); }
  return centry;
};

char VertexCompositeNtuple::GetBranchStatus(const std::string& n)
{ 
  if (!fChain_ || !(fChain_->GetBranch(n.c_str()))) return -1;
  return fChain_->GetBranchStatus(n.c_str());
};

void VertexCompositeNtuple::SetBranch(const std::string& n)
{
  if (GetBranchStatus(n) == 0) {
    fChain_->SetBranchStatus(Form("*%s*", n.c_str()), 1);
    LoadEntry(); // Needed for the first entry
  }
};

void VertexCompositeNtuple::InitTree(void)
{
  // Generate the dictionary's needed
  GenerateDictionaries();
  
  // Initialize pointers
  trigMuon1_ = 0;
  trigMuon2_ = 0;
  
  if (fChainM_.count("VertexCompositeNtuple")>0) {
    auto& fChain = fChainM_.at("VertexCompositeNtuple");

    // SET EVENT INFO BRANCHES
    if (fChain->GetBranch("RunNb"))                       fChain->SetBranchAddress("RunNb",                     &RunNb_,                      &b_RunNb                     );
    if (fChain->GetBranch("LSNb"))                        fChain->SetBranchAddress("LSNb",                      &LSNb_,                       &b_LSNb                      );
    if (fChain->GetBranch("EventNb"))                     fChain->SetBranchAddress("EventNb",                   &EventNb_,                    &b_EventNb                   );
    if (fChain->GetBranch("nPV"))                         fChain->SetBranchAddress("nPV",                       &nPV_,                        &b_nPV                       );
    if (fChain->GetBranch("bestvtxX"))                    fChain->SetBranchAddress("bestvtxX",                  &bestvtxX_,                   &b_bestvtxX                  );
    if (fChain->GetBranch("bestvtxY"))                    fChain->SetBranchAddress("bestvtxY",                  &bestvtxY_,                   &b_bestvtxY                  );
    if (fChain->GetBranch("bestvtxZ"))                    fChain->SetBranchAddress("bestvtxZ",                  &bestvtxZ_,                   &b_bestvtxZ                  );
    if (fChain->GetBranch("centrality"))                  fChain->SetBranchAddress("centrality",                &centrality_,                 &b_centrality                );
    if (fChain->GetBranch("Npixel"))                      fChain->SetBranchAddress("Npixel",                    &Npixel_,                     &b_Npixel                    );
    if (fChain->GetBranch("HFsumETPlus"))                 fChain->SetBranchAddress("HFsumETPlus",               &HFsumETPlus_,                &b_HFsumETPlus               );
    if (fChain->GetBranch("HFsumETMinus"))                fChain->SetBranchAddress("HFsumETMinus",              &HFsumETMinus_,               &b_HFsumETMinus              );
    if (fChain->GetBranch("ZDCPlus"))                     fChain->SetBranchAddress("ZDCPlus",                   &ZDCPlus_,                    &b_ZDCPlus                   );
    if (fChain->GetBranch("ZDCMinus"))                    fChain->SetBranchAddress("ZDCMinus",                  &ZDCMinus_,                   &b_ZDCMinus                  );
    if (fChain->GetBranch("Ntrkoffline"))                 fChain->SetBranchAddress("Ntrkoffline",               &Ntrkoffline_,                &b_Ntrkoffline               );
    if (fChain->GetBranch("trigPrescale"))                fChain->SetBranchAddress("trigPrescale",               trigPrescale_,               &b_trigPrescale              );
    if (fChain->GetBranch("trigHLT"))                     fChain->SetBranchAddress("trigHLT",                    trigHLT_,                    &b_trigHLT                   );
    if (fChain->GetBranch("evtSel"))                      fChain->SetBranchAddress("evtSel",                     evtSel_,                     &b_evtSel                    );

    // SET EVENT PLANE BRANCHES
    if (fChain->GetBranch("ephfpSumW"))                   fChain->SetBranchAddress("ephfpSumW",                 &ephfpSumW_,                  &b_ephfpSumW                 );
    if (fChain->GetBranch("ephfpAngle"))                  fChain->SetBranchAddress("ephfpAngle",                 ephfpAngle_,                 &b_ephfpAngle                );
    if (fChain->GetBranch("ephfpQ"))                      fChain->SetBranchAddress("ephfpQ",                     ephfpQ_,                     &b_ephfpQ                    );
    if (fChain->GetBranch("ephfmSumW"))                   fChain->SetBranchAddress("ephfmSumW",                 &ephfmSumW_,                  &b_ephfmSumW                 );
    if (fChain->GetBranch("ephfmAngle"))                  fChain->SetBranchAddress("ephfmAngle",                 ephfmAngle_,                 &b_ephfmAngle                );
    if (fChain->GetBranch("ephfmQ"))                      fChain->SetBranchAddress("ephfmQ",                     ephfmQ_,                     &b_ephfmQ                    );
    if (fChain->GetBranch("eptrackmidSumW"))                   fChain->SetBranchAddress("eptrackmidSumW",                 &eptrackmidSumW_,                  &b_eptrackmidSumW                 );
    if (fChain->GetBranch("eptrackmidAngle"))                  fChain->SetBranchAddress("eptrackmidAngle",                 eptrackmidAngle_,                 &b_eptrackmidAngle                );
    if (fChain->GetBranch("eptrackmidQ"))                      fChain->SetBranchAddress("eptrackmidQ",                     eptrackmidQ_,                     &b_eptrackmidQ                    );

    // SET CANDIDATE INFO BRANCHES
    if (fChain->GetBranch("candSize"))                    fChain->SetBranchAddress("candSize",                  &candSize_,                   &b_candSize                  );
    if (fChain->GetBranch("pT"))                          fChain->SetBranchAddress("pT",                         pT_,                         &b_pT                        );
    if (fChain->GetBranch("eta"))                         fChain->SetBranchAddress("eta",                        eta_,                        &b_eta                       );
    if (fChain->GetBranch("y"))                           fChain->SetBranchAddress("y",                          y_,                          &b_y                         );
    if (fChain->GetBranch("phi"))                         fChain->SetBranchAddress("phi",                        phi_,                        &b_phi                       );
    if (fChain->GetBranch("mass"))                        fChain->SetBranchAddress("mass",                       mass_,                       &b_mass                      );
    if (fChain->GetBranch("flavor"))                      fChain->SetBranchAddress("flavor",                     flavor_,                     &b_flavor                    );
    if (fChain->GetBranch("VtxProb"))                     fChain->SetBranchAddress("VtxProb",                    VtxProb_,                    &b_VtxProb                   );
    if (fChain->GetBranch("3DCosPointingAngle"))          fChain->SetBranchAddress("3DCosPointingAngle",         V3DCosPointingAngle_,        &b_3DCosPointingAngle        );
    if (fChain->GetBranch("3DPointingAngle"))             fChain->SetBranchAddress("3DPointingAngle",            V3DPointingAngle_,           &b_3DPointingAngle           );
    if (fChain->GetBranch("2DCosPointingAngle"))          fChain->SetBranchAddress("2DCosPointingAngle",         V2DCosPointingAngle_,        &b_2DCosPointingAngle        );
    if (fChain->GetBranch("2DPointingAngle"))             fChain->SetBranchAddress("2DPointingAngle",            V2DPointingAngle_,           &b_2DPointingAngle           );
    if (fChain->GetBranch("3DDecayLengthSignificance"))   fChain->SetBranchAddress("3DDecayLengthSignificance",  V3DDecayLengthSignificance_, &b_3DDecayLengthSignificance );
    if (fChain->GetBranch("3DDecayLength"))               fChain->SetBranchAddress("3DDecayLength",              V3DDecayLength_,             &b_3DDecayLength             );
    if (fChain->GetBranch("2DDecayLengthSignificance"))   fChain->SetBranchAddress("2DDecayLengthSignificance",  V2DDecayLengthSignificance_, &b_2DDecayLengthSignificance );
    if (fChain->GetBranch("2DDecayLength"))               fChain->SetBranchAddress("2DDecayLength",              V2DDecayLength_,             &b_2DDecayLength             );
    if (fChain->GetBranch("zDCASignificanceDaugther1"))   fChain->SetBranchAddress("zDCASignificanceDaugther1",  zDCASignificanceDaugther1_,  &b_zDCASignificanceDaugther1 );
    if (fChain->GetBranch("xyDCASignificanceDaugther1"))  fChain->SetBranchAddress("xyDCASignificanceDaugther1", xyDCASignificanceDaugther1_, &b_xyDCASignificanceDaugther1);
    if (fChain->GetBranch("HighPurityDaugther1"))         fChain->SetBranchAddress("HighPurityDaugther1",        HighPurityDaugther1_,        &b_HighPurityDaugther1       );
    if (fChain->GetBranch("NHitD1"))                      fChain->SetBranchAddress("NHitD1",                     NHitD1_,                     &b_NHitD1                    );
    if (fChain->GetBranch("pTD1"))                        fChain->SetBranchAddress("pTD1",                       pTD1_,                       &b_pTD1                      );
    if (fChain->GetBranch("pTerrD1"))                     fChain->SetBranchAddress("pTerrD1",                    pTerrD1_,                    &b_pTerrD1                   );
    if (fChain->GetBranch("EtaD1"))                       fChain->SetBranchAddress("EtaD1",                      EtaD1_,                      &b_EtaD1                     );
    if (fChain->GetBranch("PhiD1"))                       fChain->SetBranchAddress("PhiD1",                      PhiD1_,                      &b_PhiD1                     );
    if (fChain->GetBranch("chargeD1"))                    fChain->SetBranchAddress("chargeD1",                   chargeD1_,                   &b_chargeD1                  );
    if (fChain->GetBranch("dedxHarmonic2D1"))             fChain->SetBranchAddress("dedxHarmonic2D1",            dedxHarmonic2D1_,            &b_dedxHarmonic2D1           );
    if (fChain->GetBranch("zDCASignificanceDaugther2"))   fChain->SetBranchAddress("zDCASignificanceDaugther2",  zDCASignificanceDaugther2_,  &b_zDCASignificanceDaugther2 );
    if (fChain->GetBranch("xyDCASignificanceDaugther2"))  fChain->SetBranchAddress("xyDCASignificanceDaugther2", xyDCASignificanceDaugther2_, &b_xyDCASignificanceDaugther2);
    if (fChain->GetBranch("HighPurityDaugther2"))         fChain->SetBranchAddress("HighPurityDaugther2",        HighPurityDaugther2_,        &b_HighPurityDaugther2       );
    if (fChain->GetBranch("NHitD2"))                      fChain->SetBranchAddress("NHitD2",                     NHitD2_,                     &b_NHitD2                    );
    if (fChain->GetBranch("pTD2"))                        fChain->SetBranchAddress("pTD2",                       pTD2_,                       &b_pTD2                      );
    if (fChain->GetBranch("pTerrD2"))                     fChain->SetBranchAddress("pTerrD2",                    pTerrD2_,                    &b_pTerrD2                   );
    if (fChain->GetBranch("EtaD2"))                       fChain->SetBranchAddress("EtaD2",                      EtaD2_,                      &b_EtaD2                     );
    if (fChain->GetBranch("PhiD2"))                       fChain->SetBranchAddress("PhiD2",                      PhiD2_,                      &b_PhiD2                     );
    if (fChain->GetBranch("chargeD2"))                    fChain->SetBranchAddress("chargeD2",                   chargeD2_,                   &b_chargeD2                  );
    if (fChain->GetBranch("dedxHarmonic2D2"))             fChain->SetBranchAddress("dedxHarmonic2D2",            dedxHarmonic2D2_,            &b_dedxHarmonic2D2           );
    if (fChain->GetBranch("isSwap"))                      fChain->SetBranchAddress("isSwap",                     isSwap_,                     &b_isSwap                    );
    if (fChain->GetBranch("idmom_reco"))                  fChain->SetBranchAddress("idmom_reco",                 idmom_reco_,                 &b_idmom_reco                );
    if (fChain->GetBranch("matchGEN"))                    fChain->SetBranchAddress("matchGEN",                   matchGEN_,                   &b_matchGEN                  );
    if (fChain->GetBranch("PIDD1"))                       fChain->SetBranchAddress("PIDD1",                      PIDD1_,                      &b_PIDD1                     );
    if (fChain->GetBranch("PIDD2"))                       fChain->SetBranchAddress("PIDD2",                      PIDD2_,                      &b_PIDD2                     );

    // SET MUON INFO BRANCHES
    if (fChain->GetBranch("OneStMuon1"))                  fChain->SetBranchAddress("OneStMuon1",                 OneStMuon1_,                 &b_OneStMuon1                );
    if (fChain->GetBranch("PFMuon1"))                     fChain->SetBranchAddress("PFMuon1",                    PFMuon1_,                    &b_PFMuon1                   );
    if (fChain->GetBranch("GlbMuon1"))                    fChain->SetBranchAddress("GlbMuon1",                   GlbMuon1_,                   &b_GlbMuon1                  );
    if (fChain->GetBranch("trkMuon1"))                    fChain->SetBranchAddress("trkMuon1",                   trkMuon1_,                   &b_trkMuon1                  );
    if (fChain->GetBranch("tightMuon1"))                  fChain->SetBranchAddress("tightMuon1",                 tightMuon1_,                 &b_tightMuon1                );
    if (fChain->GetBranch("softMuon1"))                   fChain->SetBranchAddress("softMuon1",                  softMuon1_,                  &b_softMuon1                 );
    if (fChain->GetBranch("hybridMuon1"))                 fChain->SetBranchAddress("hybridMuon1",                hybridMuon1_,                &b_hybridMuon1               );
    if (fChain->GetBranch("HPMuon1"))                     fChain->SetBranchAddress("HPMuon1",                    HPMuon1_,                    &b_HPMuon1                   );
    if (fChain->GetBranch("trigMuon1"))                   fChain->SetBranchAddress("trigMuon1",                 &trigMuon1_,                  &b_trigMuon1                 );
    if (fChain->GetBranch("nMatchedStationD1"))           fChain->SetBranchAddress("nMatchedStationD1",          nMatchedStationD1_,          &b_nMatchedStationD1         );
    if (fChain->GetBranch("nTrackerLayerD1"))             fChain->SetBranchAddress("nTrackerLayerD1",            nTrackerLayerD1_,            &b_nTrackerLayerD1           );
    if (fChain->GetBranch("nPixelLayerD1"))               fChain->SetBranchAddress("nPixelLayerD1",              nPixelLayerD1_,              &b_nPixelLayerD1             );
    if (fChain->GetBranch("nPixelHitD1"))                 fChain->SetBranchAddress("nPixelHitD1",                nPixelHitD1_,                &b_nPixelHitD1               );
    if (fChain->GetBranch("nMuonHitD1"))                  fChain->SetBranchAddress("nMuonHitD1",                 nMuonHitD1_,                 &b_nMuonHitD1                );
    if (fChain->GetBranch("GlbTrkChiD1"))                 fChain->SetBranchAddress("GlbTrkChiD1",                GlbTrkChiD1_,                &b_GlbTrkChiD1               );
    if (fChain->GetBranch("muondZD1"))                    fChain->SetBranchAddress("muondZD1",                   muondZD1_,                   &b_muondZD1                  );
    if (fChain->GetBranch("muondXYD1"))                   fChain->SetBranchAddress("muondXYD1",                  muondXYD1_,                  &b_muondXYD1                 );
    if (fChain->GetBranch("dZD1"))                        fChain->SetBranchAddress("dZD1",                       dZD1_,                       &b_dZD1                      );
    if (fChain->GetBranch("dXYD1"))                       fChain->SetBranchAddress("dXYD1",                      dXYD1_,                      &b_dXYD1                     );    
    if (fChain->GetBranch("OneStMuon2"))                  fChain->SetBranchAddress("OneStMuon2",                 OneStMuon2_,                 &b_OneStMuon2                );
    if (fChain->GetBranch("PFMuon2"))                     fChain->SetBranchAddress("PFMuon2",                    PFMuon2_,                    &b_PFMuon2                   );
    if (fChain->GetBranch("GlbMuon2"))                    fChain->SetBranchAddress("GlbMuon2",                   GlbMuon2_,                   &b_GlbMuon2                  );
    if (fChain->GetBranch("trkMuon2"))                    fChain->SetBranchAddress("trkMuon2",                   trkMuon2_,                   &b_trkMuon2                  );
    if (fChain->GetBranch("tightMuon2"))                  fChain->SetBranchAddress("tightMuon2",                 tightMuon2_,                 &b_tightMuon2                );
    if (fChain->GetBranch("softMuon2"))                   fChain->SetBranchAddress("softMuon2",                  softMuon2_,                  &b_softMuon2                 );
    if (fChain->GetBranch("hybridMuon2"))                 fChain->SetBranchAddress("hybridMuon2",                hybridMuon2_,                &b_hybridMuon2               );
    if (fChain->GetBranch("HPMuon2"))                     fChain->SetBranchAddress("HPMuon2",                    HPMuon2_,                    &b_HPMuon2                   );
    if (fChain->GetBranch("trigMuon2"))                   fChain->SetBranchAddress("trigMuon2",                 &trigMuon2_,                  &b_trigMuon2                 );
    if (fChain->GetBranch("nMatchedStationD2"))           fChain->SetBranchAddress("nMatchedStationD2",          nMatchedStationD2_,          &b_nMatchedStationD2         );
    if (fChain->GetBranch("nTrackerLayerD2"))             fChain->SetBranchAddress("nTrackerLayerD2",            nTrackerLayerD2_,            &b_nTrackerLayerD2           );
    if (fChain->GetBranch("nPixelLayerD2"))               fChain->SetBranchAddress("nPixelLayerD2",              nPixelLayerD2_,              &b_nPixelLayerD2             );
    if (fChain->GetBranch("nPixelHitD2"))                 fChain->SetBranchAddress("nPixelHitD2",                nPixelHitD2_,                &b_nPixelHitD2               );
    if (fChain->GetBranch("nMuonHitD2"))                  fChain->SetBranchAddress("nMuonHitD2",                 nMuonHitD2_,                 &b_nMuonHitD2                );
    if (fChain->GetBranch("GlbTrkChiD2"))                 fChain->SetBranchAddress("GlbTrkChiD2",                GlbTrkChiD2_,                &b_GlbTrkChiD2               );
    if (fChain->GetBranch("muondZD2"))                    fChain->SetBranchAddress("muondZD2",                   muondZD2_,                   &b_muondZD2                  );
    if (fChain->GetBranch("muondXYD2"))                   fChain->SetBranchAddress("muondXYD2",                  muondXYD2_,                  &b_muondXYD2                 );
    if (fChain->GetBranch("dZD2"))                        fChain->SetBranchAddress("dZD2",                       dZD2_,                       &b_dZD2                      );
    if (fChain->GetBranch("dXYD2"))                       fChain->SetBranchAddress("dXYD2",                      dXYD2_,                      &b_dXYD2                     );

    // SET GEN INFO BRANCHES
    if (fChain->GetBranch("weight_gen"))                  fChain->SetBranchAddress("weight_gen",                &weight_gen_,                 &b_weight_gen                );
    if (fChain->GetBranch("candSize_gen"))                fChain->SetBranchAddress("candSize_gen",              &candSize_gen_,               &b_candSize_gen              );
    if (fChain->GetBranch("pT_gen"))                      fChain->SetBranchAddress("pT_gen",                     pT_gen_,                     &b_pT_gen                    );
    if (fChain->GetBranch("eta_gen"))                     fChain->SetBranchAddress("eta_gen",                    eta_gen_,                    &b_eta_gen                   );
    if (fChain->GetBranch("y_gen"))                       fChain->SetBranchAddress("y_gen",                      y_gen_,                      &b_y_gen                     );
    if (fChain->GetBranch("status_gen"))                  fChain->SetBranchAddress("status_gen",                 status_gen_,                 &b_status_gen                );
    if (fChain->GetBranch("PID_gen"))                     fChain->SetBranchAddress("PID_gen",                    PID_gen_,                    &b_PID_gen                   );
    if (fChain->GetBranch("DecayID_gen"))                     fChain->SetBranchAddress("DecayID_gen",                    DecayID_gen_,                    &b_DecayID_gen                   );
    if (fChain->GetBranch("MotherID_gen"))                fChain->SetBranchAddress("MotherID_gen",               MotherID_gen_,               &b_MotherID_gen              );
    if (fChain->GetBranch("RecIdx_gen"))                  fChain->SetBranchAddress("RecIdx_gen",                 RecIdx_gen_,                 &b_RecIdx_gen                );
    if (fChain->GetBranch("PIDD1_gen"))                   fChain->SetBranchAddress("PIDD1_gen",                  PIDD1_gen_,                  &b_PIDD1_gen                 );
    if (fChain->GetBranch("chargeD1_gen"))                fChain->SetBranchAddress("chargeD1_gen",               chargeD1_gen_,               &b_chargeD1_gen              );
    if (fChain->GetBranch("pTD1_gen"))                    fChain->SetBranchAddress("pTD1_gen",                   pTD1_gen_,                   &b_pTD1_gen                  );
    if (fChain->GetBranch("EtaD1_gen"))                   fChain->SetBranchAddress("EtaD1_gen",                  EtaD1_gen_,                  &b_EtaD1_gen                 );
    if (fChain->GetBranch("PhiD1_gen"))                   fChain->SetBranchAddress("PhiD1_gen",                  PhiD1_gen_,                  &b_PhiD1_gen                 );
    if (fChain->GetBranch("PIDD2_gen"))                   fChain->SetBranchAddress("PIDD2_gen",                  PIDD2_gen_,                  &b_PIDD2_gen                 );
    if (fChain->GetBranch("chargeD2_gen"))                fChain->SetBranchAddress("chargeD2_gen",               chargeD2_gen_,               &b_chargeD2_gen              );
    if (fChain->GetBranch("pTD2_gen"))                    fChain->SetBranchAddress("pTD2_gen",                   pTD2_gen_,                   &b_pTD2_gen                  );
    if (fChain->GetBranch("EtaD2_gen"))                   fChain->SetBranchAddress("EtaD2_gen",                  EtaD2_gen_,                  &b_EtaD2_gen                 );
    if (fChain->GetBranch("PhiD2_gen"))                   fChain->SetBranchAddress("PhiD2_gen",                  PhiD2_gen_,                  &b_PhiD2_gen                 );
  }
};

void VertexCompositeNtuple::Clear(void)
{
  if (fChainM_.size()==0) return;

  // CLEAR EVENT INFO VARIABLES
  if (GetBranchStatus("RunNb")==1)        RunNb_        = 0;
  if (GetBranchStatus("LSNb")==1)         LSNb_         = 0;
  if (GetBranchStatus("EventNb")==1)      EventNb_      = 0;
  if (GetBranchStatus("nPV")==1)          nPV_          = -1;
  if (GetBranchStatus("bestvtxX")==1)     bestvtxX_     = -99.;
  if (GetBranchStatus("bestvtxY")==1)     bestvtxY_     = -99.;
  if (GetBranchStatus("bestvtxZ")==1)     bestvtxZ_     = -99.;
  if (GetBranchStatus("centrality")==1)   centrality_   = -1;
  if (GetBranchStatus("Npixel")==1)       Npixel_       = -1;
  if (GetBranchStatus("HFsumETPlus")==1)  HFsumETPlus_  = -1.;
  if (GetBranchStatus("HFsumETMinus")==1) HFsumETMinus_ = -1.;
  if (GetBranchStatus("ZDCPlus")==1)      ZDCPlus_      = -1.;
  if (GetBranchStatus("ZDCMinus")==1)     ZDCMinus_     = -1.;
  if (GetBranchStatus("Ntrkoffline")==1)  Ntrkoffline_  = -1;
  if (GetBranchStatus("trigPrescale")==1) std::fill(trigPrescale_, trigPrescale_+NTRG, -9);
  if (GetBranchStatus("trigHLT")==1)      std::fill(trigHLT_, trigHLT_+NTRG, 0);
  if (GetBranchStatus("evtSel")==1)       std::fill(evtSel_, evtSel_+NSEL, 0);

  // CLEAR EVENT PLANE VARIABLES
  if (GetBranchStatus("ephfpSumW")==1)  ephfpSumW_ = -99.;
  if (GetBranchStatus("ephfmSumW")==1)  ephfmSumW_ = -99.;
  if (GetBranchStatus("eptrackmidSumW")==1)  eptrackmidSumW_ = -99.;
  if (GetBranchStatus("ephfpAngle")==1) std::fill(ephfpAngle_, ephfpAngle_+NEP, -99.);
  if (GetBranchStatus("ephfpQ")==1)     std::fill(ephfpQ_, ephfpQ_+NEP, -99.);
  if (GetBranchStatus("ephfmAngle")==1) std::fill(ephfmAngle_, ephfmAngle_+NEP, -99.);
  if (GetBranchStatus("ephfmQ")==1)     std::fill(ephfmQ_, ephfmQ_+NEP, -99.);
  if (GetBranchStatus("eptrackmidAngle")==1) std::fill(eptrackmidAngle_, eptrackmidAngle_+NEP, -99.);
  if (GetBranchStatus("eptrackmidQ")==1)     std::fill(eptrackmidQ_, eptrackmidQ_+NEP, -99.);

  // CLEAR CANDIDATE INFO VARIABLES
  const auto& nCand = (candSize_>0 ? candSize_ : NCAND);
  if (GetBranchStatus("candSize")==1)                   candSize_ = 0;
  if (GetBranchStatus("pT")==1)                         std::fill(pT_, pT_+nCand, -1.);
  if (GetBranchStatus("eta")==1)                        std::fill(eta_, eta_+nCand, -9.);
  if (GetBranchStatus("y")==1)                          std::fill(y_, y_+nCand, -9.);
  if (GetBranchStatus("phi")==1)                        std::fill(phi_, phi_+nCand, -9.);
  if (GetBranchStatus("mass")==1)                       std::fill(mass_, mass_+nCand, -1.);
  if (GetBranchStatus("flavor")==1)                     std::fill(flavor_, flavor_+nCand, 0.);
  if (GetBranchStatus("VtxProb")==1)                    std::fill(VtxProb_, VtxProb_+nCand, -1.);
  if (GetBranchStatus("3DCosPointingAngle")==1)         std::fill(V3DCosPointingAngle_, V3DCosPointingAngle_+nCand, -9.);
  if (GetBranchStatus("3DPointingAngle")==1)            std::fill(V3DPointingAngle_, V3DPointingAngle_+nCand, -9.);
  if (GetBranchStatus("2DCosPointingAngle")==1)         std::fill(V2DCosPointingAngle_, V2DCosPointingAngle_+nCand, -9.);
  if (GetBranchStatus("2DPointingAngle")==1)            std::fill(V2DPointingAngle_, V2DPointingAngle_+nCand, -9.);
  if (GetBranchStatus("3DDecayLengthSignificance")==1)  std::fill(V3DDecayLengthSignificance_, V3DDecayLengthSignificance_+nCand, -1.);
  if (GetBranchStatus("3DDecayLength")==1)              std::fill(V3DDecayLength_, V3DDecayLength_+nCand, -1.);
  if (GetBranchStatus("2DDecayLengthSignificance")==1)  std::fill(V2DDecayLengthSignificance_, V2DDecayLengthSignificance_+nCand, -1.);
  if (GetBranchStatus("2DDecayLength")==1)              std::fill(V2DDecayLength_, V2DDecayLength_+nCand, -1.);
  if (GetBranchStatus("zDCASignificanceDaugther1")==1)  std::fill(zDCASignificanceDaugther1_, zDCASignificanceDaugther1_+nCand, -1.);
  if (GetBranchStatus("xyDCASignificanceDaugther1")==1) std::fill(xyDCASignificanceDaugther1_, xyDCASignificanceDaugther1_+nCand, -1.);
  if (GetBranchStatus("HighPurityDaugther1")==1)        std::fill(HighPurityDaugther1_, HighPurityDaugther1_+nCand, 0);
  if (GetBranchStatus("NHitD1")==1)                     std::fill(NHitD1_, NHitD1_+nCand, -1.);
  if (GetBranchStatus("pTD1")==1)                       std::fill(pTD1_, pTD1_+nCand, -1.);
  if (GetBranchStatus("EtaD1")==1)                      std::fill(EtaD1_, EtaD1_+nCand, -9.);
  if (GetBranchStatus("PhiD1")==1)                      std::fill(PhiD1_, PhiD1_+nCand, -9.);
  if (GetBranchStatus("chargeD1")==1)                   std::fill(chargeD1_, chargeD1_+nCand, -9);
  if (GetBranchStatus("dedxHarmonic2D1")==1)            std::fill(dedxHarmonic2D1_, dedxHarmonic2D1_+nCand, -1.);
  if (GetBranchStatus("zDCASignificanceDaugther2")==1)  std::fill(zDCASignificanceDaugther2_, zDCASignificanceDaugther2_+nCand, -1.);
  if (GetBranchStatus("xyDCASignificanceDaugther2")==1) std::fill(xyDCASignificanceDaugther2_, xyDCASignificanceDaugther2_+nCand, -1.);
  if (GetBranchStatus("HighPurityDaugther2")==1)        std::fill(HighPurityDaugther2_, HighPurityDaugther2_+nCand, 0);
  if (GetBranchStatus("NHitD2")==1)                     std::fill(NHitD2_, NHitD2_+nCand, -1.);
  if (GetBranchStatus("pTD2")==1)                       std::fill(pTD2_, pTD2_+nCand, -1.);
  if (GetBranchStatus("EtaD2")==1)                      std::fill(EtaD2_, EtaD2_+nCand, -9.);
  if (GetBranchStatus("PhiD2")==1)                      std::fill(PhiD2_, PhiD2_+nCand, -9.);
  if (GetBranchStatus("chargeD2")==1)                   std::fill(chargeD2_, chargeD2_+nCand, -9);
  if (GetBranchStatus("dedxHarmonic2D2")==1)            std::fill(dedxHarmonic2D2_, dedxHarmonic2D2_+nCand, -1.);
  if (GetBranchStatus("isSwap")==1)                     std::fill(isSwap_, isSwap_+nCand, 0);
  if (GetBranchStatus("idmom_reco")==1)                 std::fill(idmom_reco_, idmom_reco_+nCand, -999);
  if (GetBranchStatus("matchGEN")==1)                   std::fill(matchGEN_, matchGEN_+nCand, 0);
  if (GetBranchStatus("PIDD1")==1)                      std::fill(PIDD1_, PIDD1_+nCand, -999);
  if (GetBranchStatus("PIDD2")==1)                      std::fill(PIDD2_, PIDD2_+nCand, -999);

  // CLEAR MUON INFO VARIABLES
  if (GetBranchStatus("OneStMuon1")==1)         std::fill(OneStMuon1_, OneStMuon1_+nCand, 0);
  if (GetBranchStatus("PFMuon1")==1)            std::fill(PFMuon1_, PFMuon1_+nCand, 0);
  if (GetBranchStatus("GlbMuon1")==1)           std::fill(GlbMuon1_, GlbMuon1_+nCand, 0);
  if (GetBranchStatus("trkMuon1")==1)           std::fill(trkMuon1_, trkMuon1_+nCand, 0);
  if (GetBranchStatus("tightMuon1")==1)         std::fill(tightMuon1_, tightMuon1_+nCand, 0);
  if (GetBranchStatus("softMuon1")==1)          std::fill(softMuon1_, softMuon1_+nCand, 0);
  if (GetBranchStatus("hybridMuon1")==1)        std::fill(hybridMuon1_, hybridMuon1_+nCand, 0);
  if (GetBranchStatus("HPMuon1")==1)            std::fill(HPMuon1_, HPMuon1_+nCand, 0);
  if (GetBranchStatus("trigMuon1")==1 && trigMuon1_) trigMuon1_->clear();
  if (GetBranchStatus("nMatchedStationD1")==1)  std::fill(nMatchedStationD1_, nMatchedStationD1_+nCand, -1);
  if (GetBranchStatus("nTrackerLayerD1")==1)    std::fill(nTrackerLayerD1_, nTrackerLayerD1_+nCand, -1);
  if (GetBranchStatus("nPixelLayerD1")==1)      std::fill(nPixelLayerD1_, nPixelLayerD1_+nCand, -1);
  if (GetBranchStatus("nPixelHitD1")==1)        std::fill(nPixelHitD1_, nPixelHitD1_+nCand, -1);
  if (GetBranchStatus("nMuonHitD1")==1)         std::fill(nMuonHitD1_, nMuonHitD1_+nCand, -1);
  if (GetBranchStatus("GlbTrkChiD1")==1)        std::fill(GlbTrkChiD1_, GlbTrkChiD1_+nCand, 99.);
  if (GetBranchStatus("muondZD1")==1)           std::fill(muondZD1_, muondZD1_+nCand, -99.);
  if (GetBranchStatus("muondXYD1")==1)          std::fill(muondXYD1_, muondXYD1_+nCand, -99.);
  if (GetBranchStatus("dZD1")==1)               std::fill(dZD1_, dZD1_+nCand, -99.);
  if (GetBranchStatus("dXYD1")==1)              std::fill(dXYD1_, dXYD1_+nCand, -99.);
  if (GetBranchStatus("nMatchedChamberD1")==1)  std::fill(nMatchedChamberD1_, nMatchedChamberD1_+nCand, -1);
  if (GetBranchStatus("EnergyDepositionD1")==1) std::fill(EnergyDepositionD1_, EnergyDepositionD1_+nCand, -1.);
  if (GetBranchStatus("dx1_seg")==1)            std::fill(dx1_seg_, dx1_seg_+nCand, -99.);
  if (GetBranchStatus("dy1_seg")==1)            std::fill(dy1_seg_, dy1_seg_+nCand, -99.);
  if (GetBranchStatus("dxSig1_seg")==1)         std::fill(dxSig1_seg_, dxSig1_seg_+nCand, -99.);
  if (GetBranchStatus("dySig1_seg")==1)         std::fill(dySig1_seg_, dySig1_seg_+nCand, -99.);
  if (GetBranchStatus("ddxdz1_seg")==1)         std::fill(ddxdz1_seg_, ddxdz1_seg_+nCand, -99.);
  if (GetBranchStatus("ddydz1_seg")==1)         std::fill(ddydz1_seg_, ddydz1_seg_+nCand, -99.);
  if (GetBranchStatus("ddxdzSig1_seg")==1)      std::fill(ddxdzSig1_seg_, ddxdzSig1_seg_+nCand, -99.);
  if (GetBranchStatus("ddydzSig1_seg")==1)      std::fill(ddydzSig1_seg_, ddydzSig1_seg_+nCand, -99.);
  if (GetBranchStatus("OneStMuon2")==1)         std::fill(OneStMuon2_, OneStMuon2_+nCand, 0);
  if (GetBranchStatus("PFMuon2")==1)            std::fill(PFMuon2_, PFMuon2_+nCand, 0);
  if (GetBranchStatus("GlbMuon2")==1)           std::fill(GlbMuon2_, GlbMuon2_+nCand, 0);
  if (GetBranchStatus("trkMuon2")==1)           std::fill(trkMuon2_, trkMuon2_+nCand, 0);
  if (GetBranchStatus("tightMuon2")==1)         std::fill(tightMuon2_, tightMuon2_+nCand, 0);
  if (GetBranchStatus("softMuon2")==1)          std::fill(softMuon2_, softMuon2_+nCand, 0);
  if (GetBranchStatus("hybridMuon2")==1)        std::fill(hybridMuon2_, hybridMuon2_+nCand, 0);
  if (GetBranchStatus("HPMuon2")==1)            std::fill(HPMuon2_, HPMuon2_+nCand, 0);
  if (GetBranchStatus("trigMuon2")==1 && trigMuon2_) trigMuon2_->clear();
  if (GetBranchStatus("nMatchedStationD2")==1)  std::fill(nMatchedStationD2_, nMatchedStationD2_+nCand, -1);
  if (GetBranchStatus("nTrackerLayerD2")==1)    std::fill(nTrackerLayerD2_, nTrackerLayerD2_+nCand, -1);
  if (GetBranchStatus("nPixelLayerD2")==1)      std::fill(nPixelLayerD2_, nPixelLayerD2_+nCand, -1);
  if (GetBranchStatus("nPixelHitD2")==1)        std::fill(nPixelHitD2_, nPixelHitD2_+nCand, -1);
  if (GetBranchStatus("nMuonHitD2")==1)         std::fill(nMuonHitD2_, nMuonHitD2_+nCand, -1);
  if (GetBranchStatus("GlbTrkChiD2")==1)        std::fill(GlbTrkChiD2_, GlbTrkChiD2_+nCand, 99.);
  if (GetBranchStatus("muondZD2")==1)           std::fill(muondZD2_, muondZD2_+nCand, -99.);
  if (GetBranchStatus("muondXYD2")==1)          std::fill(muondXYD2_, muondXYD2_+nCand, -99.);
  if (GetBranchStatus("dZD2")==1)               std::fill(dZD2_, dZD2_+nCand, -99.);
  if (GetBranchStatus("dXYD2")==1)              std::fill(dXYD2_, dXYD2_+nCand, -99.);
  if (GetBranchStatus("nMatchedChamberD2")==1)  std::fill(nMatchedChamberD2_, nMatchedChamberD2_+nCand, -1);
  if (GetBranchStatus("EnergyDepositionD2")==1) std::fill(EnergyDepositionD2_, EnergyDepositionD2_+nCand, -1.);
  if (GetBranchStatus("dx2_seg")==1)            std::fill(dx2_seg_, dx2_seg_+nCand, -99.);
  if (GetBranchStatus("dy2_seg")==1)            std::fill(dy2_seg_, dy2_seg_+nCand, -99.);
  if (GetBranchStatus("dxSig2_seg")==1)         std::fill(dxSig2_seg_, dxSig2_seg_+nCand, -99.);
  if (GetBranchStatus("dySig2_seg")==1)         std::fill(dySig2_seg_, dySig2_seg_+nCand, -99.);
  if (GetBranchStatus("ddxdz2_seg")==1)         std::fill(ddxdz2_seg_, ddxdz2_seg_+nCand, -99.);
  if (GetBranchStatus("ddydz2_seg")==1)         std::fill(ddydz2_seg_, ddydz2_seg_+nCand, -99.);
  if (GetBranchStatus("ddxdzSig2_seg")==1)      std::fill(ddxdzSig2_seg_, ddxdzSig2_seg_+nCand, -99.);
  if (GetBranchStatus("ddydzSig2_seg")==1)      std::fill(ddydzSig2_seg_, ddydzSig2_seg_+nCand, -99.);

  // CLEAR GEN INFO VARIABLES
  const auto& nGen = (candSize_gen_>0 ? candSize_gen_ : NGEN);
  if (GetBranchStatus("weight_gen")==1)   weight_gen_ = -99.;
  if (GetBranchStatus("candSize_gen")==1) candSize_gen_ = 0;
  if (GetBranchStatus("pT_gen")==1)       std::fill(pT_gen_, pT_gen_+nGen, -1.);
  if (GetBranchStatus("eta_gen")==1)      std::fill(eta_gen_, eta_gen_+nGen, -9.);
  if (GetBranchStatus("y_gen")==1)        std::fill(y_gen_, y_gen_+nGen, -9.);
  if (GetBranchStatus("status_gen")==1)   std::fill(status_gen_, status_gen_+nGen, -9);
  if (GetBranchStatus("PID_gen")==1)      std::fill(PID_gen_, PID_gen_+nGen, -999);
  if (GetBranchStatus("DecayID_gen")==1)  std::fill(DecayID_gen_, DecayID_gen_+nGen, -999);
  if (GetBranchStatus("MotherID_gen")==1) std::fill(MotherID_gen_, MotherID_gen_+nGen, -999);
  if (GetBranchStatus("RecIdx_gen")==1)   std::fill(RecIdx_gen_, RecIdx_gen_+nGen, -1);
  if (GetBranchStatus("PIDD1_gen")==1)    std::fill(PIDD1_gen_, PIDD1_gen_+nGen, -999);
  if (GetBranchStatus("chargeD1_gen")==1) std::fill(chargeD1_gen_, chargeD1_gen_+nGen, -9);
  if (GetBranchStatus("pTD1_gen")==1)     std::fill(pTD1_gen_, pTD1_gen_+nGen, -1.);
  if (GetBranchStatus("EtaD1_gen")==1)    std::fill(EtaD1_gen_, EtaD1_gen_+nGen, -9.);
  if (GetBranchStatus("PhiD1_gen")==1)    std::fill(PhiD1_gen_, PhiD1_gen_+nGen, -9.);
  if (GetBranchStatus("PIDD2_gen")==1)    std::fill(PIDD2_gen_, PIDD2_gen_+nGen, -999);
  if (GetBranchStatus("chargeD2_gen")==1) std::fill(chargeD2_gen_, chargeD2_gen_+nGen, -9);
  if (GetBranchStatus("pTD2_gen")==1)     std::fill(pTD2_gen_, pTD2_gen_+nGen, -1.);
  if (GetBranchStatus("EtaD2_gen")==1)    std::fill(EtaD2_gen_, EtaD2_gen_+nGen, -9.);
  if (GetBranchStatus("PhiD2_gen")==1)    std::fill(PhiD2_gen_, PhiD2_gen_+nGen, -9.);
};

void VertexCompositeNtuple::GenerateDictionaries(void)
{
  std::vector< std::string > inList = {
    "vector<vector<UChar_t>>",
  };
  const std::string& CWD = getcwd(NULL, 0);
  gSystem->mkdir((CWD+"/cpp").c_str(), kTRUE);
  gSystem->ChangeDirectory((CWD+"/cpp").c_str());
  gInterpreter->AddIncludePath(Form("%s", (CWD+"/cpp").c_str())); // Needed to find the new dictionaries
  for (const auto& d : inList) { gInterpreter->GenerateDictionary(d.c_str(), "vector"); }
  gSystem->ChangeDirectory(CWD.c_str());
};

Bool_t VertexCompositeNtuple::tightMuon1(const UInt_t& iC, const std::string& type)
{
  if      (type==""   ) { return tightMuon1()[iC]; } 
  else if (type=="Y15") {
    return ( GlbMuon1()[iC] && (GlbTrkChiD1()[iC] < 10.) &&
             (nMuonHitD1()[iC] > 0) && (nMatchedStationD1()[iC] > 1) &&
             (nPixelHitD1()[iC] > 0) && (nTrackerLayerD1()[iC] > 5) &&
             (fabs(muondXYD1()[iC]) < 0.2) && (fabs(muondZD1()[iC]) < 0.5) );
  }
  else if (type=="POG") {
    return ( GlbMuon1()[iC] && PFMuon1()[iC] && (GlbTrkChiD1()[iC] < 10.) &&
             (nMuonHitD1()[iC] > 0) && (nMatchedStationD1()[iC] > 1) &&
             (nPixelHitD1()[iC] > 0) && (nTrackerLayerD1()[iC] > 5) &&
             (fabs(muondXYD1()[iC]) < 0.2) && (fabs(muondZD1()[iC]) < 0.5) );
  }
  else { std::cout << "[ERROR] Tight MuonID is not defined for " << type << std::endl; }
  return false;
};

Bool_t VertexCompositeNtuple::tightMuon2(const UInt_t& iC, const std::string& type)
{
  if      (type==""   ) { return tightMuon2()[iC]; }
  else if (type=="Y15") {
    return ( GlbMuon2()[iC] && (GlbTrkChiD2()[iC] < 10.) &&
             (nMuonHitD2()[iC] > 0) && (nMatchedStationD2()[iC] > 1) &&
             (nPixelHitD2()[iC] > 0) && (nTrackerLayerD2()[iC] > 5) &&
             (fabs(muondXYD2()[iC]) < 0.2) && (fabs(muondZD2()[iC]) < 0.5) );
  }
  else if (type=="POG") {
    return ( GlbMuon2()[iC] && PFMuon2()[iC] && (GlbTrkChiD2()[iC] < 10.) &&
             (nMuonHitD2()[iC] > 0) && (nMatchedStationD2()[iC] > 1) &&
             (nPixelHitD2()[iC] > 0) && (nTrackerLayerD2()[iC] > 5) &&
             (fabs(muondXYD2()[iC]) < 0.2) && (fabs(muondZD2()[iC]) < 0.5) );
  }
  else { std::cout << "[ERROR] Tight MuonID is not defined for " << type << std::endl; }
  return false;
};

Bool_t VertexCompositeNtuple::hybridMuon1(const UInt_t& iC, const std::string& type)
{
  if      (type==""   ) { return (hybridMuon1()[iC] && trkMuon1()[iC]); }
  else if (type=="Y15") {
    return ( GlbMuon1()[iC] && OneStMuon1()[iC] &&
             (nPixelLayerD1()[iC] > 0) && (nTrackerLayerD1()[iC] > 5) &&
             (fabs(dXYD1()[iC]) < 0.3) && (fabs(dZD1()[iC]) < 20.) );
  }
  else if (type=="Y18") {
    return ( GlbMuon1()[iC] && trkMuon1()[iC] &&
             (nPixelLayerD1()[iC] > 0) && (nTrackerLayerD1()[iC] > 5) &&
             (fabs(dXYD1()[iC]) < 0.3) && (fabs(dZD1()[iC]) < 20.) );
  }
  else { std::cout << "[ERROR] Hybrid MuonID is not defined for " << type << std::endl; }
  return false;
};

Bool_t VertexCompositeNtuple::hybridMuon2(const UInt_t& iC, const std::string& type)
{
  if      (type==""   ) { return (hybridMuon2()[iC] && trkMuon2()[iC]); }
  else if (type=="Y15") {
    return ( GlbMuon2()[iC] && OneStMuon2()[iC] &&
             (nPixelLayerD2()[iC] > 0) && (nTrackerLayerD2()[iC] > 5) &&
             (fabs(dXYD2()[iC]) < 0.3) && (fabs(dZD2()[iC]) < 20.) );
  }
  else if (type=="Y18") {
    return ( GlbMuon2()[iC] && trkMuon2()[iC] &&
             (nPixelLayerD2()[iC] > 0) && (nTrackerLayerD2()[iC] > 5) &&
             (fabs(dXYD2()[iC]) < 0.3) && (fabs(dZD2()[iC]) < 20.) );
  }
  else { std::cout << "[ERROR] Hybrid MuonID is not defined for " << type << std::endl; }
  return false;
};

Bool_t VertexCompositeNtuple::softMuon1(const UInt_t& iC, const std::string& type)
{
  if      (type==""   ) { return softMuon1()[iC]; }
  else if (type=="POG") {
    return ( OneStMuon1()[iC] && HPMuon1()[iC] &&
             (nPixelLayerD1()[iC] > 0) && (nTrackerLayerD1()[iC] > 5) &&
             (fabs(dXYD1()[iC]) < 0.3) && (fabs(dZD1()[iC]) < 20.) );
  }
  else { std::cout << "[ERROR] Soft MuonID is not defined for " << type << std::endl; }
  return false;
};

Bool_t VertexCompositeNtuple::softMuon2(const UInt_t& iC, const std::string& type)
{
  if      (type==""   ) { return softMuon2()[iC]; }
  else if (type=="POG") {
    return ( OneStMuon2()[iC] && HPMuon2()[iC] &&
             (nPixelLayerD2()[iC] > 0) && (nTrackerLayerD2()[iC] > 5) &&
             (fabs(dXYD2()[iC]) < 0.3) && (fabs(dZD2()[iC]) < 20.) );
  }
  else { std::cout << "[ERROR] Soft MuonID is not defined for " << type << std::endl; }
  return false;
};

Double_t VertexCompositeNtuple::phiAsym(const UInt_t& iC)
{
  const auto& pT1 = pTD1()[iC];
  const auto& phi1 = PhiD1()[iC];
  const auto& pT2 = pTD2()[iC];
  const auto& phi2 = PhiD2()[iC];
  TVector3 pTV1; pTV1.SetPtEtaPhi(pT1, 0.0, phi1);
  TVector3 pTV2; pTV2.SetPtEtaPhi(pT2, 0.0, phi2);
  const auto& pTVDif = 0.5*(pTV1 - pTV2);
  const auto& pTVSum = (pTV1 + pTV2);
  return pTVSum.Angle(pTVDif);
};

#endif
