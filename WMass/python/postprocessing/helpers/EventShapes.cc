#ifndef CMGTools_WMassTools_EventShapes_h
#define CMGTools_WMassTools_EventShapes_h

#include <iostream>
#include <vector>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>

// system include files
#include <memory>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "fastjet/contrib/Njettiness.hh"
#include <DataFormats/Math/interface/deltaR.h>
//#include <fastjet> 

using namespace fastjet;

class EventShapes {

public:

  typedef TTreeReaderValue<int>   rint;
  typedef TTreeReaderArray<float> rfloats;
  typedef TTreeReaderArray<int> rints;
  typedef TLorentzVector chgparticle;
  typedef std::vector<chgparticle> chgparticles;

  EventShapes() {};
  ~EventShapes() {};

  void setChgParticles(rint *nCp, rfloats *cpPt, rfloats *cpEta, rfloats *cpPhi, rfloats *cpMass, rints *cpPdgid){
    nCp_=nCp; Cp_pt_ = cpPt; Cp_eta_ = cpEta; Cp_phi_ = cpPhi; Cp_mass_ = cpMass;  Cp_pdgid_ = cpPdgid;
  }

  void setGoodLepton(rfloats *lepPt, rfloats *lepEta, rfloats *lepPhi){
    lep_pt_ = lepPt; lep_eta_ = lepEta; lep_phi_ = lepPhi;
  }

  double deltaRCustom(float leptEta, float leptPhi, float chgEta, float chgPhi){
    
    double deltar=-10;
    
    deltar = sqrt(pow(deltaPhi(leptPhi,chgPhi),2)+pow(leptEta-chgEta,2));

    return deltar;
  }
  
  void runNjettiness(){

    std::vector<fastjet::PseudoJet> lClusterParticles;

    for(int iP = 0, nP = **nCp_; iP < nP; ++iP) {

      if(((*Cp_pdgid_)[iP]!=22) && ((*Cp_pdgid_)[iP]!=130)){ //excluding gammas and neutral hadrons (from miniAOD pdgId)

	if(deltaRCustom((*lep_eta_)[0],(*lep_phi_)[0],(*Cp_eta_)[iP],(*Cp_phi_)[iP])>0.1){
	 
	  chgparticle cp = TLorentzVector();
	  cp.SetPtEtaPhiM((*Cp_pt_)[iP],(*Cp_eta_)[iP],(*Cp_phi_)[iP],(*Cp_mass_)[iP]);
	  
	  fastjet::PseudoJet   pPart(cp.Px(),cp.Py(),cp.Pz(),cp.E());
	  lClusterParticles.push_back(pPart);
	}

      }
    }

    //njettiness variables
    fastjet::contrib::NormalizedMeasure normalizedMeasure(1.0,0.4);
    fastjet::contrib::Njettiness routine(fastjet::contrib::Njettiness::onepass_kt_axes,normalizedMeasure);
    float iTau1 = routine.getTau(1.,lClusterParticles);
    float iTau2 = routine.getTau(2.,lClusterParticles);
    float iTau3 = routine.getTau(3.,lClusterParticles);
    float iTau4 = routine.getTau(4.,lClusterParticles);  
    
    tau1_ = iTau1;
    tau2_ = iTau2;
    tau3_ = iTau3;
    tau4_ = iTau4;
  }

  chgparticles buildchgparticles(){
    
    chgparticles outputpart;

    for(int iP = 0, nP = **nCp_; iP < nP; ++iP) {

      if(((*Cp_pdgid_)[iP]!=22) && ((*Cp_pdgid_)[iP]!=130)){ //excluding gammas and neutral hadrons (from miniAOD pdgId)

	if(deltaRCustom((*lep_eta_)[0],(*lep_phi_)[0],(*Cp_eta_)[iP],(*Cp_phi_)[iP])>0.1){
	  
	  chgparticle cp = TLorentzVector();
	  cp.SetPtEtaPhiM((*Cp_pt_)[iP],(*Cp_eta_)[iP],(*Cp_phi_)[iP],(*Cp_mass_)[iP]);
	  
	  outputpart.push_back(cp);
	}
      }
    }

    return outputpart;
  }    

  const float & tau1() { return tau1_; }
  const float & tau2() { return tau2_; }
  const float & tau3() { return tau3_; }
  const float & tau4() { return tau4_; }
 

private:

  rint *nCp_ = nullptr;
  rfloats *Cp_pt_ = nullptr;
  rfloats *Cp_eta_ = nullptr;
  rfloats *Cp_phi_ = nullptr;
  rfloats *Cp_mass_ = nullptr;
  rints *Cp_pdgid_ = nullptr;

  rfloats *lep_pt_ = nullptr;
  rfloats *lep_phi_ = nullptr;
  rfloats *lep_eta_ = nullptr;
  
  float tau1_;   float tau2_;   float tau3_;   float tau4_; 

};

#endif
