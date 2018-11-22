#include "fwdet_res.h"
#include "particleidentyf.h"

#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hfwdetcand.h"
#include "hfwdetcandsim.h"

#include "hparticlecandsim.h"
#include "hparticlecand.h"
#include "hparticletool.h"

#include "hloop.h"
#include "hcategory.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <sstream>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <cassert>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";


using namespace std;

Int_t  nentries;
TCutG* leptonCut;
TFile* file;

double trackDistance(HParticleCand* track1, HParticleCand*  track2)
{
  double dist;
  HGeomVector base_1, base_2, dir_1, dir_2;
  HParticleTool p_tool;

  p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);
  p_tool.calcSegVector(track2->getZ(),track2->getR(),TMath::DegToRad()*track2->getPhi(),TMath::DegToRad()*track2->getTheta(),base_2,dir_2);
  dist=p_tool.calculateMinimumDistance(base_1,dir_1,base_2,dir_2);
  return dist;
}

double trackDistance(HParticleCand* track1, HFwDetCand*  track2)
{
  double dist;
  HGeomVector base_1, base_2, dir_1, dir_2;
  HParticleTool p_tool;

  p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);

  base_2.setX(track2->getPointX());
  base_2.setY(track2->getPointY());
  base_2.setZ(track2->getPointZ());
  dir_2.setX(track2->getDirTx());
  dir_2.setY(track2->getDirTy());
  dir_2.setZ(1);//konwencja, tak jest ustawione w fwdetstrawvec
  
  dist=p_tool.calculateMinimumDistance(base_1,dir_1,base_2,dir_2);
  return dist;
}

HGeomVector trackVertex(HParticleCand* track1, HParticleCand*  track2)
{
  HGeomVector ver;
  HGeomVector base_1, base_2, dir_1, dir_2;
  HParticleTool p_tool;

  p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);
  p_tool.calcSegVector(track2->getZ(),track2->getR(),TMath::DegToRad()*track2->getPhi(),TMath::DegToRad()*track2->getTheta(),base_2,dir_2);
  ver=p_tool.calcVertexAnalytical(base_1,dir_1,base_2,dir_2);
  return ver;
}

HGeomVector trackVertex(HParticleCand* track1, HFwDetCand*  track2)
{
  HGeomVector ver;
  HGeomVector base_1, base_2, dir_1, dir_2;
  HParticleTool p_tool;

  p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);

  base_2.setX(track2->getPointX());
  base_2.setY(track2->getPointY());
  base_2.setZ(track2->getPointZ());
  dir_2.setX(track2->getDirTx());
  dir_2.setY(track2->getDirTy());
  dir_2.setZ(1);//konwencja, tak jest ustawione w fwdetstrawvec
  
  ver=p_tool.calcVertexAnalytical(base_1,dir_1,base_2,dir_2);
  return ver;
}

Int_t getMotherIndex(HGeantKine* particle)
{
  Int_t trackID=particle->getTrack();
  HGeantKine* particleParent=particle->getParent(trackID);
  Int_t parentID=0;
  if(particleParent!=0)
    parentID=particleParent->getID();
  return parentID;
      
  //return particle->getGeneratorInfo1();
}

bool isLepton(HParticleCandSim* particle)
{
  double delta=0.1;
  double mquality=particle->getRichMatchingQuality();

  double dphi=particle->getDeltaPhi();
  double dtheta=particle->getDeltaTheta();

  if(mquality==-1 /*|| mquality>5*/)
    return false;
  //if(particle->getBeta()<(1-delta) || particle->getBeta()>(1+delta))
  //return false;
  if(!(leptonCut->IsInside(particle->getMomentum(),particle->getBeta())))
    return false;
  //if(dtheta<-0.4 || dtheta>0.4)
  //return false;
  //if(dphi*TMath::Sin(particle->getTheta())>0.4 || dphi*TMath::Sin(particle->getTheta())<-0.4)
  //return false;
  return true;
}

bool isLepton(HParticleCandSim* particle, TH2F* h1, TH2F* hphi1, TH2F* htheta1, TH2F* h2, TH2F* hphi2, TH2F* htheta2)
{
  double delta=0.1;
  double mquality=particle->getRichMatchingQuality();

  double dphi=particle->getDeltaPhi();
  double dtheta=particle->getDeltaTheta();

  double mom=particle->getMomentum();

  if(!particle->isFlagBit(kIsUsed))
    return false;
  
  if(particle->getGeantPID()==2 || particle->getGeantPID()==3)
    {
      h1->Fill(mom,particle->getBeta());
      htheta1->Fill(dtheta,mom);
      hphi1->Fill(dphi,mom);
    }
  if(isLepton(particle))
    {
      h2->Fill(mom,particle->getBeta());
      htheta2->Fill(dtheta,mom);
      hphi2->Fill(dphi,mom);
      return true;
    }
  else
    return false;
}

bool isProton(HParticleCandSim* particle, TH2F* h1, TH2F* h2)
{
  double mom=particle->getMomentum();
  //double tof=particle->getTof();
  double beta=particle->getBeta();
  double mass=TMath::Sqrt(mom*mom * (  1. / (beta*beta)  - 1. ));

  if(!particle->isFlagBit(kIsUsed))
    return false;

  if(particle->getGeantPID()==14)
    h1->Fill(particle->getMomentum(),particle->getBeta());
    
  if(mass < 1127 && mass > 650 && !isLepton(particle))
    {
      h2->Fill(particle->getMomentum(),particle->getBeta());
      return true;
    }
  
  return false;
}

bool isPion(HParticleCandSim* particle, TH2F* h1, TH2F* h2)
{
  double mom=particle->getMomentum();
  //double tof=particle->getTof();
  double beta=particle->getBeta();
  double mass=TMath::Sqrt(mom*mom * (  1. / (beta*beta)  - 1. ));

  if(!particle->isFlagBit(kIsUsed))
    return false;
  
  if(particle->getGeantPID()==8 || particle->getGeantPID()==9)
    h1->Fill(particle->getMomentum(),particle->getBeta());
  
  if(isLepton(particle))
    return false;
  
  if(mass < 240 && mass > 40)
    {
      h2->Fill(particle->getMomentum(),particle->getBeta());
      return true;
    }
  return false;
}
//#include "CreateHistos.cxx"
//#include "WriteHistos.cxx"



//--------------------------------------------------------------------

Int_t fwdet_tests(HLoop * loop, const AnaParameters & anapars)
{
  if (!loop->setInput(""))
    {                                                    // reading file structure
      std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  TStopwatch timer;
  timer.Reset();
  timer.Start();
	    


  loop->printCategories(); 

  file= new TFile("lepton_graph_cut.root");
  assert(file);
  leptonCut=(TCutG*)file->Get("CUTG");
  assert(leptonCut);
  file->Close();

  HCategory * fCatGeantKine = nullptr;
  fCatGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");
  if (!fCatGeantKine)
    {
      cout << "No catGeantKine!" << endl;
      exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    
  HCategory * fFwDetStrawCal = nullptr;
  fFwDetStrawCal = HCategoryManager::getCategory(catFwDetStrawCal, kTRUE, "catFwDetStrawCalSim");
  if (!fFwDetStrawCal)
    {
      cout << "No catFwDetStrawCal!" << endl;
      exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

  HCategory * fwDetCatSim = nullptr;
  fwDetCatSim = HCategoryManager::getCategory(catFwDetCand, kTRUE, "catFwDetCand");
  if (!fwDetCatSim)
    {
      cout << "No catFwDetCand!" << endl;
      //exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }
    
  HCategory * particleCatSim= nullptr;
  particleCatSim = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCand");
  if(!particleCatSim)
    {
      cout<< "No catParticleCandSim!"<<endl;
    }

  //************************** HISTOS ********************
    
  TH1F* hinvM_pmHpHAll=new TH1F("hinvM_pmHpHAll","hinvM_pmHpHAll",1000,1000,2200);
  TH1F* hinvM_pmHpHDist=new TH1F("hinvM_pmHpHDist","hinvM_pmHpHDist",1000,1000,2200);
  TH1F* hinvM_pmHpHDistZ=new TH1F("hinvM_pmHpHDistZ","hinvM_pmHpHDistZ",1000,1000,2200);
  TH1F* hinvM_pmHpHDistL=new TH1F("hinvM_pmHpHDistL","hinvM_pmHpHDistL",1000,1000,2200);
  TH1F* hinvM_pmHpHDistZL=new TH1F("hinvM_pmHpHDistZL","hinvM_pmHpHDistZL",1000,1000,2200);

  TH1F *hinvM_pmHpFTAll = new TH1F("hinvM_pmHpFTAll","hinvM_pmHpFTAll",1000,1000,2200);
  TH1F* hinvM_pmHpFTDist=new TH1F("hinvM_pmHpFTDist","hinvM_pmHpFTDist",1000,1000,2200);
  TH1F* hinvM_pmHpFTDistZ=new TH1F("hinvM_pmHpFTDistZ","hinvM_pmHpFTDistZ",1000,1000,2200);
  TH1F* hinvM_pmHpFTDistL=new TH1F("hinvM_pmHpFTDistL","hinvM_pmHpFTDistL",1000,1000,2200);
  TH1F* hinvM_pmHpFTDistZL=new TH1F("hinvM_pmHpFTDistZL","hinvM_pmHpFTDistZL",1000,1000,2200);

    
  TH1F* hLRmass=new TH1F("hLRmass","hLRmass",1000,1000,2200);
  TH1F* hLRmassDist=new TH1F("hLRmassDist","hLRmassDist",1000,1000,2200);
  TH1F* hLRmassDistZ=new TH1F("hLRmassDistZ","hLRmassDistZ",1000,1000,2200);
  TH1F* hLRmassDistL=new TH1F("hLRmassDistL","hLRmassDistL",1000,1000,2200);
  TH1F* hLRmassDistZL=new TH1F("hLRmassDistZL","hLRmassDistZL",1000,1000,2200);

  TH1F* hLRmassFT=new TH1F("hLRmassFT","hLRmassFT",1000,1000,2200);
  TH1F* hLRmassFTDist=new TH1F("hLRmassFTDist","hLRmassFTDist",1000,1000,2200);
  TH1F* hLRmassFTDistZ=new TH1F("hLRmassFTDistZ","hLRmassFTDistZ",1000,1000,2200);
  TH1F* hLRmassFTDistL=new TH1F("hLRmassFTDistL","hLRmassFTDistL",1000,1000,2200);
  TH1F* hLRmassFTDistZL=new TH1F("hLRmassFTDistZL","hLRmassFTDistZL",1000,1000,2200);

    
    
  TH1F* hdist_pmHpHAll=new TH1F("hdist_pmHpHAll","hdist_pmHpHAll",1000,0,300);
  TH1F* hLRdist=new TH1F("hLRdist","hLRdist",1000,0,300);
    
  TH1F* hdist_pmHpFTAll=new TH1F("hdist_pmHpFTAll","hdist_pmHpFTAll",1000,0,300);
  TH1F* hLRdistFT=new TH1F("hLRdistFT","hLRdistFT",1000,0,300);

  TH1F* hDLdistance=new TH1F("hDLdistance","hDLdistance",1000,0,800);
  TH1F* hDLdistanceFT=new TH1F("hDLdistanceFT","hDLdistanceFT",1000,0,800);
  
  TH1F *h2L1520vertex = new TH1F("h2L1520vertex","h2L1520vertex",1000,-300,300);
  TH1F *h2LHvertex = new TH1F("h2LHvertex","h2LHvertex",1000,-1000,1000);

  //**************************************************************************** 
   
  TH1F *hDLmassAll=new TH1F("hDLmassAll","hDLmassAll",1000,0,700);
  TH1F *hL1520massAll=new TH1F("hL1520massAll","hL1520massAll",1000,1200,2200);

  TH1F *hDLmassAllRL=new TH1F("hDLmassAllRL","hDLmassAllRL",1000,0,700);
  TH1F *hL1520massAllRL=new TH1F("hL1520massAllRL","hL1520massAllRL",1000,1200,2200);

  TH1F *hDLmassAllRL_L=new TH1F("hDLmassAllRL_L","hDLmassAllRL_L",1000,0,700);
  TH1F *hL1520massAllRL_L=new TH1F("hL1520massAllRL_L","hL1520massAllRL_L",1000,1200,2200);

     
  TH1F *hDLmassDist=new TH1F("hDLmassDist","hDLmassDist",1000,0,700);
  TH1F *hL1520massDist=new TH1F("hL1520massDist","hL1520massDist",1000,1200,2200);

  TH1F *hDLmassDistRL=new TH1F("hDLmassDistRL","hDLmassDistRL",1000,0,700);
  TH1F *hL1520massDistRL=new TH1F("hL1520massDistRL","hL1520massDistRL",1000,1200,2200);

  TH1F *hDLmassDistRL_L=new TH1F("hDLmassDistRL_L","hDLmassDistRL_L",1000,0,700);
  TH1F *hL1520massDistRL_L=new TH1F("hL1520massDistRL_L","hL1520massDistRL_L",1000,1200,2200);
     
  TH1F *hDLmassDistZ=new TH1F("hDLmassDistZ","hDLmassDistZ",1000,0,700);
  TH1F *hL1520massDistZ=new TH1F("hL1520massDistZ","hL1520massDistZ",1000,1200,2200);

  TH1F *hDLmassDistZRL=new TH1F("hDLmassDistZRL","hDLmassDistZRL",1000,0,700);
  TH1F *hL1520massDistZRL=new TH1F("hL1520massDistZRL","hL1520massDistZRL",1000,1200,2200);

  TH1F *hDLmassDistZRL_L=new TH1F("hDLmassDistZRL_L","hDLmassDistZRL_L",1000,0,700);
  TH1F *hL1520massDistZRL_L=new TH1F("hL1520massDistZRL_L","hL1520massDistZRL_L",1000,1200,2200);

     
  TH1F *hDLmassDistL=new TH1F("hDLmassDistL","hDLmassDistL",1000,0,700);
  TH1F *hL1520massDistL=new TH1F("hL1520massDistL","hL1520massDistL",1000,1200,2200);
  TH1F *hL1520massDistLpi0=new TH1F("hL1520massDistLpi0","hL1520massDistLpi0",1000,1200,2200);
     
  TH1F *hDLmassDistLRL=new TH1F("hDLmassDistLRL","hDLmassDistLRL",1000,0,700);
  TH1F *hL1520massDistLRL=new TH1F("hL1520massDistLRL","hL1520massDistLRL",1000,1200,2200);

  TH1F *hDLmassDistLRL_L=new TH1F("hDLmassDistLRL_L","hDLmassDistLRL_L",1000,0,700);
  TH1F *hL1520massDistLRL_L=new TH1F("hL1520massDistLRL_L","hL1520massDistLRL_L",1000,1200,2200);
     
  TH1F *hDLmassDistZL=new TH1F("hDLmassDistZL","hDLmassDistZL",1000,0,700);
  TH1F *hL1520massDistZL=new TH1F("hL1520massDistZL","hL1520massDistZL",1000,1200,2200);

  TH1F *hDLmassDistZLRL=new TH1F("hDLmassDistZLRL","hDLmassDistZLRL",1000,0,700);
  TH1F *hL1520massDistZLRL=new TH1F("hL1520massDistZLRL","hL1520massDistZLRL",1000,1200,2200);
  TH1F *hDLmassDistZLRL_L=new TH1F("hDLmassDistZLRL_L","hDLmassDistZLRL_L",1000,0,700);
  TH1F *hL1520massDistZLRL_L=new TH1F("hL1520massDistZLRL_L","hL1520massDistZLRL_L",1000,1200,2200);


     
  TH1F *hDLmassFinal=new TH1F("hDLmassFinal","hDLmassFinal",1000,0,700);
  TH1F *hL1520massFinal=new TH1F("hL1520massFinal","hL1520massFinal",1000,1200,2200);

  TH1F *hDLmassFinalRL=new TH1F("hDLmassFinalRL","hDLmassFinalRL",1000,0,700);
  TH1F *hL1520massFinalRL=new TH1F("hL1520massFinalRL","hL1520massFinalRL",1000,1200,2200);

  TH1F *hDLmassFinalRL_L=new TH1F("hDLmassFinalRL_L","hDLmassFinalRL_L",1000,0,700);
  TH1F *hL1520massFinalRL_L=new TH1F("hL1520massFinalRL_L","hL1520massFinalRL_L",1000,1200,2200);
     
     
  TH1F *hL1520massFinalRLpi0=new TH1F("hL1520massFinalRLpi0","hL1520massFinalRLpi0",1000,1200,2200);
  TH1F *hL1520massFinalRLpi0_L=new TH1F("hL1520massFinalRLpi0_L","hL1520massFinalRLpi0_L",1000,1200,2200);


  TH1F *hL1520massFinalpi0=new TH1F("hL1520massFinalpi0","hL1520massFinalpi0",1000,1200,2200);
  TH1F *hL1520massDistZLpi0=new TH1F("hL1520massDistZLpi0","hL1520massDistZLpi0",1000,1200,2200);
     

     
  //**************************************************************************** 
   
  TH1F *hDLmassFTAll=new TH1F("hDLmassFTAll","hDLmassFTAll",1000,0,700);
  TH1F *hL1520massFTAll=new TH1F("hL1520massFTAll","hL1520massFTAll",1000,1200,2200);

  TH1F *hDLmassFTAllRL=new TH1F("hDLmassFTAllRL","hDLmassFTAllRL",1000,0,700);
  TH1F *hL1520massFTAllRL=new TH1F("hL1520massFTAllRL","hL1520massFTAllRL",1000,1200,2200);
     
  TH1F *hDLmassFTDist=new TH1F("hDLmassFTDist","hDLmassFTDist",1000,0,700);
  TH1F *hL1520massFTDist=new TH1F("hL1520massFTDist","hL1520massFTDist",1000,1200,2200);

  TH1F *hDLmassFTDistRL=new TH1F("hDLmassFTDistRL","hDLmassFTDistRL",1000,0,700);
  TH1F *hL1520massFTDistRL=new TH1F("hL1520massFTDistRL","hL1520massFTDistRL",1000,1200,2200);
     
  TH1F *hDLmassFTDistZ=new TH1F("hDLmassFTDistZ","hDLmassFTDistZ",1000,0,700);
  TH1F *hL1520massFTDistZ=new TH1F("hL1520massFTDistZ","hL1520massFTDistZ",1000,1200,2200);

  TH1F *hDLmassFTDistZRL=new TH1F("hDLmassFTDistZRL","hDLmassFTDistZRL",1000,0,700);
  TH1F *hL1520massFTDistZRL=new TH1F("hL1520massFTDistZRL","hL1520massFTDistZRL",1000,1200,2200);
     
  TH1F *hDLmassFTDistL=new TH1F("hDLmassFTDistL","hDLmassFTDistL",1000,0,700);
  TH1F *hL1520massFTDistL=new TH1F("hL1520massFTDistL","hL1520massFTDistL",1000,1200,2200);
  TH1F *hL1520massFTDistLpi0=new TH1F("hL1520massFTDistLpi0","hL1520massFTDistLpi0",1000,1200,2200);

  TH1F *hDLmassFTDistLRL=new TH1F("hDLmassFTDistLRL","hDLmassFTDistLRL",1000,0,700);
  TH1F *hL1520massFTDistLRL=new TH1F("hL1520massFTDistLRL","hL1520massFTDistLRL",1000,1200,2200);
     
  TH1F *hDLmassFTDistZL=new TH1F("hDLmassFTDistZL","hDLmassFTDistZL",1000,0,700);
  TH1F *hL1520massFTDistZL=new TH1F("hL1520massFTDistZL","hL1520massFTDistZL",1000,1200,2200);

  TH1F *hDLmassFTDistZLRL=new TH1F("hDLmassFTDistZLRL","hDLmassFTDistZLRL",1000,0,700);
  TH1F *hL1520massFTDistZLRL=new TH1F("hL1520massFTDistZLRL","hL1520massFTDistZLRL",1000,1200,2200);
  TH1F *hL1520massFTDistZLRLpi0=new TH1F("hL1520massFTDistZLRLpi0","hL1520massFTDistZLRLpi0",1000,1200,2200);

  TH1F *hL1520massFTDistZLpi0=new TH1F("hL1520massFTDistZLpi0","hL1520massFTDistZLpi0",1000,1200,2200);
    
  TH1F *hDLmassFTFinal=new TH1F("hDLmassFTFinal","hDLmassFTFinal",1000,0,700);
  TH1F *hL1520massFTFinal=new TH1F("hL1520massFTFinal","hL1520massFTFinal",1000,1200,2200);

  TH1F *hDLmassFTFinalRL=new TH1F("hDLmassFTFinalRL","hDLmassFTFinalRL",1000,0,700);
  TH1F *hL1520massFTFinalRL=new TH1F("hL1520massFTFinalRL","hL1520massFTFinalRL",1000,1200,2200);

  TH1F *hL1520massFTFinalRLpi0=new TH1F("hL1520massFTFinalRLpi0","hL1520massFTFinalRLpi0",1000,1200,2200);
  TH1F *hL1520massFTFinalpi0=new TH1F("hL1520massFTFinalpi0","hL1520massFTFinalpi0",1000,1200,2200);

  TH1F *hL1520massDistZLRLpi0_L=new TH1F("hL1520massDistZLRLpi0_L","hL1520massDistZLRLpi0_L",1000,1200,2200);
  TH1F *hL1520massDistZLRLpi0=new TH1F("hL1520massDistZLRLpi0","hL1520massDistZLRLpi0",1000,1200,2200);

  TH1F *hL1520massDistZLpi0_L=new TH1F("hL1520massDistZLpi0_L","hL1520massDistZLpi0_L",1000,1200,2200);

 

  
     
  //*****************************     
  //**************************************************************************** 
   
  TH1F *hDLmassAll_epep=new TH1F("hDLmassAll_epep","hDLmassAll_epep",1000,0,700);
  TH1F *hL1520massAll_epep=new TH1F("hL1520massAll_epep","hL1520massAll_epep",1000,1200,2200);

  TH1F *hDLmassAllRL_epep=new TH1F("hDLmassAllRL_epep","hDLmassAllRL_epep",1000,0,700);
  TH1F *hL1520massAllRL_epep=new TH1F("hL1520massAllRL_epep","hL1520massAllRL_epep",1000,1200,2200);

  TH1F *hDLmassAllRL_L_epep=new TH1F("hDLmassAllRL_L_epep","hDLmassAllRL_L_epep",1000,0,700);
  TH1F *hL1520massAllRL_L_epep=new TH1F("hL1520massAllRL_L_epep","hL1520massAllRL_L_epep",1000,1200,2200);

     
  TH1F *hDLmassDist_epep=new TH1F("hDLmassDist_epep","hDLmassDist_epep",1000,0,700);
  TH1F *hL1520massDist_epep=new TH1F("hL1520massDist_epep","hL1520massDist_epep",1000,1200,2200);

  TH1F *hDLmassDistRL_epep=new TH1F("hDLmassDistRL_epep","hDLmassDistRL_epep",1000,0,700);
  TH1F *hL1520massDistRL_epep=new TH1F("hL1520massDistRL_epep","hL1520massDistRL_epep",1000,1200,2200);

  TH1F *hDLmassDistRL_L_epep=new TH1F("hDLmassDistRL_L_epep","hDLmassDistRL_L_epep",1000,0,700);
  TH1F *hL1520massDistRL_L_epep=new TH1F("hL1520massDistRL_L_epep","hL1520massDistRL_L_epep",1000,1200,2200);
     
  TH1F *hDLmassDistZ_epep=new TH1F("hDLmassDistZ_epep","hDLmassDistZ_epep",1000,0,700);
  TH1F *hL1520massDistZ_epep=new TH1F("hL1520massDistZ_epep","hL1520massDistZ_epep",1000,1200,2200);

  TH1F *hDLmassDistZRL_epep=new TH1F("hDLmassDistZRL_epep","hDLmassDistZRL_epep",1000,0,700);
  TH1F *hL1520massDistZRL_epep=new TH1F("hL1520massDistZRL_epep","hL1520massDistZRL_epep",1000,1200,2200);

  TH1F *hDLmassDistZRL_L_epep=new TH1F("hDLmassDistZRL_L_epep","hDLmassDistZRL_L_epep",1000,0,700);
  TH1F *hL1520massDistZRL_L_epep=new TH1F("hL1520massDistZRL_L_epep","hL1520massDistZRL_L_epep",1000,1200,2200);

     
  TH1F *hDLmassDistL_epep=new TH1F("hDLmassDistL_epep","hDLmassDistL_epep",1000,0,700);
  TH1F *hL1520massDistL_epep=new TH1F("hL1520massDistL_epep","hL1520massDistL_epep",1000,1200,2200);
  TH1F *hL1520massDistLpi0_epep=new TH1F("hL1520massDistLpi0_epep","hL1520massDistLpi0_epep",1000,1200,2200);
     
  TH1F *hDLmassDistLRL_epep=new TH1F("hDLmassDistLRL_epep","hDLmassDistLRL_epep",1000,0,700);
  TH1F *hL1520massDistLRL_epep=new TH1F("hL1520massDistLRL_epep","hL1520massDistLRL_epep",1000,1200,2200);

  TH1F *hDLmassDistLRL_L_epep=new TH1F("hDLmassDistLRL_L_epep","hDLmassDistLRL_L_epep",1000,0,700);
  TH1F *hL1520massDistLRL_L_epep=new TH1F("hL1520massDistLRL_L_epep","hL1520massDistLRL_L_epep",1000,1200,2200);
     
  TH1F *hDLmassDistZL_epep=new TH1F("hDLmassDistZL_epep","hDLmassDistZL_epep",1000,0,700);
  TH1F *hL1520massDistZL_epep=new TH1F("hL1520massDistZL_epep","hL1520massDistZL_epep",1000,1200,2200);

  TH1F *hDLmassDistZLRL_epep=new TH1F("hDLmassDistZLRL_epep","hDLmassDistZLRL_epep",1000,0,700);
  TH1F *hL1520massDistZLRL_epep=new TH1F("hL1520massDistZLRL_epep","hL1520massDistZLRL_epep",1000,1200,2200);
  TH1F *hDLmassDistZLRL_L_epep=new TH1F("hDLmassDistZLRL_L_epep","hDLmassDistZLRL_L_epep",1000,0,700);
  TH1F *hL1520massDistZLRL_L_epep=new TH1F("hL1520massDistZLRL_L_epep","hL1520massDistZLRL_L_epep",1000,1200,2200);


     
  TH1F *hDLmassFinal_epep=new TH1F("hDLmassFinal_epep","hDLmassFinal_epep",1000,0,700);
  TH1F *hL1520massFinal_epep=new TH1F("hL1520massFinal_epep","hL1520massFinal_epep",1000,1200,2200);

  TH1F *hDLmassFinalRL_epep=new TH1F("hDLmassFinalRL_epep","hDLmassFinalRL_epep",1000,0,700);
  TH1F *hL1520massFinalRL_epep=new TH1F("hL1520massFinalRL_epep","hL1520massFinalRL_epep",1000,1200,2200);

  TH1F *hDLmassFinalRL_L_epep=new TH1F("hDLmassFinalRL_L_epep","hDLmassFinalRL_L_epep",1000,0,700);
  TH1F *hL1520massFinalRL_L_epep=new TH1F("hL1520massFinalRL_L_epep","hL1520massFinalRL_L_epep",1000,1200,2200);
     
     
  TH1F *hL1520massFinalRLpi0_epep=new TH1F("hL1520massFinalRLpi0_epep","hL1520massFinalRLpi0_epep",1000,1200,2200);

  TH1F *hL1520massFinalRLpi0_L_epep=new TH1F("hL1520massFinalRLpi0_L_epep","hL1520massFinalRLpi0_L_epep",1000,1200,2200);


  TH1F *hL1520massFinalpi0_epep=new TH1F("hL1520massFinalpi0_epep","hL1520massFinalpi0_epep",1000,1200,2200);
  TH1F *hL1520massDistZLpi0_epep=new TH1F("hL1520massDistZLpi0_epep","hL1520massDistZLpi0_epep",1000,1200,2200);
     

     
  //**************************************************************************** 

  //**************************************************************************** 
   
  TH1F *hDLmassAll_emem=new TH1F("hDLmassAll_emem","hDLmassAll_emem",1000,0,700);
  TH1F *hL1520massAll_emem=new TH1F("hL1520massAll_emem","hL1520massAll_emem",1000,1200,2200);

  TH1F *hDLmassAllRL_emem=new TH1F("hDLmassAllRL_emem","hDLmassAllRL_emem",1000,0,700);
  TH1F *hL1520massAllRL_emem=new TH1F("hL1520massAllRL_emem","hL1520massAllRL_emem",1000,1200,2200);

  TH1F *hDLmassAllRL_L_emem=new TH1F("hDLmassAllRL_L_emem","hDLmassAllRL_L_emem",1000,0,700);
  TH1F *hL1520massAllRL_L_emem=new TH1F("hL1520massAllRL_L_emem","hL1520massAllRL_L_emem",1000,1200,2200);

     
  TH1F *hDLmassDist_emem=new TH1F("hDLmassDist_emem","hDLmassDist_emem",1000,0,700);
  TH1F *hL1520massDist_emem=new TH1F("hL1520massDist_emem","hL1520massDist_emem",1000,1200,2200);

  TH1F *hDLmassDistRL_emem=new TH1F("hDLmassDistRL_emem","hDLmassDistRL_emem",1000,0,700);
  TH1F *hL1520massDistRL_emem=new TH1F("hL1520massDistRL_emem","hL1520massDistRL_emem",1000,1200,2200);

  TH1F *hDLmassDistRL_L_emem=new TH1F("hDLmassDistRL_L_emem","hDLmassDistRL_L_emem",1000,0,700);
  TH1F *hL1520massDistRL_L_emem=new TH1F("hL1520massDistRL_L_emem","hL1520massDistRL_L_emem",1000,1200,2200);
     
  TH1F *hDLmassDistZ_emem=new TH1F("hDLmassDistZ_emem","hDLmassDistZ_emem",1000,0,700);
  TH1F *hL1520massDistZ_emem=new TH1F("hL1520massDistZ_emem","hL1520massDistZ_emem",1000,1200,2200);

  TH1F *hDLmassDistZRL_emem=new TH1F("hDLmassDistZRL_emem","hDLmassDistZRL_emem",1000,0,700);
  TH1F *hL1520massDistZRL_emem=new TH1F("hL1520massDistZRL_emem","hL1520massDistZRL_emem",1000,1200,2200);

  TH1F *hDLmassDistZRL_L_emem=new TH1F("hDLmassDistZRL_L_emem","hDLmassDistZRL_L_emem",1000,0,700);
  TH1F *hL1520massDistZRL_L_emem=new TH1F("hL1520massDistZRL_L_emem","hL1520massDistZRL_L_emem",1000,1200,2200);

     
  TH1F *hDLmassDistL_emem=new TH1F("hDLmassDistL_emem","hDLmassDistL_emem",1000,0,700);
  TH1F *hL1520massDistL_emem=new TH1F("hL1520massDistL_emem","hL1520massDistL_emem",1000,1200,2200);
  TH1F *hL1520massDistLpi0_emem=new TH1F("hL1520massDistLpi0_emem","hL1520massDistLpi0_emem",1000,1200,2200);
     
  TH1F *hDLmassDistLRL_emem=new TH1F("hDLmassDistLRL_emem","hDLmassDistLRL_emem",1000,0,700);
  TH1F *hL1520massDistLRL_emem=new TH1F("hL1520massDistLRL_emem","hL1520massDistLRL_emem",1000,1200,2200);

  TH1F *hDLmassDistLRL_L_emem=new TH1F("hDLmassDistLRL_L_emem","hDLmassDistLRL_L_emem",1000,0,700);
  TH1F *hL1520massDistLRL_L_emem=new TH1F("hL1520massDistLRL_L_emem","hL1520massDistLRL_L_emem",1000,1200,2200);
     
  TH1F *hDLmassDistZL_emem=new TH1F("hDLmassDistZL_emem","hDLmassDistZL_emem",1000,0,700);
  TH1F *hL1520massDistZL_emem=new TH1F("hL1520massDistZL_emem","hL1520massDistZL_emem",1000,1200,2200);

  TH1F *hDLmassDistZLRL_emem=new TH1F("hDLmassDistZLRL_emem","hDLmassDistZLRL_emem",1000,0,700);
  TH1F *hL1520massDistZLRL_emem=new TH1F("hL1520massDistZLRL_emem","hL1520massDistZLRL_emem",1000,1200,2200);
  TH1F *hDLmassDistZLRL_L_emem=new TH1F("hDLmassDistZLRL_L_emem","hDLmassDistZLRL_L_emem",1000,0,700);
  TH1F *hL1520massDistZLRL_L_emem=new TH1F("hL1520massDistZLRL_L_emem","hL1520massDistZLRL_L_emem",1000,1200,2200);


     
  TH1F *hDLmassFinal_emem=new TH1F("hDLmassFinal_emem","hDLmassFinal_emem",1000,0,700);
  TH1F *hL1520massFinal_emem=new TH1F("hL1520massFinal_emem","hL1520massFinal_emem",1000,1200,2200);

  TH1F *hDLmassFinalRL_emem=new TH1F("hDLmassFinalRL_emem","hDLmassFinalRL_emem",1000,0,700);
  TH1F *hL1520massFinalRL_emem=new TH1F("hL1520massFinalRL_emem","hL1520massFinalRL_emem",1000,1200,2200);

  TH1F *hDLmassFinalRL_L_emem=new TH1F("hDLmassFinalRL_L_emem","hDLmassFinalRL_L_emem",1000,0,700);
  TH1F *hL1520massFinalRL_L_emem=new TH1F("hL1520massFinalRL_L_emem","hL1520massFinalRL_L_emem",1000,1200,2200);
     
     
  TH1F *hL1520massFinalRLpi0_emem=new TH1F("hL1520massFinalRLpi0_emem","hL1520massFinalRLpi0_emem",1000,1200,2200);

  TH1F *hL1520massFinalRLpi0_L_emem=new TH1F("hL1520massFinalRLpi0_L_emem","hL1520massFinalRLpi0_L_emem",1000,1200,2200);


  TH1F *hL1520massFinalpi0_emem=new TH1F("hL1520massFinalpi0_emem","hL1520massFinalpi0_emem",1000,1200,2200);
  TH1F *hL1520massDistZLpi0_emem=new TH1F("hL1520massDistZLpi0_emem","hL1520massDistZLpi0_emem",1000,1200,2200);
     
  TH1F *hDLmassDistZL_L=new TH1F("hDLmassDistZL_L","hDLmassDistZL_L",1000,0,700);
  TH1F *hDLmassDistL_L=new TH1F("hDLmassDistL_L","hDLmassDistL_L",1000,0,700);
     
  //**************************************************************************** 

  //**************************************
   
  TH1F *hinvMass_HHemem=new TH1F("hinvMass_HHemem","hinvMass_HHemem",1000,0,700);
  TH1F *hinvMass_HHepep=new TH1F("hinvMass_HHepep","hinvMass_HHepep",1000,0,700);

  TH1F *hinvMass_HFTemem=new TH1F("hinvMass_HFTemem","hinvMass_HFTemem",1000,0,700);
  TH1F *hinvMass_HFTepep=new TH1F("hinvMass_HFTepep","hinvMass_HFTepep",1000,0,700);

  TH1F *hinvMass_emem_L=new TH1F("hinvMass_emem_L","hinvMass_emem_L",1000,0,700);
  TH1F *hinvMass_epep_L=new TH1F("hinvMass_epep_L","hinvMass_epep_L",1000,0,700);

  TH1F *hinvMass_emem=new TH1F("hinvMass_emem","hinvMass_emem",1000,0,700);
  TH1F *hinvMass_epep=new TH1F("hinvMass_epep","hinvMass_epep",1000,0,700);

  TH1F *hL1520mass_HFTemem=new TH1F("hL1520mass_HFTemem","hL1520mass_HFTemem",1000,1200,2200);
  TH1F *hL1520mass_HFTepep=new TH1F("hL1520mass_HFTepep","hL1520mass_HFTepep",1000,1200,2200);

  TH1F *hL1520mass_HHemem=new TH1F("hL1520mass_HHemem","hL1520mass_HHemem",1000,1200,2200);
  TH1F *hL1520mass_HHepep=new TH1F("hL1520mass_HHepep","hL1520mass_HHepep",1000,1200,2200);

  TH1F *hL1520mass_emem=new TH1F("hL1520mass_emem","hL1520mass_emem",1000,1200,2200);
  TH1F *hL1520mass_epep=new TH1F("hL1520mass_epep","hL1520mass_epep",1000,1200,2200);
    
  TH1F *hL1520mass_emem_L=new TH1F("hL1520mass_emem_L","hL1520mass_emem_L",1000,1200,2200);
  TH1F *hL1520mass_epep_L=new TH1F("hL1520mass_epep_L","hL1520mass_epep_L",1000,1200,2200);

     

  TH1F *hinvMass_ememOA=new TH1F("hinvMass_ememOA","hinvMass_ememOA",1000,0,700);
  TH1F *hinvMass_epepOA=new TH1F("hinvMass_epepOA","hinvMass_epepOA",1000,0,700);
  TH1F *hinvMass_epemOA=new TH1F("hinvMass_epemOA","hinvMass_epemOA",1000,0,700);
    
     
  TH1F *hZvertLam1520HH=new TH1F("hZvertLam1520HH","hZvertLam1520HH",500,-500,500);
  TH1F *hdistTgLam1520HH=new TH1F("hdistTgLam1520HH","hdistTgLam1520HH",300,0,300);
  TH1F *hZvertLam1520TgHH=new TH1F("hZvertLam1520TgHH","hZvertLam1520TgHH",500,-500,500);

  TH1F *hZvertLamHH=new TH1F("hZvertLamHH","hZvertLamHH",500,-500,500);
  TH1F *hZvertHHAll=new TH1F("hZvertHHAll","hZvertHHAll",500,-500,500);
  TH1F *hZvertHHRL=new TH1F("hZvertHHRL","hZvertHHRL",500,-500,500);

  TH1F *hdistTgLamHH=new TH1F("hdistTgLamHH","hdistTgLamHH",300,0,300);
  TH1F *hZvertLamTgHH=new TH1F("hZvertLamTgHH","hZvertLamTgHH",500,-500,500);

  TH1F *hZvertLam1520FT=new TH1F("hZvertLam1520FT","hZvertLam1520FT",500,-500,500);
  TH1F *hdistTgLam1520FT=new TH1F("hdistTgLam1520FT","hdistTgLam1520FT",300,0,300);
  TH1F *hZvertLam1520TgFT=new TH1F("hZvertLam1520TgFT","hZvertLam1520TgFT",500,-500,500);

  TH1F *hZvertLamFT=new TH1F("hZvertLamFT","hZvertLamFT",500,-500,500);
  TH1F *hZvertFTAll=new TH1F("hZvertFTAll","hZvertFTAll",500,-500,500);
  TH1F *hZvertFTRL=new TH1F("hZvertFTRL","hZvertFTRL",500,-500,500);

  TH1F *hdistTgLamFT=new TH1F("hdistTgLamFT","hdistTgLamFT",300,0,300);
  TH1F *hZvertLamTgFT=new TH1F("hZvertLamTgFT","hZvertLamTgFT",500,-500,500);
     
  //*************************************************************************
  TCanvas* cEff=new TCanvas("cEff","detection_effi");
  TH1F* hEprotons4Pi=new TH1F("hEprotons4Pi","prot_kine;theta",180,0,180);
  TH1F* hEprotonsdet=new TH1F("hEprotonsdet","prot_registered;theta",180,0,180);
  TH1F* hEprotonsEff=new TH1F("hEprotonsEff","prot_effi;theta",180,0,180);
  TH1F* hEpions4Pi=new TH1F("hEpions4Pi","pim_kine; theta",180,0,180);
  TH1F* hEpionsdet=new TH1F("hEpionsdet","pim_registered; theta",180,0,180);
  TH1F* hEpionsEff=new TH1F("hEpionsEff","pim_effi; theta;efficiency",180,0,180);
  TH1F* hEleptons4Pi=new TH1F("hEleptons4Pi","leptons_kine; theta",180,0,180);
  TH1F* hEleptonsdet=new TH1F("hEleptonsdet","leptons_registered; theta",180,0,180);
  TH1F* hEleptonsEff=new TH1F("hEleptonsEff","lept_effi; theta;efficiency",180,0,180);
  TH2F* h2Eproton4Pi=new TH2F("h2Eproton4Pi","prot_4Pi; phi; theta",100,0,360,50,0,180);
  TH2F* h2Eprotondet=new TH2F("h2Eprotondet","prot_registered; phi; theta",100,0,360,50,0,180);
  //TH2F* h2EprotonEff=new TH2F("h2EprotonEff","prot_effi_4Pi; phi; theta",100,0,360,50,0,180);
  TH2F* h2Epion4Pi=new TH2F("h2Epion4Pi","pim_4Pi; phi; theta",100,0,360,50,0,180);
  TH2F* h2Epiondet=new TH2F("h2Epiondet","pim_registered; phi; theta",100,0,360,50,0,180);
  //TH2F* h2EpionEff=new TH2F("h2EpionEff","pim_effi_4Pi; phi; theta",100,0,360,50,0,180);
  TH2F* h2Elepton4Pi=new TH2F("h2Elepton4Pi","dil_4Pi; phi; theta",100,0,360,50,0,180);
  TH2F* h2Eleptondet=new TH2F("h2Eleptondet","dil_registered; phi; theta",100,0,360,50,0,180);

   
  TCanvas* cEffFromLambda=new TCanvas("cEffFromLambda","detection_effi_Lambda1115");
  TH1F* hEFLprotons4Pi=new TH1F("hEFLprotons4Pi","prot_kine; theta",180,0,180);
  TH1F* hEFLprotonsdet=new TH1F("hEFLprotonsdet","prot_registered; theta",180,0,180);
  TH1F* hEFLprotonsEff=new TH1F("hEFLprotonsEff","prot_effi; theta;efficiency",180,0,180);
  TH1F* hEFLpions4Pi=new TH1F("hEFLprotons4Pi","prot_kine; theta",180,0,180);
  TH1F* hEFLpionsdet=new TH1F("hEFLprotonsdet","prot_registered; theta",180,0,180);
  TH1F* hEFLpionsEff=new TH1F("hEFLprotonsEff","prot_effi; theta;efficiency",180,0,180);
     
  TH2F *proton_identyf_real=new TH2F("proton_identyf_real","proton_identyf_real;p[MeV];#beta",1000,0,3000,1000,0,1.2);
  TH2F *proton_identyf_ideal=new TH2F("proton_identyf_ideal","proton_identyf_ideal;p[MeV];#beta",1000,0,3000,1000,0,1.2);
  TH2F *pion_identyf_real=new TH2F("pion_identyf_real","pion_identyf_real;p[MeV];#beta",1000,0,2000,1000,0,1.2);
  TH2F *pion_identyf_ideal=new TH2F("pion_identyf_ideal","pion_identyf_ideal;p[MeV];#beta",1000,0,2000,1000,0,1.2);
  TH2F *lepton_identyf_real=new TH2F("lepton_identyf_real","lepton_identyf_real;p[MeV];#beta",1000,0,2000,1000,0,1.2);
  TH2F *lepton_identyf_ideal=new TH2F("lepton_identyf_ideal","lepton_identyf_ideal;p[MeV];#beta",1000,0,2000,1000,0,1.2);

  TH2F *lepton_identyf_real_dphi=new TH2F("lepton_identyf_real_dphi","lepton_identyf_real_dphi;dphi;p[MeV]",200,-10,10,1000,0,2000);
  TH2F *lepton_identyf_ideal_dphi=new TH2F("lepton_identyf_ideal_dphi","lepton_identyf_ideal_dphi;dphi;p[MeV]",200,-10,10,1000,0,2000);
  TH2F *lepton_identyf_real_dtheta=new TH2F("lepton_identyf_real_dtheta","lepton_identyf_real_dtheta;dtheta;p[MeV]",200,-10,10,1000,0,2000);
  TH2F *lepton_identyf_ideal_dtheta=new TH2F("lepton_identyf_ideal_dtheta","lepton_identyf_ideal_dtheta;dtheta;p[MeV]",200,-10,10,1000,0,2000);

  
  //********************************************

  //TH1F* hinvMepemFT =new TH1F("hinvMepemFT","hinvMepemFT",1000,0,700);
  //TH1F* hinvMepemHH =new TH1F("hinvMepemHH","hinvMepemHH",1000,0,700);
     
  //**************************************************

  float ww=0;   
  Int_t entries = loop->getEntries();
  //     //setting numbers of events regarding the input number of events by the user
  if (anapars.events < entries and anapars.events >= 0 ) entries = anapars.events;


    
  TFile * output_file = TFile::Open(anapars.outfile, "RECREATE");
  output_file->cd();
    
  cout << "NEW ROOT TREE , vertex_ana" << endl;
  Int_t evnb;
    
  //  if (fChain == 0) return;
    
  vector<HParticleCandSim*> lH, pH, pimH,ep,em;
  vector<HFwDetCandSim*> pFT;
 
  int flagHHL1=0;
  int flagHHL2=0;
  int flagHHL3=0;
  int flagHHL4=0;
  
  int flagHFTL1=0;
  int flagHFTL2=0;
  int flagHFTL3=0;
  int flagHFTL4=0;
  


  for (Long_t event=0; event<entries; event++) 
    {

      loop->nextEvent(event); 
      if(!(event%10000))   cout<<"event no. "<<event<<endl;

    
      evnb=event;
      //fChain->GetEntry(event);
      //cout<< eventNum<<endl;
      lH.clear();
      pH.clear();  
      pFT.clear();  
      pimH.clear();
      em.clear();    
      ep.clear();

      flagHHL1=0;
      flagHHL2=0;
      flagHHL3=0;
      flagHHL4=0;

      
      flagHFTL1=0;
      flagHFTL2=0;
      flagHFTL3=0;
      flagHFTL4=0;

      
      int hnum=particleCatSim->getEntries();
      int fnum=fwDetCatSim->getEntries();
      int knum=fCatGeantKine->getEntries();
      //HParticleCandSim* ep=nullptr;
      //HParticleCandSim* em=nullptr;

      HParticleCandSim* pionH=nullptr;
   
      HParticleCandSim* protonH=nullptr;
      HFwDetCandSim*  protonFT=nullptr;

      HParticleCandSim* partH=nullptr;
      HFwDetCandSim*  partFT=nullptr;


      HGeantKine* kine=nullptr;
      Int_t isDilepton=0;
      Int_t isLambda=0;
      Int_t isLambda1520=0;
      Int_t leptonparticle1=-1;
      Int_t leptonparticle2=-1;
      Int_t lambdaparticle1=-1;
      Int_t lambdaparticle2=-1;
      Int_t leptonparticle=-1;
      HParticleTool tool;
	
      HGeomVector vertexL;
      HGeomVector vertexDL;
      HGeomVector vertexL1520;
      HGeomVector dirL;
      HGeomVector dirDL;
      HGeomVector dirL1520,dirL1520_1;

	
      double min_dist_dl=20; //[mm]
      double min_dist_l=20;
      double min_angle=4;//e+ e- opening angle [deg]

      HGeomVector base_Tg, dir_Tg;
			  
      base_Tg.setX(0);
      base_Tg.setY(0);
      base_Tg.setZ(-25.);
      dir_Tg.setX(0);
      dir_Tg.setY(0);
      dir_Tg.setZ(1);

      HGeomVector ver_L1520Tg, ver_L1520TgFT, baseL1520, ver_LTg, ver_LTgFT;
			 


	
      if (hnum){
	//cout<<"XXXXXXX"<<hnum<<endl;
   
	//HADES
	for (int i=0;i<hnum;i++){
	 
	  // if(pid==9)thmom9_4pi->Fill(mom,th);
	  //if(pid==14)thmom14_4pi->Fill(mom,th);
	  //if(pid==2)thmom2_4pi->Fill(mom,th);
	  //if(pid==3)thmom3_4pi->Fill(mom,th);
	  //if(pid==9)thmom9_acc->Fill(mom,th,ae);
	  //if(pid==14)thmom14_acc->Fill(mom,th,ae);
	  partH=HCategoryManager::getObject(partH, particleCatSim,i);
	
	  //lep=HCategoryManager::getObject(lep, particleCatSim,i);
	  //pion=HCategoryManager::getObject(pion, particleCatSim,i);
	  //protonFT=HCategoryManager::getObject(protonFT,fwDetCatSim,i);

	  float nerbyFit=0;
	  float nerbyUnFit=0;
	  
	  //nerbyFit=fabs(partH->getAngleToNearbyFittedInner())*TMath::RadToDeg(); 
	  //nerbyUnFit=fabs(partH->getAngleToNearbyUnfittedInner())*TMath::RadToDeg(); 
	  nerbyFit=fabs(partH->getAngleToNearbyFittedInner()); 
	  nerbyUnFit=fabs(partH->getAngleToNearbyUnfittedInner()); 

	  //cout<<"---> "<<nerbyFit<<" "<<nerbyUnFit<<endl;
	  if(partH->getRichMatchingQuality()!=-1 && partH->isFlagBit(kIsUsed) && nerbyFit>4.  && nerbyUnFit>4. )
	    {
	      //if(partH->getRichMatchingQuality()!=-1 && partH->isFlagBit(kIsUsed)){
	      //if(partH->getRichMatchingQuality()!=-1 && partH->isFlagBit(kIsUsed) && nerbyFit>2.  && nerbyUnFit>2. ){
	      int flagdil=0;
	   
	      if(//partH->getGeantPID()==2
		 isLepton(partH,lepton_identyf_ideal,lepton_identyf_ideal_dphi,lepton_identyf_ideal_dtheta,lepton_identyf_real,lepton_identyf_real_dphi,lepton_identyf_real_dtheta)
		 && partH->getCharge()==1
		 ) //e+
		{
		  partH->calc4vectorProperties(HPhysicsConstants::mass(2));
		  ep.push_back(partH);  
		  hEleptonsdet->Fill(partH->getTheta());
		  h2Eleptondet->Fill(partH->getPhi(),partH->getTheta());
		  //flagdil=1;
		  //h2IIleptonsInAcceptance->Fill(partH->getMomentum(),partH->getTheta());
		  //cout<<"e+ "<<evnb<<" "<<partH->getGeantParentPID()<<" "<<partH->getGeantParentTrackNum()<<" "<<partH->getGeantGrandParentPID()<<endl;
		  //dobre pary: 0,-1
		}
	   		  
	      if(//partH->getGeantPID()==3
		 isLepton(partH,lepton_identyf_ideal,lepton_identyf_ideal_dphi,lepton_identyf_ideal_dtheta,lepton_identyf_real,lepton_identyf_real_dphi,lepton_identyf_real_dtheta)
		 && partH->getCharge()==-1
		 ) //e-
		{
		  partH->calc4vectorProperties(HPhysicsConstants::mass(3));
		  em.push_back(partH);  
		  hEleptonsdet->Fill(partH->getTheta());
		  //if(flagdil)
		  h2Eleptondet->Fill(partH->getPhi(),partH->getTheta());
		  //cout<<"e- "<<evnb<<" "<<partH->getGeantParentPID()<<" "<<partH->getGeantParentTrackNum()<<" "<<partH->getGeantGrandParentPID()<<endl;
		}
	    }
	  
	  if(//partH->getGeantPID()==9 && partH->isFlagBit(kIsUsed)
	     isPion(partH,pion_identyf_ideal,pion_identyf_real) && partH->getCharge()==-1
	     )// proton->getChi2()<10)
	    {  
	      partH->calc4vectorProperties(HPhysicsConstants::mass(9));
	      pimH.push_back(partH);  

	      hEpionsdet->Fill(partH->getTheta());
	      h2Epiondet->Fill(partH->getPhi(),partH->getTheta());

	      if(partH->getGeantParentPID()==18) hEFLpionsdet->Fill(partH->getTheta());
	   	      
	    }

	  if(//partH->getGeantPID()==14 && partH->isFlagBit(kIsUsed)
	     isProton(partH,proton_identyf_ideal,proton_identyf_real) && partH->getCharge()==1
	     )// proton->getChi2()<10)
	    {
	   
	      partH->calc4vectorProperties(HPhysicsConstants::mass(14));
	      pH.push_back(partH);  


	      hEprotonsdet->Fill(partH->getTheta());
	      h2Eprotondet->Fill(partH->getPhi(),partH->getTheta());
    
	      if(partH->getGeantParentPID()==18) hEFLpionsdet->Fill(partH->getTheta());
	    }
	}		 
      }	 
      //FT
	
      if (fnum){
	// cout<<"0:: "<<fnum<<endl;	  
	for (int i=0;i<fnum;i++){

	  partFT=HCategoryManager::getObject(partFT,fwDetCatSim,i);

	  //if(partFT->getGeantPID()==14 || partFT->getGeantPID()==9 || partFT->getGeantPID()==11 ){
	  //if(partFT->getGeantPID()==14 && partFT->getGeantParentTrackNum()==0){

	  partFT->calc4vectorProperties(HPhysicsConstants::mass(14));
	  pFT.push_back(partFT);  
	  //cout<<" "<<partFT->getTheta()<<endl;

	  if(partFT->getGeantPID()==14){
	   
	    hEprotonsdet->Fill(partFT->getTheta());
	    h2Eprotondet->Fill(partFT->getPhi(),partFT->getTheta());
    
	    if(partFT->getGeantParentPID()==18) hEFLpionsdet->Fill(partFT->getTheta());
	  }

	  if(partFT->getGeantPID()==9){
		
	    hEpionsdet->Fill(partFT->getTheta());
	    h2Epiondet->Fill(partFT->getPhi(),partFT->getTheta());
		
	    if(partFT->getGeantParentPID()==18) hEFLpionsdet->Fill(partFT->getTheta());
	  }
	      
	  //}
	    
	}
	  
      }
      //**************
      //HADES-HADES
      //*************
      //cout<<"::::"<<pH.size()<<" "<<pimH.size()<<endl;

      if (pH.size()>=1 && pimH.size()>=1)
	{
	  for (int k=0;k<pH.size();k++){
	    for (int j=0;j<pimH.size();j++){

	      //pimH[j]->calc4vectorProperties(HPhysicsConstants::mass(pimH[j]->getGeantPID()));
	      //pH[k]->calc4vectorProperties(HPhysicsConstants::mass(14));
	      ww=0;
	      ww=pimH[j]->getGeantGenweight();
	      //cout<<ww<<endl;
	      double lambdaM=(*pH[k]+*pimH[j]).M();
	      double lambdaD=trackDistance(pH[k],pimH[j]);
	      //cout<<"----------------->>> "<<lambdaM<<" "<<lambdaD<<endl;

	      vertexL=trackVertex(pimH[j],pH[k]);
	      dirL.setXYZ((*pimH[j]+*pH[k]).X(),(*pimH[j]+*pH[k]).Y(),(*pimH[j]+*pH[k]).Z());

	      ver_LTg=tool.calcVertexAnalytical(base_Tg,dir_Tg,vertexL,dirL);
	      double distLamZ=tool.calculateMinimumDistance(base_Tg, dir_Tg,vertexL,dirL);
			  
			  
	      hZvertHHAll->Fill(vertexL.Z());//***
	      //hdistTgLamHH->Fill(distLamZ);
	      //hZvertLamTgHH->Fill(ver_LTg.Z());
	      
	      hinvM_pmHpHAll->Fill(lambdaM,ww);
	      hdist_pmHpHAll->Fill(lambdaD,ww);
	      
	      if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18)
		{
		  hZvertHHRL->Fill(vertexL.Z());
		  hLRdist->Fill(lambdaD,ww);
		  hLRmass->Fill(lambdaM,ww);
		     
		}

	      flagHHL1=0;
	      flagHHL2=0;
	      flagHHL3=0;
	      flagHHL4=0;

	      if (lambdaD<min_dist_l)
		{

		  flagHHL1=1;
		  hinvM_pmHpHDist->Fill(lambdaM,ww);
		  if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18)
		    hLRmassDist->Fill(lambdaM,ww);
		  //hLHmassDist->Fill(lambdaM,ww);

		  h2LHvertex->Fill(vertexL.Z(),TMath::Sqrt(vertexL.X()*vertexL.X()+vertexL.Y()*vertexL.Y()));
		}

	      if (lambdaD<min_dist_l && vertexL.Z()>0.)
		{
		  flagHHL2=1;
		  hinvM_pmHpHDistZ->Fill(lambdaM,ww);
		  if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18) hLRmassDistZ->Fill(lambdaM,ww);

		}
	  
	      if (lambdaM>1105 && lambdaM<1125 && lambdaD<min_dist_l)
		{
		  flagHHL3=1;

		  hinvM_pmHpHDistL->Fill(lambdaM,ww);
		  if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18) hLRmassDistL->Fill(lambdaM,ww);
		
		  hZvertLamHH->Fill(vertexL.Z());
		  hdistTgLamHH->Fill(distLamZ);
		  hZvertLamTgHH->Fill(ver_LTg.Z());

		}
	      if (lambdaM>1105 && lambdaM<1125 && lambdaD<min_dist_l && vertexL.Z()>0.)
		{
		  flagHHL4=1;
	      
		  hinvM_pmHpHDistZL->Fill(lambdaM,ww);
		  if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18) hLRmassDistZL->Fill(lambdaM,ww);
		}

	      if (ep.size() || em.size()){
		for (int s=0;s<ep.size();s++){
		  for (int ss=0;ss<ep.size();ss++){

		    TLorentzVector lvLambda=*pH[k]+*pimH[j];
		    TLorentzVector lvDiLepton=*ep[s]+*ep[ss];
		    double oa = tool.getOpeningAngle(ep[s],ep[ss]);
		    double mass_1520=(lvLambda+lvDiLepton).M();

		    int flagDil2=1;

		    //if(ep[s]->getGeantParentTrackNum()==0 && ep[ss]->getGeantParentTrackNum()==0 && ep[ss]->getGeantGrandParentPID()==-1 && ep[s]->getGeantGrandParentPID()==-1) flagDil2=1;
		    //if(ep[s]->getGeantParentTrackNum()== ep[ss]->getGeantParentTrackNum()) flagDil2=1;

		    
		    //***************************************

		    if(oa>min_angle)
		      {
			hinvMass_epepOA->Fill(lvDiLepton.M(),ww);  
			hDLmassAll_epep->Fill(lvDiLepton.M(),ww);
			hL1520massAll_epep->Fill(mass_1520,ww);

			if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18)
			  {
			    hDLmassAllRL_epep->Fill(lvDiLepton.M(),ww);
			    hL1520massAllRL_epep->Fill(mass_1520,ww);
			    if(flagDil2)
			      {
				hDLmassAllRL_L_epep->Fill(lvDiLepton.M(),ww);
				hL1520massAllRL_L_epep->Fill(mass_1520,ww);
			
			      }
			  }

			if(flagHFTL1){
			  hDLmassDist_epep->Fill(lvDiLepton.M(),ww);
			  hL1520massDist_epep->Fill(mass_1520,ww);

			  if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){
			    hDLmassDistRL_epep->Fill(lvDiLepton.M(),ww);
			    hL1520massDistRL_epep->Fill(mass_1520,ww);

			    if(flagDil2){
			      hDLmassDistRL_L_epep->Fill(lvDiLepton.M(),ww);
			      hL1520massDistRL_L_epep->Fill(mass_1520,ww);

			    }
      
			  }
			}

			if(flagHFTL2){
			  hDLmassDistZ_epep->Fill(lvDiLepton.M(),ww);
			  hL1520massDistZ_epep->Fill(mass_1520,ww);

			  if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){
			    hDLmassDistZRL_epep->Fill(lvDiLepton.M(),ww);
			    hL1520massDistZRL_epep->Fill(mass_1520,ww);

			    if(flagDil2){
			      hDLmassDistZRL_L_epep->Fill(lvDiLepton.M(),ww);
			      hL1520massDistZRL_L_epep->Fill(mass_1520,ww);
			  
			    }

			

			  }
			}


		    
			if(flagHFTL3){

			  hDLmassDistL_epep->Fill(lvDiLepton.M(),ww);
			  hL1520massDistL_epep->Fill(mass_1520,ww);
			  if(lvDiLepton.M()>140. && lvDiLepton.M()>420.)hL1520massDistLpi0_epep->Fill(mass_1520,ww);

			  if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){

			    hDLmassDistLRL_epep->Fill(lvDiLepton.M(),ww);
			    hL1520massDistLRL_epep->Fill(mass_1520,ww);

			    if(flagDil2){
			      hDLmassDistLRL_L_epep->Fill(lvDiLepton.M(),ww);
			      hL1520massDistLRL_L_epep->Fill(mass_1520,ww);
			    }
			
			  }



			}

		    		    
		      }//oa

		    //***************************************

	    
		    if(flagHHL4 && oa>min_angle){

		      hDLmassDistZL_epep->Fill(lvDiLepton.M(),ww);

		      hL1520mass_HHepep->Fill(mass_1520,ww);		      
		      if(mass_1520>1450 && mass_1520<1550)hinvMass_HHepep->Fill(lvDiLepton.M(),ww);

		      if(flagDil2){

			hL1520mass_epep_L->Fill(mass_1520,ww);		      
			if(mass_1520>1450 && mass_1520<1550)hinvMass_epep_L->Fill(lvDiLepton.M(),ww);
			if(lvDiLepton.M()>140. && lvDiLepton.M()<420.)
			  hL1520massDistZLpi0_epep->Fill(mass_1520,ww);
		      

		      }

		      
		    }
		    
		  }
		}
		
		for (int s=0;s<em.size();s++){
		  for (int ss=0;ss<em.size();ss++){

		    TLorentzVector lvLambda=*pH[k]+*pimH[j];
		    TLorentzVector lvDiLepton=*em[s]+*em[ss];
		    double oa = tool.getOpeningAngle(em[s],em[ss]);
		    double mass_1520=(lvLambda+lvDiLepton).M();

		    int flagDil1=1;
		    //if(em[s]->getGeantParentTrackNum()==0 && em[ss]->getGeantParentTrackNum()==0 && em[ss]->getGeantGrandParentPID()==-1 && em[s]->getGeantGrandParentPID()==-1) flagDil1=1;
		    //if(em[s]->getGeantParentTrackNum()== em[ss]->getGeantParentTrackNum()) flagDil1=1;

		    
		    //***************************************

		    if(oa>min_angle){
		      hinvMass_ememOA->Fill(lvDiLepton.M(),ww);  
		      
		      hDLmassAll_emem->Fill(lvDiLepton.M(),ww);
		      hL1520massAll_emem->Fill(mass_1520,ww);

		      if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){

			hDLmassAllRL_emem->Fill(lvDiLepton.M(),ww);
			hL1520massAllRL_emem->Fill(mass_1520,ww);
		      
			if(flagDil1){
			  hDLmassAllRL_L_emem->Fill(lvDiLepton.M(),ww);
			  hL1520massAllRL_L_emem->Fill(mass_1520,ww);
			
			}
		      
		      }

		      if(flagHFTL1){

			hDLmassDist_emem->Fill(lvDiLepton.M(),ww);
			hL1520massDist_emem->Fill(mass_1520,ww);

			if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){
			  hDLmassDistRL_emem->Fill(lvDiLepton.M(),ww);
			  hL1520massDistRL_emem->Fill(mass_1520,ww);

			  if(flagDil1){
			    hDLmassDistRL_L_emem->Fill(lvDiLepton.M(),ww);
			    hL1520massDistRL_L_emem->Fill(mass_1520,ww);

			  }
      
			}

		      }

		      if(flagHFTL2){
			hDLmassDistZ_emem->Fill(lvDiLepton.M(),ww);
			hL1520massDistZ_emem->Fill(mass_1520,ww);

			if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){
			  hDLmassDistZRL_emem->Fill(lvDiLepton.M(),ww);
			  hL1520massDistZRL_emem->Fill(mass_1520,ww);

			  if(flagDil1){
			    hDLmassDistZRL_L_emem->Fill(lvDiLepton.M(),ww);
			    hL1520massDistZRL_L_emem->Fill(mass_1520,ww);
			  
			  }

			

			}
		      }


		    
		      if(flagHFTL3){

			hDLmassDistL_emem->Fill(lvDiLepton.M(),ww);
			hL1520massDistL_emem->Fill(mass_1520,ww);
			if(lvDiLepton.M()>140.&& lvDiLepton.M()<420.)hL1520massDistLpi0_emem->Fill(mass_1520,ww);

			if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){

			  hDLmassDistLRL_emem->Fill(lvDiLepton.M(),ww);
			  hL1520massDistLRL_emem->Fill(mass_1520,ww);

			  if(flagDil1){
			    hDLmassDistLRL_L_emem->Fill(lvDiLepton.M(),ww);
			    hL1520massDistLRL_L_emem->Fill(mass_1520,ww);
			  }
			
			}



		      }

		    		    
		    }//oa

		    //***************************************
		    
		    if(flagHHL4 && oa>min_angle){
		      hDLmassDistZL_emem->Fill(lvDiLepton.M(),ww);

		      hL1520mass_emem->Fill(mass_1520,ww);		      
		      if(mass_1520>1450 && mass_1520<1550)hinvMass_emem->Fill(lvDiLepton.M(),ww);

		      if(flagDil1){

			hL1520mass_emem_L->Fill(mass_1520,ww);		      
			if(mass_1520>1450 && mass_1520<1550)hinvMass_emem_L->Fill(lvDiLepton.M(),ww);

			//hL1520massFinalpi0_emem
			if(lvDiLepton.M()>140. && lvDiLepton.M()<420.)
			  hL1520massDistZLpi0_emem->Fill(mass_1520,ww);

			
		      }
		      
		      
		    }
		    
		  }
		}
		
	      }
	

	      
	      if (ep.size() && em.size()){
		//cout<<":::::::::::epem "<<ep.size()<<" "<<em.size()<<endl;
		for (int s=0;s<ep.size();s++){
		  for (int ss=0;ss<em.size();ss++){

 
		    //em[ss]->calc4vectorProperties(HPhysicsConstants::mass(em[ss]->getGeantPID()));
		    //ep[s]->calc4vectorProperties(HPhysicsConstants::mass(ep[s]->getGeantPID()));
		    int flagDil=0;
		    
		    //if(ep[s]->getGeantParentTrackNum()==0 && em[ss]->getGeantParentTrackNum()==0 && ep[s]->getGeantGrandParentPID()==-1 && em[ss]->getGeantGrandParentPID()==-1) flagDil=1;
		    if(ep[s]->getGeantParentTrackNum()== em[ss]->getGeantParentTrackNum() /*&& em[s]->getGeantParentPID()==7*/)
		      flagDil=1;
		    //cout<<"yyy "<<   ep[s]->getGeantParentTrackNum()<<" "<< em[ss]->getGeantParentTrackNum()<<endl;
		   
		    
		    TLorentzVector lvLambda=*pH[k]+*pimH[j];
		    TLorentzVector lvDiLepton=*ep[s]+*em[ss];

		    vertexDL=trackVertex(ep[s],em[ss]);
		    dirDL.setXYZ((*ep[s]+*em[ss]).X(),(*ep[s]+*em[ss]).Y(),(*ep[s]+*em[ss]).Z());

		    

		    double mass_1520=(lvLambda+lvDiLepton).M();
		    //min dist between dilepton and Lam1115
		    double distance_1520=tool.calculateMinimumDistance(vertexL,dirL,vertexDL,dirDL);
		    
		    //double invMdilLam=lvLambda.M();
		    double invMepem= lvDiLepton.M(); 

		    int oaFlag=0;		  
		    double oa = tool.getOpeningAngle(ep[s],em[ss]);
		    double dilTrDist=trackDistance(ep[s],em[ss]);
		    if(oa>min_angle)oaFlag=1;


		    if(oaFlag){
		      hinvMass_epemOA->Fill(invMepem,ww);  

		      hDLmassAll->Fill(invMepem,ww);
		      hL1520massAll->Fill(mass_1520,ww);
		      //if(lambdaD<min_dist_l){
		      //hDLmassDist->Fill(invMepem,ww);
		      //}
		      if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){
			hDLmassAllRL->Fill(invMepem,ww);
			hL1520massAllRL->Fill(mass_1520,ww);

			if(flagDil){
			  hDLmassAllRL_L->Fill(invMepem,ww);
			  hL1520massAllRL_L->Fill(mass_1520,ww);
		      
			}
		      }
		    
		      if(flagHHL1){
			hDLmassDist->Fill(invMepem,ww);
			hL1520massDist->Fill(mass_1520,ww);
			if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){
			  hDLmassDistRL->Fill(invMepem,ww);
			  hL1520massDistRL->Fill(mass_1520,ww);

			  if(flagDil){
			    hDLmassDistRL_L->Fill(invMepem,ww);
			    hL1520massDistRL_L->Fill(mass_1520,ww);

			  }
			
			}
		      }


		      if(flagHHL2){
			hDLmassDistZ->Fill(invMepem,ww);
			hL1520massDistZ->Fill(mass_1520,ww);

			if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){

			  hDLmassDistZRL->Fill(invMepem,ww);
			  hL1520massDistZRL->Fill(mass_1520,ww);

			  if(flagDil){
			    hDLmassDistZRL_L->Fill(invMepem,ww);
			    hL1520massDistZRL_L->Fill(mass_1520,ww);
			  
			  }

			
			}

		      }


		      if(flagHHL3){
		      
			hDLmassDistL->Fill(invMepem,ww);
			hL1520massDistL->Fill(mass_1520,ww);
			if(flagDil){
			  hDLmassDistL_L->Fill(invMepem,ww);
			}
			if(invMepem>140. && invMepem<420.)hL1520massDistLpi0->Fill(mass_1520,ww);

			if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){

			  hDLmassDistLRL->Fill(invMepem,ww);
			  hL1520massDistLRL->Fill(mass_1520,ww);
			  if(flagDil){
			    hDLmassDistLRL_L->Fill(invMepem,ww);
			    hL1520massDistLRL_L->Fill(mass_1520,ww);
			  }
			
			}

		      }
		    
		      if(flagHHL4)
			{
			  hDLmassDistZL->Fill(invMepem,ww);
			  hL1520massDistZL->Fill(mass_1520,ww);
			  if(flagDil)
			    {
			      hDLmassDistZL_L->Fill(invMepem,ww);
			    }
			  if(invMepem>140. && invMepem<420.)
			    {
			      hL1520massDistZLpi0->Fill(mass_1520,ww);
			      if(flagDil)
				{
				  hL1520massDistZLpi0_L->Fill(mass_1520,ww);
				  //cout<<"xxxx: "<<mass_1520<<endl;
				}
			    }
			  if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18){

			    hDLmassDistZLRL->Fill(invMepem,ww);
			    hL1520massDistZLRL->Fill(mass_1520,ww);
			    if(invMepem>140. && invMepem<420.)hL1520massDistZLRLpi0->Fill(mass_1520,ww);
		   
			    if(flagDil)
			      {
				hDLmassDistZLRL_L->Fill(invMepem,ww);
				hL1520massDistZLRL_L->Fill(mass_1520,ww);
				if(invMepem>140. && invMepem<420.)hL1520massDistZLRLpi0_L->Fill(mass_1520,ww);
			      }
			  }
		      
			  //if(oa>min_angle){
			  // hL1520massDistOAL->Fill(mass_1520,ww);
			  
			  if(mass_1520>1450 && mass_1520<1550)
			    {

			      hDLmassFinal->Fill(invMepem,ww);
			      hL1520massFinal->Fill(mass_1520,ww);
			      if(invMepem>140. && invMepem<420.)
				hL1520massFinalpi0->Fill(mass_1520,ww);

			      if(pimH[j]->getGeantParentPID()==18 && pH[k]->getGeantParentPID()==18)
				{
				  hDLmassFinalRL->Fill(invMepem,ww);
				  hL1520massFinalRL->Fill(mass_1520,ww);
				  if(invMepem>140. && invMepem<420.)hL1520massFinalRLpi0->Fill(mass_1520,ww);

				  if(flagDil)
				    {
				      hDLmassFinalRL_L->Fill(invMepem,ww);
				      hL1520massFinalRL_L->Fill(mass_1520,ww);
				      if(invMepem>140. && invMepem<420.)
					{
					  hL1520massFinalRLpi0_L->Fill(mass_1520,ww);
					  //cout<<"xxxxx"<<endl;
					}
				
				    }
				}
			      h2L1520vertex->Fill(vertexL1520.Z(),TMath::Sqrt(vertexL1520.X()*vertexL1520.X()+vertexL1520.Y()*vertexL1520.Y()));
			  
			      TLorentzVector lvLam1520=lvLambda+lvDiLepton;
			      //HParticleCandSim *tr;
			      HParticleTool p_tool, tool,tool1;

			      double distance_1520=tool1.calculateMinimumDistance(vertexL,dirL,vertexDL,dirDL);
		
			      vertexL1520=tool.calcVertexAnalytical(vertexL,dirL,vertexDL,dirDL);
			      dirL1520.setXYZ((*ep[s]+*em[ss]+*pimH[j]+*pH[k]).X(),(*ep[s]+*em[ss]+*pimH[j]+*pH[k]).Y(),(*ep[s]+*em[ss]+*pimH[j]+*pH[k]).Z());
			      //dirL1520.setXYZ(lvLam1520.X(),lvLam1520.Y(),lvLam1520.Z());
			
			      //p_tool.calcSegVector(tr->getZ(),tr->getR(),TMath::DegToRad()*tr->getPhi(),TMath::DegToRad()*tr->getTheta(),baseL1520,dirL1520_1);
			  	  
			      //ver_L1520Tg = p_tool.calcVertexAnalytical(base_Tg,dir_Tg,baseL1520,dirL1520);
			      ver_L1520Tg=tool.calcVertexAnalytical(base_Tg,dir_Tg,vertexL1520,dirL1520);

			      //if(vertexL1520.Z()>0){
			      //cout<<"----------------->>1 "<<s<<" "<<ss<<" "<<vertexL1520.Z()<<endl;
			      //cout<<"----------------->>2 "<<s<<" "<<ss<<" "<<ver_L1520Tg.Z()<<endl;
			  		  
			      double distLam1520Z=tool.calculateMinimumDistance(base_Tg, dir_Tg,vertexL1520,dirL1520);
			      //cout<<"distLam1520Z: "<<distLam1520Z<<endl;
			  
			      hZvertLam1520HH->Fill(vertexL1520.Z());
			      hZvertLam1520TgHH->Fill(ver_L1520Tg.Z());
			      hdistTgLam1520HH->Fill(distLam1520Z);
			  
			      //}
			  
			  
			      //if(invMdilLam>1400 && invMdilLam<1700)hDLmassDistOALcut->Fill(invMepem);	    
			  
			  
			    }
			}//flagHHL4
		    }//oa
	
		  }
		}
	      }
	    }
	    
	  }//end of Hades-Hades
	}
      //************** HADES - FT ********************
      //cout<<"::: "<<pimH.size()<<" "<<pFT.size()<<endl;
	
      if (pimH.size()>=1 && pFT.size()>=1){
	for (int k=0;k<pimH.size();k++){
	  for (int j=0;j<pFT.size();j++){

	    ww=0.;
	    ww=pimH[k]->getGeantGenweight();
	    //cout<<ww<<endl;
	    double lambdaM=(*pimH[k]+*pFT[j]).M();
	    double lambdaD=trackDistance(pimH[k],pFT[j]);
	      
	    vertexL=trackVertex(pimH[k],pFT[j]);
	    dirL.setXYZ((*pimH[k]+*pFT[j]).X(),(*pimH[k]+*pFT[j]).Y(),(*pimH[k]+*pFT[j]).Z());
			      
	    hinvM_pmHpFTAll->Fill(lambdaM,ww);
	    hdist_pmHpFTAll->Fill(lambdaD,ww);

	    hZvertFTAll->Fill(vertexL.Z());
		    
	    //cout<<"::: "<<pimH.size()<<" "<<pFT.size()<<endl;
	    //cout<<"::: "<<lambdaM<<endl;

	    /*
	    //recognize proton from FT origin------------
	    //int tn=pFT[j]->getTrack();
	    int tn;
	    int kineID=-1;
	    int kineparentID=-1;
	    for(int p=0;p<knum;p++)
	    {
	    kine=HCategoryManager::getObject(kine, fCatGeantKine,p);
	    int kineT=kine->getTrack();
	    if(kineT==tn)
	    {
	    kineID=kine->getID();
	    kineparentID=getMotherIndex(kine);
	    break;
	    }
	    }
	    if(kineID==14 && kineparentID==18)
	    hFDdistanceLambda->Fill(lambdaD);

	    if(kineID==14 && kineparentID!=18)
	    hFDdistanceBg->Fill(lambdaD);
	    //---------------------------------------
	    */


	    if(pFT[j]->getGeantParentPID()==18 && pimH[k]->getGeantParentPID()==18)//proton and pion from L(1115)
	      {
		hZvertFTRL->Fill(vertexL.Z());
		hLRmassFT->Fill(lambdaM,ww);
		hLRdistFT->Fill(lambdaD,ww);
	      }
	    flagHFTL1=0;
	    flagHFTL2=0;
	    flagHFTL3=0;
	    flagHFTL4=0;
    
	    if(lambdaD<min_dist_l)
	      {
		flagHFTL1=1;
		hinvM_pmHpFTDist->Fill(lambdaM,ww);
		//h2Lvertex->Fill(vertexL.Z(),TMath::Sqrt(vertexL.X()*vertexL.X()+vertexL.Y()*vertexL.Y()));
		if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
		  hLRmassFTDist->Fill(lambdaM,ww);
	      }

	    if (lambdaD<min_dist_l && vertexL.Z()>0.)
	      {
		flagHFTL2=1;
		hinvM_pmHpFTDistZ->Fill(lambdaM,ww);
		if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)hLRmassFTDistZ->Fill(lambdaM,ww);

	      }
	    
	    ver_LTg=tool.calcVertexAnalytical(base_Tg,dir_Tg,vertexL,dirL);
	    double distLamZ=tool.calculateMinimumDistance(base_Tg, dir_Tg,vertexL,dirL);

	    if (lambdaM>1105 && lambdaM<1125 && lambdaD<min_dist_l)
	      {
		flagHFTL3=1;
		hinvM_pmHpFTDistL->Fill(lambdaM,ww);
		if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18) hLRmassFTDistL->Fill(lambdaM,ww);
		      
		hZvertLamFT->Fill(vertexL.Z());
		hZvertLamTgFT->Fill(ver_LTg.Z());
		hdistTgLamFT->Fill(distLamZ);
		      		      
	      }
	    if (lambdaM>1105 && lambdaM<1125 && lambdaD<min_dist_l && vertexL.Z()>0.)
	      {
		      
		flagHFTL4=1;
		      
		hinvM_pmHpFTDistZL->Fill(lambdaM,ww);
		if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18) hLRmassFTDistZL->Fill(lambdaM,ww);
		      
	      }

	    if (ep.size() || em.size()){
	      for (int s=0;s<ep.size();s++){
		for (int ss=0;ss<ep.size();ss++){
		  TLorentzVector lvLambda=*pFT[j]+*pimH[k];
		  TLorentzVector lvDiLepton=*ep[s]+*ep[ss];
		  double oa = tool.getOpeningAngle(ep[s],ep[ss]);
		  double mass_1520=(lvLambda+lvDiLepton).M();
		  int flagDil2=1;
		    

		  if(oa>min_angle)
		    {
		      hinvMass_epepOA->Fill(lvDiLepton.M(),ww);  
		      hDLmassAll_epep->Fill(lvDiLepton.M(),ww);
		      hL1520massAll_epep->Fill(mass_1520,ww);

		      if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
			{
			  hDLmassAllRL_epep->Fill(lvDiLepton.M(),ww);
			  hL1520massAllRL_epep->Fill(mass_1520,ww);
		      
			  if(flagDil2)
			    {
			      hDLmassAllRL_L_epep->Fill(lvDiLepton.M(),ww);
			      hL1520massAllRL_L_epep->Fill(mass_1520,ww);
			
			    }
			}

		      if(flagHFTL1)
			{
			  hDLmassDist_epep->Fill(lvDiLepton.M(),ww);
			  hL1520massDist_epep->Fill(mass_1520,ww);

			  if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
			    {
			      hDLmassDistRL_epep->Fill(lvDiLepton.M(),ww);
			      hL1520massDistRL_epep->Fill(mass_1520,ww);

			      if(flagDil2)
				{
				  hDLmassDistRL_L_epep->Fill(lvDiLepton.M(),ww);
				  hL1520massDistRL_L_epep->Fill(mass_1520,ww);
				}
			    }
			}

		      if(flagHFTL2)
			{
			  hDLmassDistZ_epep->Fill(lvDiLepton.M(),ww);
			  hL1520massDistZ_epep->Fill(mass_1520,ww);

			  if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
			    {
			      hDLmassDistZRL_epep->Fill(lvDiLepton.M(),ww);
			      hL1520massDistZRL_epep->Fill(mass_1520,ww);

			      if(flagDil2)
				{
				  hDLmassDistZRL_L_epep->Fill(lvDiLepton.M(),ww);
				  hL1520massDistZRL_L_epep->Fill(mass_1520,ww);
				}
			    }
			}
		      if(flagHFTL3)
			{
			  hDLmassDistL_epep->Fill(lvDiLepton.M(),ww);
			  hL1520massDistL_epep->Fill(mass_1520,ww);
			  if(lvDiLepton.M()>140. && lvDiLepton.M()<420.)hL1520massDistLpi0_epep->Fill(mass_1520,ww);

			  if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
			    {
			      hDLmassDistLRL_epep->Fill(lvDiLepton.M(),ww);
			      hL1520massDistLRL_epep->Fill(mass_1520,ww);

			      if(flagDil2)
				{
				  hDLmassDistLRL_L_epep->Fill(lvDiLepton.M(),ww);
				  hL1520massDistLRL_L_epep->Fill(mass_1520,ww);
				}
			    }
			}
    		    
		    }//oa

		  //***************************************

		    
		  if(flagHFTL4 && oa>min_angle)
		    {
		      hDLmassDistZL_epep->Fill(lvDiLepton.M(),ww);
		      
		      hL1520mass_epep->Fill(mass_1520,ww);
		      if(mass_1520>1450 && mass_1520<1550)
			hinvMass_epep->Fill(lvDiLepton.M(),ww);

		      if(flagDil2)
			{
			  hL1520mass_epep_L->Fill(mass_1520,ww);
			  if(mass_1520>1450 && mass_1520<1550)hinvMass_epep_L->Fill(lvDiLepton.M(),ww);
			  if(lvDiLepton.M()>140. && lvDiLepton.M()<420.)
			    hL1520massDistZLpi0_epep->Fill(mass_1520,ww);
			}

		    }
		    
		}
	      }
	

	      for (int s=0;s<em.size();s++)
		{
		  for (int ss=0;ss<em.size();ss++)
		    {
		      TLorentzVector lvLambda=*pFT[j]+*pimH[k];
		      TLorentzVector lvDiLepton=*em[s]+*em[ss];
		      double oa = tool.getOpeningAngle(em[s],em[ss]);
		      double mass_1520=(lvLambda+lvDiLepton).M();
		    
		      int flagDil1=1;
		   
		      //if(em[s]->getGeantParentTrackNum()==0 && em[ss]->getGeantParentTrackNum()==0 && em[s]->getGeantGrandParentPID()==-1 && em[ss]->getGeantGrandParentPID()==-1) flagDil1=1;
		      //if(em[s]->getGeantParentTrackNum()== em[ss]->getGeantParentTrackNum()) flagDil1=1;

		      //***************************************

		      if(oa>min_angle)
			{
			  hinvMass_ememOA->Fill(lvDiLepton.M(),ww);  
			  hDLmassAll_emem->Fill(lvDiLepton.M(),ww);
			  hL1520massAll_emem->Fill(mass_1520,ww);

			  if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
			    {
			      hDLmassAllRL_emem->Fill(lvDiLepton.M(),ww);
			      hL1520massAllRL_emem->Fill(mass_1520,ww);
		      
			      if(flagDil1)
				{
				  hDLmassAllRL_L_emem->Fill(lvDiLepton.M(),ww);
				  hL1520massAllRL_L_emem->Fill(mass_1520,ww);
				}
			    }

			  if(flagHFTL1)
			    {
			      hDLmassDist_emem->Fill(lvDiLepton.M(),ww);
			      hL1520massDist_emem->Fill(mass_1520,ww);

			      if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
				{
				  hDLmassDistRL_emem->Fill(lvDiLepton.M(),ww);
				  hL1520massDistRL_emem->Fill(mass_1520,ww);

				  if(flagDil1)
				    {
				      hDLmassDistRL_L_emem->Fill(lvDiLepton.M(),ww);
				      hL1520massDistRL_L_emem->Fill(mass_1520,ww);
				    }
				}
			    }

			  if(flagHFTL2)
			    {
			      hDLmassDistZ_emem->Fill(lvDiLepton.M(),ww);
			      hL1520massDistZ_emem->Fill(mass_1520,ww);
		      
			      if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
				{
				  hDLmassDistZRL_emem->Fill(lvDiLepton.M(),ww);
				  hL1520massDistZRL_emem->Fill(mass_1520,ww);

				  if(flagDil1)
				    {
				      hDLmassDistZRL_L_emem->Fill(lvDiLepton.M(),ww);
				      hL1520massDistZRL_L_emem->Fill(mass_1520,ww);
			  
				    }
				}
			    }

			  if(flagHFTL3)
			    {
			      hDLmassDistL_emem->Fill(lvDiLepton.M(),ww);
			      hL1520massDistL_emem->Fill(mass_1520,ww);
			      if(lvDiLepton.M()>140. && lvDiLepton.M()<420.)hL1520massDistLpi0_emem->Fill(mass_1520,ww);

			      if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18){

				hDLmassDistLRL_emem->Fill(lvDiLepton.M(),ww);
				hL1520massDistLRL_emem->Fill(mass_1520,ww);

				if(flagDil1){
				  hDLmassDistLRL_L_emem->Fill(lvDiLepton.M(),ww);
				  hL1520massDistLRL_L_emem->Fill(mass_1520,ww);
				}
			
			      }
			    }
			}//oa
		      //***************************************

		      if(flagHFTL4 && oa>min_angle)
			{
			  hDLmassDistZL_emem->Fill(lvDiLepton.M(),ww);
			  hL1520mass_emem->Fill(mass_1520,ww);		      
			  if(mass_1520>1450 && mass_1520<1550)
			    hinvMass_emem->Fill(lvDiLepton.M(),ww);
		 
			  if(flagDil1)
			    {

			      hL1520mass_emem_L->Fill(mass_1520,ww);		      
			      if(mass_1520>1450 && mass_1520<1550)hinvMass_emem_L->Fill(lvDiLepton.M(),ww);

			      if(lvDiLepton.M()>140. && lvDiLepton.M()<420.)
				hL1520massDistZLpi0_emem->Fill(mass_1520,ww);
			    }
			}
		    }
		}
	    }
	
	    for (int s=0;s<ep.size();s++)
	      {
		for (int ss=0;ss<em.size();ss++)
		  {	     
		    int flagDil=0;

		    //if(ep[s]->getGeantParentTrackNum()==0 && em[ss]->getGeantParentTrackNum()==0 && ep[s]->getGeantGrandParentPID()==-1 && em[ss]->getGeantGrandParentPID()==-1) flagDil=1;
		    if(ep[s]->getGeantParentTrackNum()==em[ss]->getGeantParentTrackNum() /*&& em[s]->getGeantParentPID()==7*/)
		      flagDil=1;
			      
		    TLorentzVector lvLambda=*pFT[j]+*pimH[k];
		    TLorentzVector lvDiLepton=*ep[s]+*em[ss];

		    vertexDL=trackVertex(ep[s],em[ss]);
		    dirDL.setXYZ((*ep[s]+*em[ss]).X(),(*ep[s]+*em[ss]).Y(),(*ep[s]+*em[ss]).Z());

			    
		    double mass_1520=(lvLambda+lvDiLepton).M();
		    //min dist between dilepton and Lam1115
		    double distance_1520=tool.calculateMinimumDistance(vertexL,dirL,vertexDL,dirDL);
		    double invMepem= lvDiLepton.M(); 
		    double oa = tool.getOpeningAngle(ep[s],em[ss]);
		    double dilTrDist=trackDistance(ep[s],em[ss]);
		    int oaFlag=0;
		    
		    if(oa>min_angle)
		      oaFlag=1;
		    
		    if(oaFlag)
		      {
			hinvMass_epemOA->Fill(invMepem,ww);  
			hDLdistanceFT->Fill(distance_1520);
			hDLmassFTAll->Fill(invMepem,ww);
			hL1520massFTAll->Fill(mass_1520,ww);

			if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
			  {
			    hDLmassFTAllRL->Fill(invMepem,ww);
			    hL1520massFTAllRL->Fill(mass_1520,ww);

			    if(flagDil)
			      {
				hDLmassAllRL_L->Fill(invMepem,ww);
				hL1520massAllRL_L->Fill(mass_1520,ww);
			      }
			  }

			if(flagHFTL1)
			  {
			    hDLmassFTDist->Fill(invMepem,ww);
			    hL1520massFTDist->Fill(mass_1520,ww);

			    if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
			      {
				hDLmassFTDistRL->Fill(invMepem,ww);
				hL1520massFTDistRL->Fill(mass_1520,ww);

				if(flagDil)
				  {
				    hDLmassDistRL_L->Fill(invMepem,ww);
				    hL1520massDistRL_L->Fill(mass_1520,ww);

				  }
			      }
			  }

			if(flagHFTL2)
			  {
			    hDLmassFTDistZ->Fill(invMepem,ww);
			    hL1520massFTDistZ->Fill(mass_1520,ww);

			    if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
			      {
				hDLmassFTDistZRL->Fill(invMepem,ww);
				hL1520massFTDistZRL->Fill(mass_1520,ww);

				if(flagDil)
				  {
				    hDLmassDistZRL_L->Fill(invMepem,ww);
				    hL1520massDistZRL_L->Fill(mass_1520,ww);
			  
				  }
			      }
			  }

			if(flagHFTL3){

			  hDLmassFTDistL->Fill(invMepem,ww);
			  hL1520massFTDistL->Fill(mass_1520,ww);
			  if(flagDil)
			    {
			      hDLmassDistL_L->Fill(invMepem,ww);
			    }
			  if(invMepem>140. && invMepem<420.)
			    hL1520massFTDistLpi0->Fill(mass_1520,ww);

			  if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
			    {
			      hDLmassFTDistLRL->Fill(invMepem,ww);
			      hL1520massFTDistLRL->Fill(mass_1520,ww);

			      if(flagDil)
				{
				  hDLmassDistLRL_L->Fill(invMepem,ww);
				  hL1520massDistLRL_L->Fill(mass_1520,ww);
				}
			    }
			}
		    
			if(flagHFTL4)
			  {
			    hDLmassFTDistZL->Fill(invMepem,ww);
			    hL1520massFTDistZL->Fill(mass_1520,ww);
			    if(flagDil)
			      {
				hDLmassDistZL_L->Fill(invMepem,ww);
			      }
			    if(invMepem>140. && invMepem<420.)
			      {
				hL1520massDistZLpi0->Fill(mass_1520,ww);
				if(flagDil)hL1520massDistZLpi0_L->Fill(mass_1520,ww);
			      }
		      
			    if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
			      {
				hDLmassFTDistZLRL->Fill(invMepem,ww);
				hL1520massFTDistZLRL->Fill(mass_1520,ww);
				if(invMepem>140. && invMepem<420.)hL1520massFTDistZLRLpi0->Fill(mass_1520,ww);
				if(flagDil)
				  {
				    hDLmassDistZLRL_L->Fill(invMepem,ww);
				    hL1520massDistZLRL_L->Fill(mass_1520,ww);
				    if(invMepem>140. && invMepem<420.)hL1520massDistZLRLpi0_L->Fill(mass_1520,ww);
				  }
			      }
		   
			
			    if(mass_1520>1450 && mass_1520<1550)
			      {
				hL1520massFTFinal->Fill(mass_1520,ww);
				hDLmassFTFinal->Fill(invMepem,ww);
				if(invMepem>140. && invMepem<420.)hL1520massFTFinalpi0->Fill(mass_1520,ww);

				if(pimH[k]->getGeantParentPID()==18 && pFT[j]->getGeantParentPID()==18)
				  {
				    hL1520massFTFinalRL->Fill(mass_1520,ww);
				    hDLmassFTFinalRL->Fill(invMepem,ww);
				    if(invMepem>140. && invMepem<420.)hL1520massFTFinalRLpi0->Fill(mass_1520,ww);

				    if(flagDil)
				      {
					hDLmassFinalRL_L->Fill(invMepem,ww);
					hL1520massFinalRL_L->Fill(mass_1520,ww);
					if(invMepem>140. && invMepem<420.){hL1520massFinalRLpi0_L->Fill(mass_1520,ww);
					  //cout<<"xxxxx "<<mass_1520<<" "<<ww<<endl;
			
					}
				      }
				  }
			      }

			    //h2L1520vertex->Fill(vertexL1520.Z(),TMath::Sqrt(vertexL1520.X()*vertexL1520.X()+vertexL1520.Y()*vertexL1520.Y()));
			  
			    TLorentzVector lvLam1520=lvLambda+lvDiLepton;
			    //HParticleCandSim *tr;
			    HGeomVector base_Tg, dir_Tg, ver_L1520Tg, baseL1520, vertexL1520;
			    HParticleTool p_tool, tool,tool1;


			    //tr=ep[s]+em[ss]+pimH[j]+pH[k];

			    double distance_1520=tool1.calculateMinimumDistance(vertexL,dirL,vertexDL,dirDL);
		
			    vertexL1520=tool.calcVertexAnalytical(vertexL,dirL,vertexDL,dirDL);
			    //dirL1520.setXYZ((*ep[s]+*em[ss]+*pimH[j]+*pH[k]).X(),(*ep[s]+*em[ss]+*pimH[j]+*pH[k]).Y(),(*ep[s]+*em[ss]+*pimH[j]+*pH[k]).Z());
		
			    //dirL1520.setXYZ(lvLam1520.X(),lvLam1520.Y(),lvLam1520.Z());
			    //p_tool.calcSegVector(tr->getZ(),tr->getR(),TMath::DegToRad()*tr->getPhi(),TMath::DegToRad()*tr->getTheta(),baseL1520,dirL1520_1);
			  	  
			    //ver_L1520Tg = p_tool.calcVertexAnalytical(base_Tg,dir_Tg,baseL1520,dirL1520);

			    ver_L1520Tg=tool.calcVertexAnalytical(base_Tg,dir_Tg,vertexL1520,dirL1520);
		     
			    //if(vertexL1520.Z()>0){
			    //cout<<"----------------->FT>>1 "<<s<<" "<<ss<<" "<<vertexL1520.Z()<<endl;
			    //cout<<"----------------->>2 "<<s<<" "<<ss<<" "<<ver_L1520Tg.Z()<<endl;
			  		  
			    double distLam1520Z=tool.calculateMinimumDistance(base_Tg, dir_Tg,vertexL1520,dirL1520);
			    //cout<<"distLam1520Z: "<<distLam1520Z<<endl;

			    hZvertLam1520FT->Fill(vertexL1520.Z());
			    hZvertLam1520TgFT->Fill(ver_L1520Tg.Z());
			    hdistTgLam1520FT->Fill(distLam1520Z);
			  
			    //}
			    //if(invMdilLam>1400 && invMdilLam<1700)hDLmassDistOALcut->Fill(invMepem);	    
			  }
		      }//oa
			    	      
		  }//em
			    
	      }//ep
			  
			  

	  }
	   	    
	    
	}
      }
	
      //**************************************************************
      //kine analysis*****************************************
      for(int p=0;p<knum;p++)
	{
	  kine=HCategoryManager::getObject(kine, fCatGeantKine,p);
	  int kineID=kine->getID();
	  int mech=kine->getMechanism();
	  int kineparentID=getMotherIndex(kine);
	  HGeomVector lambdaVertex;

	  //hPRmotherindex->Fill(kineparentID);

	  //if(kineparentID==18)hPRLambdaDoughter->Fill(kineID);

	  if(mech==0 && (kineID==2 || kineID==3))
	    {
	      //h2IIleptonsfromPV->Fill(kine->getTotalMomentum(),kine->getThetaDeg());

	      hEleptons4Pi->Fill(kine->getThetaDeg());
	      h2Elepton4Pi->Fill(kine->getPhiDeg(),kine->getThetaDeg());
		
	    }
	  /*
	    int fld=0;
	    if(mech==0 && kineID==2)
	    {
		
	    hEleptons4Pi->Fill(kine->getThetaDeg());
	    h2Elepton4Pi->Fill(kine->getPhiDeg(),kine->getThetaDeg());
	    //fld=1;
	    }
	    
	    if(mech==0 && kineID==3)
	    {

	    //if(fld)
	    h2Elepton4Pi->Fill(kine->getPhiDeg(),kine->getThetaDeg());
	    
	    }
	  */
	  if(kineID==14 && kineparentID==18)//proton from lambda
	    {
	      //h2IIprotons->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
	      //hEprotons4Pi->Fill(kine->getThetaDeg());
	      //h2Eproton4Pi->Fill(kine->getPhiDeg(),kine->getThetaDeg());
	      hEFLprotons4Pi->Fill(kine->getThetaDeg());
	    }
	  if(kineID==14 && mech==0)
	    {//proton from primary vertex
	      //if(kineID==14){
	      //h2IIprotons->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
	      hEprotons4Pi->Fill(kine->getThetaDeg());
	      h2Eproton4Pi->Fill(kine->getPhiDeg(),kine->getThetaDeg());
	    }
	  if(kineID==9 && kineparentID==18)//Pi- from Lambda
	    {
	      kine->getVertex(lambdaVertex);
	      //h2IIpions->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
	      hEFLpions4Pi->Fill(kine->getThetaDeg());
		
	      //hEIZpionSim->Fill(lambdaVertex.getZ());
	    }

	  if(kineID==9 && mech==0){//pim from primary vertex
	    //if(kineID==9){
	    hEpions4Pi->Fill(kine->getThetaDeg());
	    h2Epion4Pi->Fill(kine->getPhiDeg(),kine->getThetaDeg());
	  }
	    
	    
	  //if(kineID==14 && kine->getThetaDeg()<6.5)h2FDsimProtons->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
	}
      //******************************************************end kine
     

	
      //  WriteHistos();

     
    }

  //***************************
  hinvM_pmHpHAll->Write();
  hinvM_pmHpFTAll->Write();
    
  hinvM_pmHpHDist->Write();
  hinvM_pmHpFTDist->Write();
 
  hinvM_pmHpHDistZ->Write();
  hinvM_pmHpFTDistZ->Write();
 
  hinvM_pmHpHDistL->Write();
  hinvM_pmHpFTDistL->Write();

  hinvM_pmHpHDistZL->Write();
  hinvM_pmHpFTDistZL->Write();

  hLRmass->Write();
  hLRmassFT->Write();
   
  hLRmassDist->Write();
  hLRmassFTDist->Write();

  hLRmassDistZ->Write();
  hLRmassFTDistZ->Write();

  hLRmassDistL->Write();
  hLRmassFTDistL->Write();
    
  hLRmassDistZL->Write();
  hLRmassFTDistZL->Write();
   
      
  hDLmassAll->Write();
  hDLmassFTAll->Write();
  hL1520massAll->Write();
  hL1520massFTAll->Write();

  hDLmassAllRL->Write();
  hDLmassFTAllRL->Write();
  hL1520massAllRL->Write();
  hL1520massFTAllRL->Write();
    
  hDLmassDist->Write();
  hDLmassFTDist->Write();
  hL1520massDist->Write();
  hL1520massFTDist->Write();

  hDLmassDistRL->Write();
  hDLmassFTDistRL->Write();
  hL1520massDistRL->Write();
  hL1520massFTDistRL->Write();
     
  hDLmassDistZ->Write();
  hDLmassFTDistZ->Write();
  hL1520massDistZ->Write();
  hL1520massFTDistZ->Write();

  hDLmassDistZRL->Write();
  hDLmassFTDistZRL->Write();
  hL1520massDistZRL->Write();
  hL1520massFTDistZRL->Write();
     
  hDLmassDistL->Write();
  hDLmassFTDistL->Write();
  hL1520massDistL->Write();
  hL1520massDistLpi0->Write();

  hL1520massFTDistL->Write();
  hL1520massFTDistLpi0->Write();
    
  hDLmassDistLRL->Write();
  hDLmassFTDistLRL->Write();
  hL1520massDistLRL->Write();
  hL1520massFTDistLRL->Write();
     
  hDLmassDistZL->Write();
  hDLmassFTDistZL->Write();
  hL1520massDistZL->Write();
  hL1520massFTDistZL->Write();

  hDLmassDistZLRL->Write();
  hDLmassFTDistZLRL->Write();
  hL1520massDistZLRL->Write();
  hL1520massFTDistZLRL->Write();
     
  hDLmassFinal->Write();
  hDLmassFTFinal->Write();
  hL1520massFinal->Write();
  hL1520massFTFinal->Write();

  hDLmassFinalRL->Write();
  hDLmassFTFinalRL->Write();
  hL1520massFinalRL->Write();
  hL1520massFTFinalRL->Write();


  hL1520massDistZLpi0->Write();
  //hL1520massFTDistZLpi0->Write();
    
  hL1520massFinalpi0->Write();
  hL1520massFTFinalpi0->Write();
  hL1520massFinalRLpi0->Write();
  hL1520massFTFinalRLpi0->Write();


  hL1520massDistZLpi0_L->Write();

    
  //***************
  hdist_pmHpHAll->Write();
  hdist_pmHpFTAll->Write();
   
  hLRdist->Write();
  hLRdistFT->Write();
   
  hDLdistance->Write();
  hDLdistanceFT->Write();

  h2L1520vertex->Write(); 
  h2LHvertex ->Write();
     
   
   
  
  //***********************

  
  hZvertLam1520HH->Write();
  hZvertLam1520FT->Write();
  
  hdistTgLam1520FT->Write();
  hdistTgLam1520HH->Write();

  //hZvertLam1520HH ->Write();
  hdistTgLam1520HH ->Write();
  hZvertLam1520TgHH ->Write();
 
  hZvertHHAll ->Write();
  hZvertLamHH ->Write();
  hZvertHHRL->Write();
  
  hZvertFTAll ->Write();
  hZvertLamFT ->Write();
  hZvertFTRL->Write();

  hdistTgLamHH ->Write();
  hZvertLamTgHH ->Write();
  hdistTgLamFT ->Write();
  hZvertLamTgFT ->Write();
   

  hZvertLam1520FT ->Write();
  hdistTgLam1520FT ->Write();
  hZvertLam1520TgFT ->Write();
  

  hinvMass_HHemem->Write();  
  hinvMass_HHepep->Write();  
  hinvMass_HFTemem->Write();  
  hinvMass_HFTepep->Write();  

  hinvMass_epemOA->Write();  
  hinvMass_ememOA->Write();  
  hinvMass_epepOA->Write();  

 
  hL1520mass_HHemem->Write();
  hL1520mass_HHepep->Write();
  hL1520mass_HFTemem->Write();
  hL1520mass_HFTepep->Write();

  hL1520mass_epep->Write();
  hL1520mass_emem->Write();
  
 
  hDLmassAllRL_L->Write();
  hL1520massAllRL_L->Write();
  hDLmassDistRL_L->Write();
  hL1520massDistRL_L->Write();
  hDLmassDistZRL_L->Write();
  hL1520massDistZRL_L->Write();
  hDLmassDistLRL_L->Write();
  hL1520massDistLRL_L->Write();
  hDLmassDistZLRL_L->Write();
  hL1520massDistZLRL_L->Write();
  hL1520massDistZLRLpi0_L->Write();
  hDLmassFinalRL_L->Write();
  hL1520massFinalRL_L->Write();
  hL1520massFinalRLpi0_L->Write();
  hinvMass_epep_L->Write();
  hinvMass_emem_L->Write();

  hDLmassDistZL_L->Write();
  hDLmassDistL_L->Write();
  
  //**************************************************************************** 
 
  hDLmassAll_epep->Write();
  hL1520massAll_epep->Write();
  hDLmassAllRL_epep->Write();
  hL1520massAllRL_epep->Write();
  hDLmassAllRL_L_epep->Write();
  hL1520massAllRL_L_epep->Write();
  hDLmassDist_epep->Write();
  hL1520massDist_epep->Write();
  hDLmassDistRL_epep->Write();
  hL1520massDistRL_epep->Write();
  hDLmassDistRL_L_epep->Write();
  hL1520massDistRL_L_epep->Write();
  hDLmassDistZ_epep->Write();
  hL1520massDistZ_epep->Write();
  hDLmassDistZRL_epep->Write();
  hL1520massDistZRL_epep->Write();
  hDLmassDistZRL_L_epep->Write();
  hL1520massDistZRL_L_epep->Write();
  hDLmassDistL_epep->Write();
  hL1520massDistL_epep->Write();
  hL1520massDistLpi0_epep->Write();
  hDLmassDistLRL_epep->Write();
  hL1520massDistLRL_epep->Write();
  hDLmassDistLRL_L_epep->Write();
  hL1520massDistLRL_L_epep->Write();
  hDLmassDistZL_epep->Write();
  hL1520massDistZL_epep->Write();
  hDLmassDistZLRL_epep->Write();
  hL1520massDistZLRL_epep->Write();
  hDLmassDistZLRL_L_epep->Write();
  hL1520massDistZLRL_L_epep->Write();
  hDLmassFinal_epep->Write();
  hL1520massFinal_epep->Write();
  hDLmassFinalRL_epep->Write();
  hL1520massFinalRL_epep->Write();
  hDLmassFinalRL_L_epep->Write();
  hL1520massFinalRL_L_epep->Write();
  hL1520massFinalRLpi0_epep->Write();
  hL1520massFinalRLpi0_L_epep->Write();
  hL1520massFinalpi0_epep->Write();
  hL1520massDistZLpi0_epep->Write();
    
  
  //********************************
  hDLmassAll_emem->Write();
  hL1520massAll_emem->Write();
  hDLmassAllRL_emem->Write();
  hL1520massAllRL_emem->Write();
  hDLmassAllRL_L_emem->Write();
  hL1520massAllRL_L_emem->Write();
  hDLmassDist_emem->Write();
  hL1520massDist_emem->Write();
  hDLmassDistRL_emem->Write();
  hL1520massDistRL_emem->Write();
  hDLmassDistRL_L_emem->Write();
  hL1520massDistRL_L_emem->Write();
  hDLmassDistZ_emem->Write();
  hL1520massDistZ_emem->Write();
  hDLmassDistZRL_emem->Write();
  hL1520massDistZRL_emem->Write();
  hDLmassDistZRL_L_emem->Write();
  hL1520massDistZRL_L_emem->Write();
  hDLmassDistL_emem->Write();
  hL1520massDistL_emem->Write();
  hL1520massDistLpi0_emem->Write();
  hDLmassDistLRL_emem->Write();
  hL1520massDistLRL_emem->Write();
  hDLmassDistLRL_L_emem->Write();
  hL1520massDistLRL_L_emem->Write();
  hDLmassDistZL_emem->Write();
  hL1520massDistZL_emem->Write();
  hDLmassDistZLRL_emem->Write();
  hL1520massDistZLRL_emem->Write();
  hDLmassDistZLRL_L_emem->Write();
  hL1520massDistZLRL_L_emem->Write();
  hDLmassFinal_emem->Write();
  hL1520massFinal_emem->Write();
  hDLmassFinalRL_emem->Write();
  hL1520massFinalRL_emem->Write();
  hDLmassFinalRL_L_emem->Write();
  hL1520massFinalRL_L_emem->Write();
  hL1520massFinalRLpi0_emem->Write();
  hL1520massFinalRLpi0_L_emem->Write();
  hL1520massFinalpi0_emem->Write();
  hL1520massDistZLpi0_emem->Write();

  proton_identyf_ideal->Write();
  proton_identyf_real->Write();
  pion_identyf_ideal->Write();
  pion_identyf_real->Write();
  lepton_identyf_ideal->Write();
  lepton_identyf_real->Write();

  lepton_identyf_ideal_dphi->Write();
  lepton_identyf_real_dphi->Write();
  lepton_identyf_ideal_dtheta->Write();
  lepton_identyf_real_dtheta->Write();
  /*
    TCanvas* c_identyf=new TCanvas("c_identyf","c_identyf");
    c_identyf->Divide(2,3);
    c_identyf->cd(1);
    proton_identyf_ideal->Draw();
    c_identyf->cd(2);
    proton_identyf_real->Draw();
    c_identyf->cd(3);
    pion_identyf_ideal->Draw();
    c_identyf->cd(4);
    pion_identyf_real->Draw();
    c_identyf->cd(5);
    lepton_identyf_real->Draw();
    c_identyf->cd(6);
    lepton_identyf_ideal->Draw();

    c_identyf->Write();
    //********************************
    /*
    cEff->Divide(4,3);
    cEff->cd(1);
    hEprotons4Pi->Draw();
    hEprotonsdet->SetLineColor(kRed);
    hEprotonsdet->Draw("SAME");
    //cEff->cd(2);
    //hEprotonsdet->Draw();
    cEff->cd(2);
    hEprotonsEff->Divide(hEprotonsdet,hEprotons4Pi);
    hEprotonsEff->Draw();
    //cEff->cd(4);
    //h2EprotonEff->Divide(h2Eprotondet,h2Eproton4Pi);
    //h2EprotonEff->Draw("COLZ");
    cEff->cd(3);
    //h2EprotonEff->Divide(h2Eprotondet,h2Eproton4Pi);
    h2Eprotondet->Draw("COLZ");
    cEff->cd(4);
    //h2EprotonEff->Divide(h2Eprotondet,h2Eproton4Pi);
    h2Eproton4Pi->Draw("COLZ");


    cEff->cd(5);
    hEpions4Pi->Draw();
    hEpionsdet->SetLineColor(kRed);
    hEpionsdet->Draw("SAME");
    cEff->cd(6);
    //hEpionsdet->Draw();
    hEpionsEff->Divide(hEpionsdet,hEpions4Pi);
    hEpionsEff->Draw();
    cEff->cd(7);
    h2Epiondet->Draw("colz");
    cEff->cd(8);
    h2Epion4Pi->Draw("colz");
    //h2EpionEff->Divide(h2Epiondet,h2Epion4Pi);
    //h2EpionEff->Draw("COLZ");

    cEff->cd(9);
    hEleptons4Pi->Draw();
    hEleptonsdet->SetLineColor(kRed);
    hEleptonsdet->Draw("SAME");
    cEff->cd(10);
    hEleptonsEff->Divide(hEleptonsdet,hEleptons4Pi);
    hEleptonsEff->Draw();
    cEff->cd(11);
    h2Eleptondet->Draw("colz");
    cEff->cd(12);
    h2Elepton4Pi->Draw("colz");

    cEff->Write();
  */
  //*****************************
  /*
    cEffFromLambda->Divide(3,2);
    cEffFromLambda->cd(1);
    hEFLprotons4Pi->Draw();
    cEffFromLambda->cd(2);
    hEFLprotonsdet->Draw();
    cEffFromLambda->cd(3);
    hEFLprotonsEff->Divide(hEFLprotonsdet,hEFLprotons4Pi);
    hEFLprotonsEff->Draw();

    cEffFromLambda->cd(4);
    hEFLpions4Pi->Draw();
    cEffFromLambda->cd(5);
    hEFLpionsdet->Draw();
    cEffFromLambda->cd(6);
    hEFLpionsEff->Divide(hEFLpionsdet,hEFLpions4Pi);
    hEFLpionsEff->Draw();

    cEffFromLambda->Write();
  */

  //hinvMepemHH->Write();
  //hinvMepemFT->Write();
      
}



