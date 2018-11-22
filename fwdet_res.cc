
#include "fwdet_res.h"
#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hfwdetcand.h"
#include "hfwdetcandsim.h"
#include <TCanvas.h>
#include <TStyle.h>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;



HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);
HGeomVector calcPrimVertex_Track_Mother(const std::vector<HParticleCand *>cands, const HGeomVector & beamVector, const HGeomVector & DecayVertex, const HGeomVector & dirMother, int trackA_num, int trackB_num);

Int_t fwdet_tests(HLoop * loop, const AnaParameters & anapars)
{
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
    if (!loop->setInput(""))
    {                                                    // reading file structure
        std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    //////////////////////////////////////////////////////////////////////////////
    //      Fast tree builder for creating of ntuples                            //
    //////////////////////////////////////////////////////////////////////////////

    loop->printCategories();    // print all categories found in input + status
    //     loop->printChain();            // print all files in the chain
    //     loop->Print();

    //     HEventHeader * header = loop->getEventHeader();
    HCategory * fCatGeantKine = nullptr;
    fCatGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");

    if (!fCatGeantKine)
    {
        cout << "No catGeantKine!" << endl;
        exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    // HCategory * fFwDetStrawCal = nullptr;
    // fFwDetStrawCal = HCategoryManager::getCategory(catFwDetStrawCal, kTRUE, "catFwDetStrawCalSim");

    // if (!fFwDetStrawCal)
    // {
    //  cout << "No catFwDetStrawCal!" << endl;
    //  exit(EXIT_FAILURE);  // do you want a brute force exit ?
    // }

    HCategory * fCatVectorCand = nullptr;
    fCatVectorCand = HCategoryManager::getCategory(catVectorCand, kTRUE, "catVectorCand");

    // HCategory * fCatVectorCandSim = nullptr;
    // fCatVectorCandSim = HCategoryManager::getCategory(catVectorCand, kTRUE, "catVectorCand");

    if (!fCatVectorCand)
    {
        cout << "No catVectorCand!" << endl;
	//exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    //
    Int_t entries = loop->getEntries();
    //     //setting numbers of events regarding the input number of events by the user
    if (anapars.events < entries and anapars.events >= 0 ) entries = anapars.events;

    //     // specify output file
    TFile * output_file = TFile::Open(anapars.outfile, "RECREATE");
    output_file->cd();
    //
    cout << "NEW ROOT TREE " << endl;
    //
    //crete histograms
    
    //event loop *************************************************
    //*********************************************************
    for (Int_t i = 0; i < entries; i++)                   
    {
      if(i % 1409 ==0)
	cout<<"Event No. "<<i<<endl;
        /*Int_t nbytes =*/  loop->nextEvent(i);         // get next event. categories will be cleared before
        
        int geant_hKine_cnt = fCatGeantKine->getEntries();

        HGeantKine * gKine = nullptr;
	HVectorCandSim* fwdetstrawvec = nullptr;
	//HVectorCand* fwdetstrawvecSim= new HVectorCand;
	
        // kine
	for (int j = 0; j < geant_hKine_cnt; j++)
        {
            gKine = HCategoryManager::getObject(gKine, fCatGeantKine, j);
	   
	 }

	//resolution reconstruction
    } // end eventloop
	//***********************************************************************************
	
    //Drawing

    	
    output_file->Close();
    cout << "writing root tree done" << endl;

    timer.Stop();
    timer.Print();

    return 0;
}
