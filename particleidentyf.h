#ifndef particleidentyf_hpp
#define particleidentyf_hpp
#include "hades.h"
#include "hloop.h"
#include "htaskset.h"
#include "TString.h"
#include <iostream>

#include "TStopwatch.h"

#include "TH1I.h"
#include "TH2F.h"

#include "hruntimedb.h"
#include "hrecevent.h"
#include "hreconstructor.h"
#include "hcategorymanager.h"
#include "hcategory.h"


#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hfwdetcand.h"
#include "hfwdetcandsim.h"

//--------category definitions---------
#include "hparticledef.h"
#include "hstartdef.h"
#include "hgeantdef.h"
#include "hpiontrackerdef.h"
//-------------------------------------

//-------objects-----------------------
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticleevtinfo.h"
#include "hstart2hit.h"
#include "hgeantkine.h"
#include "heventheader.h"
#include "hpiontrackertrack.h"
#include "hparticletracksorter.h"

#include "hstart2cal.h"
//-------------------------------------
#include "hphysicsconstants.h"
#include "hparticletool.h"

#include "henergylosscorrpar.h"

#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <string>


double trackDistance(HParticleCand* track1, HParticleCand*  track2);
double trackDistance(HParticleCand* track1, HFwDetCand*  track2);
HGeomVector trackVertex(HParticleCand* track1, HParticleCand*  track2);
HGeomVector trackVertex(HParticleCand* track1, HFwDetCand*  track2);
Int_t getMotherIndex(HGeantKine* particle);
bool isLepton(HParticleCand* particle);

#endif
