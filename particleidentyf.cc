#include "particleidentyf.h"

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


bool isLepton(HParticleCand* particle)
{
  double delta=0.05;
  double mquality=particle->getRichMatchingQuality();

  double dphi=particle->getDeltaPhi();
  double dtheta=particle->getDeltaTheta();
      
          
  if(particle->isFlagBit(kIsUsed))
    {
      if(mquality==-1 || mquality>5)
	return false;
      if(particle->getBeta()<(1-delta) || particle->getBeta()>(1+delta))
	return false;
      // if(dtheta<-0.4 || dtheta>0.4)
      // return false;
      // if(dphi*TMath::Sin(particle->getTheta())>0.4 || dphi*TMath::Sin(particle->getTheta())<-0.4)
      // return false;
    }
  else
    return false;
     
  return true;
}



