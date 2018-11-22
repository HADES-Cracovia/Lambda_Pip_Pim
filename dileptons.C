#include <sstream>
#include <string>
#include <fstream>


int dileptons()
{
  const int fn=8;//number of files in "infile"
  std::ifstream infile("files_list.dat");//list of histograms
  std::string file;
  std::string directory="/lustre/nyx/hades/user/iciepal/Lambda1520_ic/withSec2/";//directory comon for all files
  int n=0;
    
  float scale[]={
    1840.,
    300. ,
    100.,
    20.,
    7.,
    130.*7.8e-5,
    130./5.34,
    200,
  };//ub
  
  TH1F *hDLmassDist[fn];
  TH1F *hDLmassFTDist[fn];
  TH1F *hDLmassDistL[fn];;
  TH1F *hDLmassFTDistL[fn];
  TH1F *hDLmassDistZL[fn];
  TH1F *hDLmassFTDistZL[fn];
  TH1F *hDLmassFinal[fn];
  TH1F *hDLmassFTFinal[fn];
  TH1F *hinvMass_HHepep[fn];
  TH1F *hinvMass_HFTepep[fn];
  TH1F *hL1520mass_HHepep[fn];
  TH1F *hL1520mass_HFTepep[fn];
  TH1F *hinvMass_HHemem[fn];
  TH1F *hinvMass_HFTemem[fn];
  TH1F *hL1520mass_HHemem[fn];
  TH1F *hL1520mass_HFTemem[fn];
  
  //read file wth list of histograms
  if(!infile)
    {
      cout<<"Can't open the file"<<endl;
      return 1;
    }

  //loop over files names and clone histograms*********************************************************
  while (std::getline(infile,file))
    {
      cout<<"file number: "<<n<<" ";
      cout<<file<<endl;
      
      const char* file_name=(directory+file).c_str();
      TFile* hist_file=new TFile(file_name);
      if(hist_file->IsZombie())
	{
	  cout<<"can't open a file!"<<endl<<file_name<<endl;
	  n++;
	  continue;
	}
      hist_file->cd();
      hDLmassDist[n]= (TH1F*)hist_file->Get("hDLmassDist")->Clone();
      hDLmassFTDist[n]= (TH1F*)hist_file->Get("hDLmassFTDist")->Clone();
      hDLmassDistL[n]= (TH1F*)hist_file->Get("hDLmassDistL")->Clone();
      hDLmassFTDistL[n]= (TH1F*)hist_file->Get("hDLmassFTDistL")->Clone();
      hDLmassDistZL[n]= (TH1F*)hist_file->Get("hDLmassDistZL")->Clone();
      hDLmassFTDistZL[n]= (TH1F*)hist_file->Get("hDLmassFTDistZL")->Clone();

      hDLmassFinal[n]= (TH1F*)hist_file->Get("hDLmassFinal")->Clone();
      hDLmassFTFinal[n]= (TH1F*)hist_file->Get("hDLmassFTFinal")->Clone();

      hinvMass_HHepep[n]= (TH1F*)hist_file->Get("hinvMass_HHepep")->Clone();
      hinvMass_HFTepep[n]= (TH1F*)hist_file->Get("hinvMass_HFTepep")->Clone();
      hL1520mass_HHepep[n]= (TH1F*)hist_file->Get("hL1520mass_HHepep")->Clone();
      hL1520mass_HFTepep[n]= (TH1F*)hist_file->Get("hL1520mass_HFTepep")->Clone();
      hinvMass_HHemem[n]= (TH1F*)hist_file->Get("hinvMass_HHemem")->Clone();
      hinvMass_HFTemem[n]= (TH1F*)hist_file->Get("hinvMass_HFTemem")->Clone();
      hL1520mass_HHemem[n]= (TH1F*)hist_file->Get("hL1520mass_HHemem")->Clone();
      hL1520mass_HFTemem[n]= (TH1F*)hist_file->Get("hL1520mass_HFTemem")->Clone();
      n++;
    }
  //end of reading histograms*************************************************************************

  //Print resoults************************************************************************************
  for(int k=0; k<n;k++)
    {
      //set colors
      hDLmassDist[k]->SetLineColor(k+1);
      hDLmassFTDist[k]->SetLineColor(k+1);
      hDLmassDistL[k]->SetLineColor(k+1);
      hDLmassFTDistL[k]->SetLineColor(k+1);
      hDLmassDistZL[k]->SetLineColor(k+1);
      hDLmassFTDistZL[k]->SetLineColor(k+1);

      hDLmassFinal[k]->SetLineColor(k+1);
      hDLmassFTFinal[k]->SetLineColor(k+1);

      hinvMass_HHepep[k]->SetLineColor(k+1);
      hinvMass_HFTepep[k]->SetLineColor(k+1);
      hL1520mass_HHepep[k]->SetLineColor(k+1);
      hL1520mass_HFTepep[k]->SetLineColor(k+1);
      hinvMass_HHemem[k]->SetLineColor(k+1);
      hinvMass_HFTemem[k]->SetLineColor(k+1);
      hL1520mass_HHemem[k]->SetLineColor(k+1);
      hL1520mass_HFTemem[k]->SetLineColor(k+1);

      hinvMass_HHepep[k]->SetLineStyle(2);
      hinvMass_HFTepep[k]->SetLineStyle(2);
      hL1520mass_HHepep[k]->SetLineStyle(2);
      hL1520mass_HFTepep[k]->SetLineStyle(2);
      hinvMass_HHemem[k]->SetLineStyle(5);
      hinvMass_HFTemem[k]->SetLineStyle(5);
      hL1520mass_HHemem[k]->SetLineStyle(5);
      hL1520mass_HFTemem[k]->SetLineStyle(5);

      //sum FW with HADES
      hDLmassDist[k]->Add(hDLmassFTDist[k]);
      hDLmassDistL[k]->Add(hDLmassFTDistL[k]);
      hDLmassDistZL[k]->Add(hDLmassFTDistZL[k]);
      hDLmassFinal[k]->Add(hDLmassFTFinal[k]);

      hinvMass_HHepep[k]->Add(hinvMass_HFTepep[k]);
      hL1520mass_HHepep[k]->Add(hL1520mass_HFTepep[k]);
      hinvMass_HHemem[k]->Add(hinvMass_HFTemem[k]);
      hL1520mass_HHemem[k]->Add(hL1520mass_HFTemem[k]);
      
      //rebin distance histograms
      int bins=20;
      
      hDLmassDist[k]->Rebin(bins);
      hDLmassDistL[k]->Rebin(bins);
      hDLmassDistZL[k]->Rebin(bins);
      hDLmassFinal[k]->Rebin(bins);

      hinvMass_HHepep[k]->Rebin(bins);
      hL1520mass_HHepep[k]->Rebin(bins);
      hinvMass_HHemem[k]->Rebin(bins);
      hL1520mass_HHemem[k]->Rebin(bins);

      hDLmassDist[k]->Scale(scale[k]);
      hDLmassDistL[k]->Scale(scale[k]);
      hDLmassDistZL[k]->Scale(scale[k]);
      hDLmassFinal[k]->Scale(scale[k]);

      hinvMass_HHepep[k]->Scale(scale[k]);
      hL1520mass_HHepep[k]->Scale(scale[k]);
      hinvMass_HHemem[k]->Scale(scale[k]);
      hL1520mass_HHemem[k]->Scale(scale[k]);
      
    }

  TLegend *legend = new TLegend(0.1,0.2,0.99,0.9);
  legend->AddEntry(hDLmassDist[5],"p K+ #Lambda(1520)[#Lambda(1115) e+ e-] 130#mub","l");

  legend->AddEntry(hDLmassDist[2],"p K+ #Lambda(1115) #pi^{0}  100#mub","l");
  legend->AddEntry(hDLmassDist[3],"p K+ #Lambda(1115) 2#pi^{0} 20#mub","l");
  legend->AddEntry(hDLmassDist[4],"p K+ #Lambda(1115) 3#pi^{0} 7#mub","l");

  legend->AddEntry(hDLmassDist[0],"p p #pi^{+} #pi^{-} #pi^{0} 1840#mub","l");
  legend->AddEntry(hDLmassDist[1],"p p #pi^{+} #pi^{-} 2#pi^{0} 300#mub","l");
  legend->AddEntry(hDLmassDist[7],"p n 2#pi^{+} #pi^{-} #pi^{0} 200#mub","l");
  legend->AddEntry(hDLmassDist[6],"L1520 decays 130#mub","l");

  TCanvas *cEpEm = new TCanvas("L1520dileptons","L1520dileptons",1400,500);

  cEpEm->Divide(3,2);
  double ymin=1e-3; //min value for y axis  

  TH1F *background_sum;
  int x_sig=5; //no. of signall channel

  
  cEpEm->cd(1);
  gPad->SetLogy();
  hDLmassDist[0]->GetYaxis()->SetRangeUser(ymin,10e5);
  background_sum=(TH1F*)hDLmassDist[0]->Clone();

  for(int x=0;x<n;x++)
    {
      hDLmassDist[x]->Draw("same");
      if(x!=x_sig && x>0)
	background_sum->Add(hDLmassDist[x]);
    }
  background_sum->Draw("same");
  
  cEpEm->cd(2);
  gPad->SetLogy();
  hDLmassDistL[0]->GetYaxis()->SetRangeUser(ymin,10e4);
  for(int x=0;x<n;x++)
    hDLmassDistL[x]->Draw("same");
  

  cEpEm->cd(3);
  gPad->SetLogy();
  hDLmassDistZL[0]->GetYaxis()->SetRangeUser(ymin,10e3);
  for(int x=0;x<n;x++)
    hDLmassDistZL[x]->Draw("same");
  

  cEpEm->cd(4);
  gPad->SetLogy();
  hDLmassFinal[0]->GetYaxis()->SetRangeUser(1e-6,10e3);
  for(int x=0;x<n;x++)
    hDLmassFinal[x]->Draw("same");
  
  cEpEm->cd(5);
  gPad->SetLogy();
  hL1520mass_HHepep[0]->GetYaxis()->SetRangeUser(1e-6,10e2);
  hL1520mass_HHemem[0]->GetYaxis()->SetRangeUser(1e-6,10e2);
  for(int x=0;x<n;x++)
    {
      hL1520mass_HHepep[x]->Draw("same");
      hL1520mass_HHemem[x]->Draw("same");
    }
  cEpEm->cd(6);
  legend->Draw();
  
  //end of printing results***************************************************************************
  return 1;
}



