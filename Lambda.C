#include <sstream>
#include <string>
#include <fstream>


int Lambda()
{
  const int fn=8;//number of files in "infile"
  std::ifstream infile("files_list.dat");//list of histograms
  std::string file;
  std::string directory="/lustre/nyx/hades/user/iciepal/Lambda1520_ic/";//directory comon for all files
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

  TH1F *hL1520massDist[fn];
  TH1F *hL1520massFTDist[fn];
  TH1F *hL1520massDistL[fn];;
  TH1F *hL1520massFTDistL[fn];
  TH1F *hL1520massDistZL[fn];
  TH1F *hL1520massFTDistZL[fn];
  TH1F *hL1520massFinal[fn];
  TH1F *hL1520massFTFinal[fn];
  
  TH1F *hL1520massDistZLRL[fn]; 
  TH1F *hL1520massFTDistZLRL[fn];
  TH1F *hL1520massFinalRL[fn];
  TH1F *hL1520massFTFinalRL[fn];
  
  TH1F *hL1520massDistZLpi0[fn];
  TH1F *hL1520massFTDistZLpi0[fn];
  TH1F *hL1520massFinalRLpi0[fn];
  TH1F *hL1520massFTFinalRLpi0[fn];

  TH1F *hL1520mass_HHemem[fn];
  TH1F *hL1520mass_HHepep[fn];
  TH1F *hL1520mass_HFTemem[fn];
  TH1F *hL1520mass_HFTepep[fn];

  TH1F *hL1520massFinalRLpi0_L[fn];
  TH1F *hL1520massDistZLRLpi0_L[fn];
  
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
      hL1520massDist[n]= (TH1F*)hist_file->Get("hL1520massDist")->Clone();
      hL1520massFTDist[n]= (TH1F*)hist_file->Get("hL1520massFTDist")->Clone();
      hL1520massDistL[n]= (TH1F*)hist_file->Get("hL1520massDistL")->Clone();
      hL1520massFTDistL[n]= (TH1F*)hist_file->Get("hL1520massFTDistL")->Clone();
      hL1520massDistZL[n]= (TH1F*)hist_file->Get("hL1520massDistZL")->Clone();
      hL1520massFTDistZL[n]= (TH1F*)hist_file->Get("hL1520massFTDistZL")->Clone();

      hL1520massDistZLRL[n]= (TH1F*)hist_file->Get("hL1520massDistZLRL")->Clone();
      hL1520massFTDistZLRL[n]= (TH1F*)hist_file->Get("hL1520massFTDistZLRL")->Clone();
      hL1520massFinalRL[n]= (TH1F*)hist_file->Get("hL1520massFinalRL")->Clone();
      hL1520massFTFinalRL[n]= (TH1F*)hist_file->Get("hL1520massFTFinalRL")->Clone();
      hL1520massFinal[n]= (TH1F*)hist_file->Get("hL1520massFinal")->Clone();
      hL1520massFTFinal[n]= (TH1F*)hist_file->Get("hL1520massFTFinal")->Clone();

      hL1520massDistZLpi0[n]= (TH1F*)hist_file->Get("hL1520massDistZLpi0")->Clone();
      hL1520massFTDistZLpi0[n]= (TH1F*)hist_file->Get("hL1520massFTDistZLpi0")->Clone();
      hL1520massFinalRLpi0[n]= (TH1F*)hist_file->Get("hL1520massFinalRLpi0")->Clone();
      hL1520massFTFinalRLpi0[n]= (TH1F*)hist_file->Get("hL1520massFTFinalRLpi0")->Clone();

      hL1520mass_HHemem[n]= (TH1F*)hist_file->Get("hL1520mass_HHemem")->Clone();
      hL1520mass_HHepep[n]= (TH1F*)hist_file->Get("hL1520mass_HHepep")->Clone();
      hL1520mass_HFTemem[n]= (TH1F*)hist_file->Get("hL1520mass_HFTemem")->Clone();
      hL1520mass_HFTepep[n]= (TH1F*)hist_file->Get("hL1520mass_HFTepep")->Clone();


      hL1520massFinalRLpi0_L[n]=(TH1F*)hist_file->Get("hL1520massFinalRLpi0_L")->Clone();
      hL1520massDistZLRLpi0_L[n]=(TH1F*)hist_file->Get("hL1520massDistZLRLpi0_L")->Clone();
	
     n++;
    }
  //end of reading histograms*************************************************************************

  //Print resoults************************************************************************************
  for(int k=0; k<n;k++)
    {
      //set colors
      hL1520massDist[k]->SetLineColor(k+1);
      hL1520massFTDist[k]->SetLineColor(k+1);
      hL1520massDistL[k]->SetLineColor(k+1);
      hL1520massFTDistL[k]->SetLineColor(k+1);
      hL1520massDistZL[k]->SetLineColor(k+1);
      hL1520massFTDistZL[k]->SetLineColor(k+1);

      hL1520massDistZLRL[k]->SetLineColor(k+1);
      hL1520massFTDistZLRL[k]->SetLineColor(k+1);
      hL1520massFinalRL[k]->SetLineColor(k+1);
      hL1520massFTFinalRL[k]->SetLineColor(k+1);
      hL1520massFinal[k]->SetLineColor(k+1);
      hL1520massFTFinal[k]->SetLineColor(k+1);

      hL1520massDistZLpi0[k]->SetLineColor(k+1);
      hL1520massFTDistZLpi0[k]->SetLineColor(k+1);
      hL1520massFinalRLpi0[k]->SetLineColor(k+1);
      hL1520massFTFinalRLpi0[k]->SetLineColor(k+1);

      hL1520mass_HHemem[k]->SetLineColor(k+1);
      hL1520mass_HHepep[k]->SetLineColor(k+1);
      hL1520mass_HFTemem[k]->SetLineColor(k+1);
      hL1520mass_HFTepep[k]->SetLineColor(k+1);
     
      //sum FW with HADES
      hL1520massDist[k]->Add(hL1520massFTDist[k]);
      hL1520massDistL[k]->Add(hL1520massFTDistL[k]);
      hL1520massDistZL[k]->Add(hL1520massFTDistZL[k]);
      hL1520massFinal[k]->Add(hL1520massFTFinal[k]);

      hL1520massDistZLRL[k]->Add(hL1520massFTDistZLRL[k]);
      hL1520massFinalRL[k]->Add(hL1520massFTFinalRL[k]);
      hL1520massDistZLpi0[k]->Add(hL1520massFTDistZLpi0[k]);
      hL1520massFinalRLpi0[k]->Add(hL1520massFTFinalRLpi0[k]);
      
      hL1520mass_HHemem[k]->Add(hL1520mass_HFTemem[k]);
      hL1520mass_HHepep[k]->Add(hL1520mass_HFTepep[k]);
      
      //rebin histograms
      int bins=25;
      
      hL1520massDist[k]->Rebin(bins);
      hL1520massDistL[k]->Rebin(bins);
      hL1520massDistZL[k]->Rebin(bins);
      hL1520massFinal[k]->Rebin(bins);
      
      hL1520massDistZLRL[k]->Rebin(bins);
      hL1520massFinalRL[k]->Rebin(bins);
      hL1520massDistZLpi0[k]->Rebin(bins);
      hL1520massFinalRLpi0[k]->Rebin(bins);

      hL1520mass_HHemem[k]->Rebin(bins);
      hL1520mass_HHepep[k]->Rebin(bins);
      
      hL1520massDist[k]->Scale(scale[k]);
      hL1520massDistL[k]->Scale(scale[k]);
      hL1520massDistZL[k]->Scale(scale[k]);
      hL1520massFinal[k]->Scale(scale[k]);
      
      hL1520massDistZLRL[k]->Scale(scale[k]);
      hL1520massFinalRL[k]->Scale(scale[k]);
      hL1520massDistZLpi0[k]->Scale(scale[k]);
      hL1520massFinalRLpi0[k]->Scale(scale[k]);

      hL1520mass_HHemem[k]->Scale(scale[k]);
      hL1520mass_HHepep[k]->Scale(scale[k]);
      
      
      /*hL1520massFinalRLpi0_L[k]->SetLineStyle(2);
      hL1520massDistZLRLpi0_L[k]->SetLineStyle(2);
      hL1520massFinalRLpi0_L[k]->Rebin(bins);
      hL1520massDistZLRLpi0_L[k]->Rebin(bins);
      hL1520massFinalRLpi0_L[k]->Scale(scale[k]);
      hL1520massDistZLRLpi0_L[k]->Scale(scale[k]);
      */    
    }  
  TLegend *legend = new TLegend(0.1,0.2,0.99,0.9);
  legend->AddEntry(hL1520massDist[5],"p K+ #Lambda(1520)[#Lambda(1115) e+ e-] 130#mub","l");

  legend->AddEntry(hL1520massDist[2],"p K+ #Lambda(1115) #pi^{0}  100#mub","l");
  legend->AddEntry(hL1520massDist[3],"p K+ #Lambda(1115) 2#pi^{0} 20#mub","l");
  legend->AddEntry(hL1520massDist[4],"p K+ #Lambda(1115) 3#pi^{0} 7#mub","l");

  legend->AddEntry(hL1520massDist[0],"p p #pi^{+} #pi^{-} #pi^{0} 1840#mub","l");
  legend->AddEntry(hL1520massDist[1],"p p #pi^{+} #pi^{-} 2#pi^{0} 300#mub","l");
  legend->AddEntry(hL1520massDist[7],"p n 2#pi^{+} #pi^{-} #pi^{0} 200#mub","l");
  legend->AddEntry(hL1520massDist[6],"L1520 decays 130#mub","l");

  TCanvas *cEpEm = new TCanvas("cEpEm","cEpEm");

  cEpEm->Divide(4,3);
  double ymin=1e-3; //min value for y axis  
  
  cEpEm->cd(1);
  gPad->SetLogy();
  hL1520massDist[0]->GetYaxis()->SetRangeUser(ymin,10e6);
  for(int x=0;x<n;x++)
    hL1520massDist[x]->Draw("same");
      
  cEpEm->cd(2);
  gPad->SetLogy();
  hL1520massDistL[0]->GetYaxis()->SetRangeUser(ymin,10e6);
  for(int x=0;x<n;x++)
    hL1520massDistL[x]->Draw("same");
  

  cEpEm->cd(3);
  gPad->SetLogy();
  hL1520massDistZL[0]->GetYaxis()->SetRangeUser(ymin,10e6);
  for(int x=0;x<n;x++)
    hL1520massDistZL[x]->Draw("same");
  

  cEpEm->cd(4);
  gPad->SetLogy();
  hL1520massFinal[0]->GetYaxis()->SetRangeUser(1e-6,10e5);
  for(int x=0;x<n;x++)
    hL1520massFinal[x]->Draw("same");
  
  cEpEm->cd(5);
  gPad->SetLogy();
  hL1520massDistZLpi0[0]->GetYaxis()->SetRangeUser(1e-6,10e5);
  for(int x=0;x<n;x++)
    hL1520massDistZLpi0[x]->Draw("same");
  //hL1520massDistZLpi0_L[5]->Draw("same");

  cEpEm->cd(6);
  gPad->SetLogy();
  hL1520massFinalRLpi0[0]->GetYaxis()->SetRangeUser(1e-6,10e5);
  for(int x=0;x<n;x++)
    hL1520massFinalRLpi0[x]->Draw("same");
  //hL1520massFinalRLpi0_L[5]->Draw("SAME");
	 
  cEpEm->cd(7);
  gPad->SetLogy(7);
  hL1520mass_HHemem[0]->GetYaxis()->SetRangeUser(1e-6,10e4);
  for(int x=0;x<n;x++)
    hL1520mass_HHemem[x]->Draw("same");

  cEpEm->cd(8);
  gPad->SetLogy(8);
  hL1520mass_HHepep[0]->GetYaxis()->SetRangeUser(1e-6,10e4);
  for(int x=0;x<n;x++)
    hL1520mass_HHepep[x]->Draw("same");

  
  cEpEm->cd(9);
  legend->Draw();
  
  //end of printing results***************************************************************************
  return 1;
}
