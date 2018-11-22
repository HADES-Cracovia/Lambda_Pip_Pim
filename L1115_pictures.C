#include <sstream>
#include <string>
#include <fstream>

void sethist(TH1F* hist, int color=1, int rebin=1, int style=1, double scale=1)
{
  hist->SetLineColor(color);
  hist->Rebin(rebin);
  hist->SetLineStyle(style);
  hist->Scale(scale);
  hist->SetLineWidth(3);
}

int calcBg(TH1F* &histBG, TH1F** hists, int channels, int signal)
{
  if(signal!=1)
    histBG=(TH1F*)hists[0]->Clone("background");
  else
    histBG=(TH1F*)hists[1]->Clone("background");

  histBG->Reset();
  
  for(int i=0;i<channels;i++)
    {
      if(i==signal-1)
	continue;
      else
	histBG->Add(hists[i]);
    }
  return 1;
}

int L1115_pictures()
{
  const int fn=8;//number of files in "infile"
  std::ifstream infile("files_list_k.dat");//list of histograms
  std::string file;
  std::string directory="/lustre/nyx/hades/user/iciepal/Lambda1520_ic/withSec4/";//directory comon for all files
  int n=0;

  
  float scale[]={
    1840. ,
    130.*7.8e-5,
    130./5.34,
    300,
    300,
    43,
    10,
    7
  };//ub

  TH1F *hinvM_pmHpHDist[fn];
  TH1F *hinvM_pmHpFTDist[fn];
  
  //TH1F *hDLmassFinalRL_L[fn];
  // TH1F *hDLmassDistZLRL_L[fn];
  //TH1F *hDLmassDistZL[fn];
  //TH1F *hDLmassFTDistZL[fn];
  //TH1F *hinvM_pmHpHDist_emem[fn];
  //TH1F *hinvM_pmHpHDist_epep[fn];
 


  
  TH1F *hL1520mass_background;
  TH1F *hL1520mass_background_L;
  //TH1F *hDLmass_background;
  //TH1F *hDLmass_background_L;
  
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
      
      hinvM_pmHpHDist[n]= (TH1F*)hist_file->Get("hinvM_pmHpHDist")->Clone();
      hinvM_pmHpFTDist[n]= (TH1F*)hist_file->Get("hinvM_pmHpFTDist")->Clone();
       
      n++;
    }
  //end of reading histograms*************************************************************************

  //Print resoults************************************************************************************
  for(int k=0; k<n;k++)
    {
      int bins=1;
      //set colors
          
      sethist(hinvM_pmHpHDist[k],k+1,bins,1,scale[k]);
      sethist(hinvM_pmHpFTDist[k],k+1,bins,1,scale[k]);
      
      
      //sum FW with HADES
      
      hinvM_pmHpHDist[k]->Add(hinvM_pmHpFTDist[k]);
      
      
    }
  
  //calculate sum of all background channels
  calcBg(hL1520mass_background, hinvM_pmHpHDist, 8, 2);
  //calcBg(hDLmass_background, hDLmassDistZL, 8, 2);
  
  cout<<"all histograms set"<<endl;
  
  TLegend *legend = new TLegend(0.1,0.2,0.99,0.9);
  legend->AddEntry(hinvM_pmHpHDist[1],"p K+ #Lambda(1520)[#Lambda(1115) e+ e-] 130#mub","l");
  //legend->AddEntry(hL1520massFinalRLpi0_L[1],"p K+ #Lambda(1520)[#Lambda(1115) e+ e-] 130#mub - true signal","l");
  legend->AddEntry(hinvM_pmHpHDist[5],"p K+ #Lambda(1115) #pi^{0}  100#mub","l");
  legend->AddEntry(hinvM_pmHpHDist[6],"p K+ #Lambda(1115) 2#pi^{0} 20#mub","l");
  legend->AddEntry(hinvM_pmHpHDist[7],"p K+ #Lambda(1115) 3#pi^{0} 7#mub","l");
  legend->AddEntry(hinvM_pmHpHDist[0],"p p #pi^{+} #pi^{-} #pi^{0} 1840#mub","l");
  legend->AddEntry(hinvM_pmHpHDist[4],"p p #pi^{+} #pi^{-} 2#pi^{0} 300#mub","l");
  legend->AddEntry(hinvM_pmHpHDist[3],"p n 2#pi^{+} #pi^{-} #pi^{0} 200#mub","l");
  legend->AddEntry(hinvM_pmHpHDist[2],"L1520 decays 130#mub","l");


  //Draw all things
  TCanvas *cPictures = new TCanvas("cPictures","cPictures");

  cPictures->Divide(2);
  double ymin=1e-0; //min value for y axis  
  
  
  cPictures->cd(1);
  gPad->SetLogy();
  hinvM_pmHpHDist[0]->GetYaxis()->SetRangeUser(ymin,10e7);
  hinvM_pmHpHDist[0]->GetYaxis()->SetTitle("#sigma *  counts");
  hinvM_pmHpHDist[0]->GetXaxis()->SetTitle("Inv_M [MeV]");
  for(int x=0;x<n;x++)
    {
      hinvM_pmHpHDist[x]->Draw("same");
    }
  //hL1520mass_background->Draw("same");

  
  cPictures->cd(2);
  legend->Draw();
  
  double x_min=1100;
  double x_max=1130;

  double signal=hinvM_pmHpHDist[1]->Integral(hinvM_pmHpHDist[1]->FindBin(x_min),hinvM_pmHpHDist[1]->FindBin(x_max));
  double background=hL1520mass_background->Integral(hL1520mass_background->FindBin(x_min),hL1520mass_background->FindBin(x_max));
  
  cout<<"*********** Final raport ***********"<<endl;
  cout<<"Signal window: "<<x_min<<" - "<<x_max<<"MeV"<<endl;
  cout<<"signal: "<<signal<<" background: "<<background<<endl;
  cout<<"S/B ratio :"<< signal/background << endl;
  cout<<"significance (S/Sqrt(S+B))"<<signal/TMath::Sqrt(signal+background)<<endl;
    
  //end of printing results***************************************************************************
  return 1;
}
