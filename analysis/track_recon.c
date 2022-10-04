#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"

void track_recon(int event_num){
	TCanvas *c1 = new TCanvas();
	TCanvas *c2 = new TCanvas();	
	TH2D* dHist = new TH2D("dhist", "dhist", 500, 0, 500,500, 0, 500);
	TH2D* tHist = (TH2D*)dHist->Clone(); 
	TFile* f_in = new TFile("../test.root");
	TTree *t = (TTree*) f_in->Get("diffusion output");
	int max_dep=5000;
	int event_id;
	int nDep_pix;
	int x_pix[nDep_pix], y_pix[nDep_pix], z_pos[nDep_pix], z_toa[nDep_pix], depC_pix[nDep_pix]; 
	float mean_z;
	
	cout<<t->GetEntries()<<" aasdf"<<endl;
	t->SetBranchAddress("event_id", &event_id); 
	t->SetBranchAddress("nDep_pix", &nDep_pix); 
	t->SetBranchAddress("x_pix", &x_pix);
	t->SetBranchAddress("y_pix", &y_pix);
	t->SetBranchAddress("z_toa", &z_toa);
	t->SetBranchAddress("depC_pix", &depC_pix);
	t->SetBranchAddress("mean_z", &mean_z);

	dHist->GetXaxis()->SetTitle("X Pixel");
	dHist->GetYaxis()->SetTitle("Y Pixel");	
	tHist->SetName("thist");
	tHist->SetTitle("thist");
	t->GetEntry(event_num);
	int xmin = TMath::MinElement(nDep_pix, x_pix);
	int ymin = TMath::MinElement(nDep_pix, y_pix);
	int xmax = TMath::MaxElement(nDep_pix, x_pix);
	int ymax = TMath::MaxElement(nDep_pix, y_pix);
	
	int cdep = TMath::MaxElement(nDep_pix, depC_pix);
	//find idxs of max depositions and z-positions
	
	for(int i=0; i<nDep_pix; i++){
		if (depC_pix[i]== cdep){
			cout<<x_pix[i] << " "<< y_pix[i]<< endl;
		}	
		tHist->Fill(x_pix[i], y_pix[i], z_toa[i]);
		dHist->Fill(x_pix[i], y_pix[i], depC_pix[i]);
	}
	
	cout<<"Mean Z (um): " << mean_z/1000 <<endl;
	c1->cd();
	tHist->Draw("COLZ");
	c2->cd();
	dHist->Draw("COLZ");
}	
