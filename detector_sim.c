#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h" 
#include <algorithm>
int extracthist(TH2D* pix_hist, vector<int> bins, vector<int> t_vec, int x_pix[], int y_pix[], int depC[], int z_toa[]){
	int nbinsx = pix_hist->GetNbinsX();
	int nbinsy = pix_hist->GetNbinsY();
	//cout<<nbinsz<<endl;
	//cout<<pix_hist->GetZaxis()->GetXmin()<<endl;	
	int nDep=0;
	unordered_set<int> used_bins;
	for(int i=0; i<bins.size(); i++){
		if(used_bins.find(bins[i]) == used_bins.end()){
			used_bins.insert(bins[i]);

			int bin_val = pix_hist->GetBinContent(bins[i]);
			int x,y,z;
			pix_hist->GetBinXYZ(bins[i],x,y,z);
			z_toa[nDep]=t_vec[i];
			x_pix[nDep]=x;
			y_pix[nDep]=y;
			depC[nDep]=bin_val;
			nDep++;
			//cout<<x<<" "<<y<<" "<<bin_val<<" "<<z<<endl;
		}		
	}	
	cout<<nDep<<endl;
	return nDep;
}
double calc_toa(double z_pos){
	//drift time t= z/(mu_h * E)
	//field is 50v/um
	//mobility is 0.144 cm^2/(volts * seconds)
	//z position in um
	//1.389 nanoseconds per micron		
	return 1.389 * z_pos;
}


void detector_sim(int start_num, int end_num){
	//Apply diffusion script to all events in file from idx start_num to end_num
	//TCanvas *c1 = new TCanvas();
	
	int dWidth = 300;
	TRandom3 *diffusion = new TRandom3(801);
	TRandom3 *z_random = new TRandom3(999);
	TH2D* dHist = new TH2D("dhist", "dhist", dWidth*2, -dWidth, dWidth, dWidth*2, -dWidth, dWidth);
	TH2D* tHist = new TH2D("thist", "thist", dWidth*2, -dWidth, dWidth, dWidth*2, -dWidth, dWidth);
	
	TFile* f = new TFile("/data/damic/selena/nix-sim-results/3MeV_0dl.root");

	TTree* t= (TTree*) f->Get("out");
	:ouble diff_factor = .00323;
	int nEvents = t->GetEntries();
	int n=23000;
	Int_t event_id;
	Double_t depX[n], depZ[n], depY[n], depT[n];
	Int_t nDep;
	Int_t depC[8*n];
	int totC[8];
	//cout<<nEvents<<endl;
	t->SetBranchAddress("eventID", &event_id);
	t->SetBranchAddress("depX_true", &depX);
	t->SetBranchAddress("depY_true", &depY);
	t->SetBranchAddress("depZ_true", &depZ);
	t->SetBranchAddress("depT_true", &depT);
	t->SetBranchAddress("nDep", &nDep);
	t->SetBranchAddress("depC", &depC);
	t->SetBranchAddress("totC", &totC);	
	
	//Setup Output
	TFile f_out("test.root", "recreate");
	TTree *t_out = new TTree("diffusion output", "Diffusion test");
	int max_charge=26000;

	int nDep_pix;
	int x_pix[max_charge], y_pix[max_charge], z_toa[max_charge], depC_pix[max_charge]; 
	float mean_z;

	t_out->Branch("event_id", &event_id, "event_id/I");
	t_out->Branch("nDep_pix", &nDep_pix, "nDep_pix/I");
	t_out->Branch("x_pix", &x_pix, "x_pix[nDep_pix]/I");
	t_out->Branch("y_pix", &y_pix, "y_pix[nDep_pix]/I");
	t_out->Branch("z_toa", &z_toa, "z_toa[nDep_pix]/I");
	t_out->Branch("depC_pix", &depC_pix, "depC_pix[nDep_pix]/I");
	t_out->Branch("mean_z", &mean_z, "mean_z/F");

	
	for (int i=start_num; i<end_num; i++){
		t->GetEntry(i);
		double sum_x = 0;
		double sum_y = 0;
		double sum_z_0 = 0;
		
		double max_x=TMath::MaxElement(nDep,depX);
		double max_y=TMath::MaxElement(nDep,depY);


		double min_x=TMath::MinElement(nDep,depX);
		double min_y=TMath::MinElement(nDep,depY);
		double min_z = 0;	
		//set bounds for hist to save time in extract depositions
		
		for(int i =0; i<nDep; i++){
			sum_x += depX[i];
			sum_y += depY[i];
			sum_z_0 += depZ[i];
			if (min_z>depZ[i]){
				min_z=depZ[i];
			}

		}
		
		double mean_z_0 = sum_z_0/nDep;
		double sum_z=0;
		double mean_x = sum_x/nDep;
		double mean_y = sum_y/nDep;
		int n_charge=0;
		//Fill Histograms with diffused Track
	
		//calc hist bounds	
		int xmin_pix = TMath::Nint((min_x-mean_x)/10)-1;
		int xmax_pix = TMath::Nint((max_x-mean_x)/10)+1;
		int ymin_pix = TMath::Nint((min_y-mean_y)/10)-1;
		int ymax_pix = TMath::Nint((max_y-mean_y)/10)+1;

		dHist->SetBins(xmax_pix-xmin_pix, xmin_pix, xmax_pix, ymax_pix-ymin_pix, ymin_pix, ymax_pix);
		double zOffset = z_random->Rndm()*5000;

		bool skip_track=false;
		vector<int> pix_bins;
		vector<int> toa;	
		for(int j=0; j<nDep; j++){
			//cout<<depZ[j]<<" ";	
			int highFieldIdx = 8*j + 7;
			double xPos, yPos;
			xPos = (depX[j]-mean_x)/10.;
			yPos = (depY[j]-mean_y)/10.;
				
			//Need to consider z position relative to bottom of detector, so add positive values for tracks in negative z-plane
			if (min_z<0){
				depZ[j]= depZ[j] + (-1*min_z) + zOffset;
			}else{
				depZ[j]=depZ[j] - mean_z_0 + zOffset;// + (z_random->Rndm())*5000;
			}
			if (depZ[j]<0 || depZ[j]>5000){
				skip_track=true;
				cout<<"here"<<endl;
				break;
			}


			sum_z += depZ[j];
				
			//calculate sigma for different depths
			double sigma = diff_factor * sqrt(depZ[j]);
		 	int numCharge = depC[highFieldIdx];
			
			for (int k =0; k<numCharge; k++){

				Double_t x_1 = diffusion->Gaus(xPos, sigma);
				Double_t y_1 = diffusion->Gaus(yPos, sigma);
				toa.push_back(TMath::Floor(calc_toa(depZ[j])/5));
				pix_bins.push_back(dHist->Fill(x_1, y_1));
				
			}
		}
		mean_z = sum_z/nDep;
		if (!skip_track){
			//dHist->Draw();
			nDep_pix=extracthist(dHist, pix_bins,toa,x_pix, y_pix, depC_pix, z_toa); 		
			t_out->Fill();
		}
		if (mean_z>5000 || skip_track){
			continue;
		}
	
	}
	f_out.cd();
	t_out->Write();

}
	
