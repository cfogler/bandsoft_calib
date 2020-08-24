#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"

#include "constants.h"
#include "calib_helper.h"
#include "bandhit.h"

// For processing data

using namespace std;
shiftsReader shifts;

// Neutron range of interest
const double SIGNAL_MIN = 5.8;
const double SIGNAL_MAX = 16;
const double BACKGROUND_MIN = -5;
const double BACKGROUND_MAX = 2;

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

int main(int argc, char ** argv){
	// Set style
	gStyle->SetOptFit(1);

	if (argc < 3){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./code [outputTxtfile] [outputPDFfile] [inputDatafiles]\n";
		return -1;
	}

	// Histograms for every run that sum over all bars
	const int nRuns = argc-3;
	TH1D ** ToF_spec 	= new TH1D*[nRuns];
	TH1D ** ToF_bkg 	= new TH1D*[nRuns];
	TF1 ** ToF_fits		= new TF1*[nRuns];
	TF1 ** ToF_fits_it	= new TF1*[nRuns];
	TCanvas ** cRun		= new TCanvas*[nRuns];
	for(int run = 0; run < nRuns; run++){
		ToF_spec[run] = new TH1D( Form("ToF_spec_%i",run), Form("ToF_spec_%i",run), 500, -5, 45);
		ToF_bkg[run]  = new TH1D( Form("ToF_bkg_%i",run),  Form("ToF_bkg_%i",run),  500, -5, 45);
	}

	// Load the bar shifts that were calculated at the initial bar code
	shifts.LoadInitBarFadc("../include/FADC_pass1v0_initbar.txt");
	shifts.LoadInitRunFadc("../include/FADC_pass1v0_initrun.txt");
	double * FADC_INITBAR = (double*) shifts.getInitBarFadc();
	double * FADC_INITRUN = (double*) shifts.getInitRunFadc();

	double CountsPerCharge[100000] = {0.};
	ofstream out_file;
	out_file.open(argv[1]);
	TCanvas * c0 = new TCanvas("c0","c0",900,900);
	TString openname = string(argv[2]) + "(";
	c0 -> Print(openname);

	// Loop over all the files that are given to me to get the best statistics per bar
	for( int i = 3 ; i < argc ; i++ ){
		TFile * inFile = new TFile(argv[i]);
		if (!(inFile->GetListOfKeys()->Contains("neutrons"))){
			cerr << "File has no entries\n";
			return -2;
		}
		TTree * inTree = (TTree*)inFile->Get("neutrons");

		// Initialize the input branches
		int Runno		= 0;
		double Ebeam		= 0;
		double gated_charge	= 0;
		double livetime		= 0;
		double starttime	= 0;
		int nMult		= 0;
		bandhit * nHit = new bandhit[maxNeutrons];
		inTree->SetBranchAddress("Runno"		,&Runno			);
		inTree->SetBranchAddress("Ebeam"		,&Ebeam			);
		inTree->SetBranchAddress("gated_charge"		,&gated_charge		);
		inTree->SetBranchAddress("livetime"		,&livetime		);
		inTree->SetBranchAddress("starttime"		,&starttime		);
		//	Neutron branches:
		inTree->SetBranchAddress("nMult"		,&nMult			);
		inTree->SetBranchAddress("nHit"			,&nHit			);

		// Start working on one of the files, looping over all of the events
		cout << "Working on file: " << argv[i] << "\n";
		for( int ev = 0 ; ev < inTree->GetEntries() ; ev++ ){
			if( ev % 1000000 == 0 ) cout << "\ton event " << ev << "\n";

			// Clear all branches before getting the entry from tree
			gated_charge	= 0;
			livetime	= 0;
			starttime 	= 0;
			nMult		= 0;
			for( int thishit = 0; thishit < maxNeutrons ; thishit++)
				nHit[thishit].Clear();
			
			inTree->GetEntry(ev);

			if( nMult != 1 ) continue;
			if( nHit[0].getStatus() != 0 ) continue;
			if( nHit[0].getTofFadc() == 0 ) continue;
			if( nHit[0].getEdep() < 2200*5 ) continue;


			double time = nHit[0].getTofFadc();
			double dL_n = sqrt( pow(nHit[0].getX(),2) + pow(nHit[0].getY(),2) + pow(nHit[0].getZ(),2) )/100.;


			double tof = time;// - FADC_INITBAR[(int)nHit[0].getBarID()] - FADC_INITRUN[Runno];
			tof = tof/dL_n;

			if( tof >= SIGNAL_MIN && tof <= SIGNAL_MAX ){
				ToF_spec[i-3]->Fill(tof);
				ToF_spec[i-3]->SetTitle(Form("ToF_spec_%i",Runno));
			}
			else if( tof >= BACKGROUND_MIN && tof <= BACKGROUND_MAX ){
				ToF_bkg[i-3]->Fill(tof);
				ToF_bkg[i-3]->SetTitle(Form("ToF_bkg_%i",Runno));

			}



		} // end loop over events

		double integral = ToF_spec[i-3]->Integral();
		double integral_bkg = ToF_bkg[i-3]->Integral();
		cRun[i-3] = new TCanvas(Form("Run %i",(i-3)),Form("Run %i",(i-3)),900,900);
		gated_charge /= 1000; // mC
		out_file << (Runno) << " "  
			<< integral << " " << gated_charge << " " << integral/gated_charge << " " << integral_bkg/gated_charge << "\n";

		cRun[i-3]->cd(1);
		ToF_spec[i-3]->Draw();
		ToF_bkg[i-3]->Draw("same");
		cRun[i-3]->Update();
		cRun[i-3]->Modified();
		cRun[i-3] -> Print(argv[2]);


		inFile->Close();
	}// end loop over files

	TString endname = string(argv[2]) + ")";
	c0 -> Print(endname);
	out_file.close();



	return 0;
}
