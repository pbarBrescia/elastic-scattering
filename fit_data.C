#include "get_stdout.h"
// #include <TMinuit.h>
// #include <TObjArray.h>
// #include <TObjString.h>
// #include <TStopwatch.h>
// #include <TString.h>
// #include <TCanvas.h>
// #include <TGraph.h>
// #include <TGraphErrors.h>
// #include <TROOT.h>
// #include <TStyle.h>

// #include <vector>
// #include <cstdlib>
// #include <fstream>


using namespace std;

// Int_t ndata;
vector<Double_t> x;
vector<Double_t> data;
vector<Double_t> error;
Double_t mom, at, zt, charge;
string opt;

Double_t csline(Double_t x,Double_t *par) {
	// run an external command and get its output in a TString
	// dependence on momentum (see Lee-Wong paper)
	Double_t dep = 1.0;//cosh(sqrt(14.04+mom)-sqrt(14.04))/cosh(sqrt(7.92+mom)-sqrt(7.92));
	TString s;
	// TString s = get_stdout(Form("./bin/antip_scan %f %f %f %f %f %f %f 1.25 %f %f ang %f",mom,at,zt,par[0]*dep,par[1],par[2],par[2],par[3],par[3],x)); // 4 free param
	if(opt == "ang") s = get_stdout(Form("./bin/antip_scan %f %f %f %f %f %f %f 1.25 %f %f %s %f %f",mom,at,zt,par[0]*dep,par[1],par[2],par[3],par[4],par[5],opt.c_str(),x,charge));
	else if(opt == "mom") s = get_stdout(Form("./bin/antip_scan %f %f %f %f %f %f %f 1.25 %f %f %s %f %f",x,at,zt,par[0]*dep,par[1],par[2],par[3],par[4],par[5],opt.c_str(),999.,charge));
	// split the TString using space as delimiter
	// Example taken from here: https://root-forum.cern.ch/t/split-tstring-by-delimeter-in-root-c/18228/2
	TObjArray *ts = s.Tokenize(" ");
	// get the 1st element in a Double_t and returns to the program
	Double_t ret = atof( ((TObjString *) (ts->At(1)) )->String());
	return ret;
}

Double_t csline(Double_t* x,Double_t *par) {
	// run an external command and get its output in a TString
	// dependence on momentum (see Lee-Wong paper)
	Double_t dep = 1.0;//cosh(sqrt(14.04+mom)-sqrt(14.04))/cosh(sqrt(7.92+mom)-sqrt(7.92));
	TString s;
	// TString s = get_stdout(Form("./bin/antip_scan %f %f %f %f %f %f %f 1.25 %f %f ang %f",mom,at,zt,par[0]*dep,par[1],par[2],par[2],par[3],par[3],x[0])); // 4 free param
	if(opt == "ang") s = get_stdout(Form("./bin/antip_scan %f %f %f %f %f %f %f 1.25 %f %f %s %f %f",mom,at,zt,par[0]*dep,par[1],par[2],par[3],par[4],par[5],opt.c_str(),x[0],charge));
	else if(opt == "mom") s = get_stdout(Form("./bin/antip_scan %f %f %f %f %f %f %f 1.25 %f %f %s %f %f",x[0],at,zt,par[0]*dep,par[1],par[2],par[3],par[4],par[5],opt.c_str(),999.,charge));
	// split the TString using space as delimiter
	// Example taken from here: https://root-forum.cern.ch/t/split-tstring-by-delimeter-in-root-c/18228/2
	TObjArray *ts = s.Tokenize(" ");
	// get the 1st element in a Double_t and returns to the program
	Double_t ret = atof( ((TObjString *) (ts->At(1)) )->String());
	return ret;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   Int_t i;

//calculate chisquare
   Double_t chisq = 0;
   Double_t delta,den1,den2;
   Int_t ndata = (unsigned) x.size();
   for (i=0; i<ndata; i++) {
   	 // den1 = error[i];
   	 // den2 = derivative(x[i],par)*0.2;//*(data[i]/csline(x[i],par))*0.2;
     delta  = (data[i]-csline(x[i],par))/error[i];
     chisq += delta*delta;
     // chisq += delta*delta/(den1*den1+den2*den2);
   }
   f = chisq;
}

Double_t error_prop(Int_t npar, TVirtualFitter *gMinuit, TF1 *func, Double_t *xx)
{
	Double_t errf = 0.,sqerrf = 0.;
	for(int j = 0; j < npar; j++)
	{
		Double_t covmatdiag = gMinuit->GetCovarianceMatrixElement(j,j);
		// Double_t parerror = gMinuit->GetParError(j);
		errf += pow((func->GradientPar(j,xx)),2.)*covmatdiag;
		// sqerrf += errf*errf;
	}

	return sqrt(errf);
}

int fit_data(string datafile="data/ca40_1798_err.dat",Double_t Charge = -1., Double_t Mom = 599.91, Double_t At = 40.078, Double_t Zt = 20.0, string Opt = "ang"){

	auto start = chrono::system_clock::now();
	time_t start_time = chrono::system_clock::to_time_t(start);
	cout << "Started at " << ctime(&start_time) << endl;
	// record the execution time
	TStopwatch t;

	// Graphics option
	// true = plot directly the fit results vs data
	// false = no plots (default)
	// WARNING: CHOOSING TRUE HAS GREAT COMPUTATIONAL TIME COSTS!
	bool WithGraphics = false;

	// Classic option
	// true = use "Fit" method from TGraph(Errors) (default)
	// false = minimization with TMinuit and saves
	//		   the results in a TF1.
	// N.B: false opt cannot do some things (calculation of 95% CI)
	bool ClassicFit = true;

	// open the file with exp. data
	ifstream in(datafile.c_str(), ifstream::in );

	charge = Charge; // projectile charge
	mom = Mom; // N.B: for option "mom" must be the initial one (e.g. 50.0 -> from 50.0 MeV/c to ...)
	at = At;
	zt = Zt;
	opt = Opt; // ang or mom

	string pp;

	cout << "================ I N F O ===============" << endl;
	cout << "Data from file " << datafile.c_str() << endl;
	if(opt == "mom") pp = "p_start = ";
	else if(opt == "ang") pp = "p = \t"; 
	cout << pp.c_str() << mom << " MeV/c \n" 
		 << "A_t = \t" << at << " amu \n"
		 << "Z_t = \t" << zt << "\n"
		 << "Opt = \t" << opt.c_str() << "\n"
		 << "q_p = \t" << charge << " e \n";
	cout << "========================================" << endl;

	Int_t i = 0;
	const Int_t npar = 6;

	while(!in.eof())
	{
		Double_t x_i, data_i, error_i;
		in >> x_i >> data_i >> error_i;
		x.push_back(x_i);
		data.push_back(data_i);
		error.push_back(error_i);
		i++;
	}
	// ndata=i;
	in.close();

   	Double_t arglist[10];

   	static vector<Double_t> vstart;
   	static vector<Double_t> step;
// Set starting values and step sizes for parameters
   	if(npar==6)
   	{
   	 	vstart={50.0, 50.0, 1.25, 1.25, 0.5, 0.5};
   	 	step={0.001 , 0.001 , 0.001 , 0.001, 0.001, 0.001};
   	}
   	else if(npar==4)
   	{ 	
   		vstart={50.0, 50.0, 1.25, 0.5};
   		step={0.001, 0.001, 0.001, 0.001};
   	}

   	// set x-axis max and min for TF1
   	Double_t xmax, xmin;
   	if(opt == "ang") {xmin = 1.; xmax = 180.;} // deg
   	else if (opt == "mom") {xmin = mom; xmax = 500.;} // MeV/c

   	TF1* func = new TF1("func",csline,xmin,xmax,npar);
   	func->SetNpx(500); // std are 100
   	func->SetParameters(&vstart[0]);

   	if(npar==6) func->SetParNames("U0", "W0", "R0r", "R0i", "A0r", "A0i");
   	else if(npar==4) func->SetParNames("U0", "W0", "R0", "A0");

	func->SetParLimits(0,0.0,5000.);
	func->SetParLimits(1,0.0,5000.);
	func->SetParLimits(2,0.8,2.5);
	func->SetParLimits(3,0.8,2.5);
	func->SetParLimits(4,0.2,0.8);
   	func->SetParLimits(5,0.2,0.8);

   	// set output filename
   	string outname = datafile;
	outname.erase(outname.begin(),outname.begin()+5); // delete the initial part of the path ("data/")
	outname.erase(outname.end()-4,outname.end()); // delete format (.dat or .txt)

   	// data in TGraphErrors
   	TGraphErrors * g3 = new TGraphErrors(datafile.c_str(),"%lg %lg %lg");
	// g3->SetTitle(Form("#bar{p} + <nucleus> @ %.2f MeV/c",mom));
	func->SetNDF(g3->GetN()-npar);

   	string typeFitter = "Genetic";
   	TVirtualFitter *gMinuit = TVirtualFitter::Fitter(0,npar);
   	Double_t vliminf[6] = {0,0,0,0,0,0};
	Double_t vlimsup[6] = {0,0,0,0,0,0};
   	if(typeFitter != "Minuit")
   	{
   		ROOT::Math::MinimizerOptions::SetDefaultMinimizer(typeFitter.c_str());
   		ROOT::Math::MinimizerOptions::SetDefaultTolerance(1.); // edm_max = 2*tolerance*1.e-3
   		ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1); // print only the final result + warnings
   	}
   	else if(typeFitter == "Minuit")
   	{
	   	TVirtualFitter::SetDefaultFitter(typeFitter.c_str());
	   	TVirtualFitter* actMinuit = (TVirtualFitter*) (TVirtualFitter::GetFitter());
   		for(int j = 0; j < npar; j++) actMinuit->SetParameter(j, func->GetParName(j), func->GetParameter(j), (func->GetParameter(j))*3./10., vliminf[j], vlimsup[j]);
   	}
   	gMinuit->SetFCN(fcn);
   	
   	if(ClassicFit) // fit with Fit method of TGraph(Errors)
	{
		TFitResultPtr fitres = g3->Fit(func,"U B S"); // user function and save fit result in the pointer
		TF1* fitfct = (TF1*) func->Clone("fitfct");
		fitfct->SetParameters(fitres->GetParams());
		fitfct->SetParErrors(fitres->GetErrors());
		fitfct->SetChisquare(fitres->Chi2());
		fitfct->SetNDF((g3->GetN())-(func->GetNpar()));
   		TFile* fff = new TFile(Form("results/fitres_%s_000.root",outname.c_str()),"recreate","fit results a TFitResultPtr");
   		fitres->Write("fitresptr");
   		func->Write("fitfct");
   		fff->Close();
   	}
   	else // fit with TVirtualFitter using TMinuit
   	{
	   	TFile* f = new TFile(Form("results/fit_%s_00.root",outname.c_str()),"recreate","file with fit results");
	// first two par free, others fixed
	   	gMinuit->FixParameter(2);
	   	gMinuit->FixParameter(3);
	   	gMinuit->FixParameter(4);
	   	gMinuit->FixParameter(5);
	  
	   	arglist[0] = 2; // print level 2 -> progresses
	   	gMinuit->ExecuteCommand("SET PRINT",arglist,1); 
	   	// gMinuit->ExecuteCommand("SET STR",arglist,1);
	   	string min = "MIGRAD";

	// Now ready for minimization step
	    arglist[0] = 5000; // max iterations
	    arglist[1] = 1.e-3; // Convergence at EDM < 1e-3 * arglist[1]
	    cout << "Executing MIGRAD Minimization..." << endl;
	// minimize the first two par, then release one-by-one 
	// the others and minimize again
	    gMinuit->ExecuteCommand(min.c_str(),arglist,2);
	    gMinuit->ReleaseParameter(2);
	    gMinuit->ExecuteCommand(min.c_str(),arglist,2);
	    gMinuit->ReleaseParameter(3);
	    gMinuit->ExecuteCommand(min.c_str(),arglist,2);
	    gMinuit->ReleaseParameter(4);
	    gMinuit->ExecuteCommand(min.c_str(),arglist,2);
	    gMinuit->ReleaseParameter(5);
	 	gMinuit->ExecuteCommand(min.c_str(),arglist,2);

	// Print results
	    Double_t minParams[npar], parErrors[npar];
	    gMinuit->SetErrorDef(1.); //-> 1-sigma 
	    for(int j = 0; j < npar; j++)
	    {
	    	minParams[j] = gMinuit->GetParameter(j);
	    	parErrors[j] = gMinuit->GetParError(j);
	    }

	    double chi2, edm, errdef;
	  	int nvpar, nparx;
	  	gMinuit->GetStats(chi2,edm,errdef,nvpar,nparx);

	// save infos of fit in TF1 and TMatrixDSym
		func->SetParameters(minParams);
		func->SetParErrors(parErrors);
		func->SetChisquare(chi2);
		int ndf = x.size()-nvpar;
		func->SetNDF(ndf);
		TMatrixDSym* covmat = new TMatrixDSym(npar,gMinuit->GetCovarianceMatrix()); 

		TGraphErrors* gerr = new TGraphErrors(500);

		if(WithGraphics) // plots (if activated)
		{
			gerr->SetTitle(Form("#bar{p} + <nucleus> @ %.2f MeV/c",mom));

			Double_t sig=1.-0.025;
			Double_t t95 = ROOT::Math::tdistribution_quantile(sig,ndf-1);
			for(int j = 0; j < 500; j++)
			{
				Double_t xx[1] = {(j+1)*(1./500.)*60.};
				gerr->SetPoint(j,xx[0],func->Eval(xx[0]));
				Double_t errf = error_prop(npar,gMinuit,func,xx);
				gerr->SetPointError(j,0,t95*errf);
			}

			gStyle->SetOptFit(1111);
			TCanvas* c1 = new TCanvas("c1","c1",800,800);
			gPad->SetLogy();
			c1->cd();
			gerr->GetYaxis()->SetRangeUser(1e-3,1e6);
			gerr->SetFillStyle(1001);
			gerr->SetFillColorAlpha(kRed,0.5);
			gerr->Draw("AE3");
			g3->SetMarkerStyle(20);
			g3->SetMarkerSize(0.5);
			g3->Draw("PSAME");

			c1->SaveAs(Form("results/%s_%.3f_plots.root",outname.c_str(),minParams[0]));

			g3->Write("gdata");
			gerr->Write("gerr");
		}  

		cout << '\a'; // sound at the end

		func->Write("fitfct");
		covmat->Write("covmat");
		f->Close(); 
	}

	t.Stop();
	t.Print("u");

	return 0;
}