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
vector<Double_t> x;// = new Double_t[ndata];
vector<Double_t> data;// = new Double_t[ndata]; 
vector<Double_t> error;// = new Double_t[ndata];
Double_t mom, at, zt;

Double_t csline(Double_t x,Double_t *par) {
	// run an external command and get its output in a TString
	// TString s = get_stdout(Form("./p_scan.sh -i 599.91 40.078 20.0 40.5 111.0 1.1 1.1 1.25 0.63 %f ang %f",par[0]/*, par[1], par[2], par[3], par[4], par[5]*/,x[0]));
	TString s = get_stdout(Form("./bin/antip_scan %f %f %f %f %f %f %f 1.25 %f %f ang %f",mom,at,zt,par[0],par[1],par[2],par[3],par[4],par[5],x));
	// split the TString using space as delimiter
	// Example taken from here: https://root-forum.cern.ch/t/split-tstring-by-delimeter-in-root-c/18228/2
	TObjArray *ts = s.Tokenize(" ");
	// get the 1st element in a Double_t and returns to the program
	Double_t ret = atof( ((TObjString *) (ts->At(1)) )->String());
	return ret;
}
Double_t csline(Double_t* x,Double_t *par) {
	// run an external command and get its output in a TString
	// TString s = get_stdout(Form("./p_scan.sh -i 599.91 40.078 20.0 40.5 111.0 1.1 1.1 1.25 0.63 %f ang %f",par[0]/*, par[1], par[2], par[3], par[4], par[5]*/,x[0]));
	TString s = get_stdout(Form("./bin/antip_scan %f %f %f %f %f %f %f 1.25 %f %f ang %f",mom,at,zt,par[0],par[1],par[2],par[3],par[4],par[5],x[0]));
	// split the TString using space as delimiter
	// Example taken from here: https://root-forum.cern.ch/t/split-tstring-by-delimeter-in-root-c/18228/2
	TObjArray *ts = s.Tokenize(" ");
	// get the 1st element in a Double_t and returns to the program
	Double_t ret = atof( ((TObjString *) (ts->At(1)) )->String());
	return ret;
}

// Double_t derivative(Double_t x, Double_t *par, Int_t npar)
// {
// 	Double_t h[6] = {0.};
// 	h[npar] = 1e-4;
// 	Double_t parplus2[6], parplus[6], parminus[6], parminus2[6];
// 	for(int j=0; j<6; j++)
// 	{
// 		parplus[j] = par[j]+h[j];
// 		parplus2[j] = par[j]+2*h[j];
// 		parminus[j] = par[j]-h[j];
// 		parminus2[j] = par[j]-2*h[j]; 
// 	}
// 	Double_t D = (-csline(x,parplus2)+8*csline(x,parplus)-8*csline(x,parminus)+csline(x,parminus2))/(12*h[npar]);
// 	return D;
// }

// Double_t derivative(Double_t x, Double_t *par)
// {
// 	Double_t h = 1e-4;
// 	Double_t D = (-csline(x+2*h,par)+8*csline(x+h,par)-8*csline(x-h,par)+csline(x-2*h,par))/(12*h);
// 	return D;
// }

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

Double_t error_prop(Int_t npar, TVirtualFitter *gMinuit, TF1 *func, Double_t *x, Double_t eps = 0.01)
{
	Double_t errf = 0.;
	Double_t covmat[6][6];
	for(int j = 0; j < npar; j++)
	{
		for(int k = 0; k < npar; k++)
		{
			covmat[j][k] = gMinuit->GetCovarianceMatrixElement(j,k);
			errf += (func->GradientPar(j,x,eps))*(func->GradientPar(k,x,eps))*covmat[j][k];
		}
	}

	return sqrt(errf);
}

int fit_data(string datafile="data/ca40_1798_err.dat", Double_t Mom = 599.91, Double_t At = 40.078, Double_t Zt = 20.0){

	// this is just to record the execution time
	TStopwatch t;

	// string datapath = datafile;//"data/ca40_1798_err.dat";
	ifstream in(datafile.c_str(), ifstream::in );

	mom = Mom;
	at = At;
	zt = Zt;

	cout << "Data from file " << datafile.c_str() << endl;
	cout << "p = " << mom << " MeV/c \n"
		 << "A_t =  " << at << " amu \n"
		 << "Z_t = " << zt << "\n";

	Int_t i = 0;

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

// Set starting values and step sizes for parameters
   	static Double_t vstart[6] = {50.0, 50.0, 1.2, 1.2, 0.5, 0.5}; //if free
   	// static Double_t vstart[6] = {40.5, 111.0, 1.1, 1.1, 0.5, 0.5}; // if fixed
   	static Double_t step[6] = {0.05 , 0.05 , 0.0012 , 0.012, 0.0005, 0.0005};
   	
   	TF1* func = new TF1("func",csline,1.,60.,6);
   	func->SetParameters(vstart);
   	func->SetParNames("U0", "W0", "R0r", "R0i", "A0r", "A0i");

   	TVirtualFitter::SetDefaultFitter("Minuit");
   	TVirtualFitter *gMinuit = TVirtualFitter::Fitter(0,6);

   	for(int j = 0; j < 6; j++) gMinuit->SetParameter(j, func->GetParName(j), func->GetParameter(j), 0.001, 0,0);
   	gMinuit->SetFCN(fcn);
  
   	arglist[0] = 2;
   	gMinuit->ExecuteCommand("SET PRINT",arglist,1);

// Now ready for minimization step
    arglist[0] = 1500;
    arglist[1] = 1.;
    gMinuit->ExecuteCommand("MIGRAD", arglist,2);

// Print results
    Double_t minParams[6], parErrors[6];
    for(int j = 0; j < 6; j++)
    {
    	minParams[j] = gMinuit->GetParameter(j);
    	parErrors[j] = gMinuit->GetParError(j);
    }

    double chi2, edm, errdef;
  	int nvpar, nparx;
  	gMinuit->GetStats(chi2,edm,errdef,nvpar,nparx);

	TGraphErrors * g3 = new TGraphErrors("data/ca40_1798_err.dat","%lg %lg %lg");
	g3->SetTitle("#bar{p} + ^{40}Ca @ 179.8 MeV (p=599.91 MeV/c)");

	gStyle->SetOptFit(1);
	TCanvas * c = new TCanvas("c", "c", 800,800);

	func->SetParameters(minParams);
	func->SetParErrors(parErrors);
	func->SetChisquare(chi2);
	int ndf = x.size()-nvpar;
	func->SetNDF(ndf);

	TGraphErrors* gerr = new TGraphErrors(150);
	// cout << endl;
	// parErrors->as independent parameters!
	// Must implement the correct formula for
	// error propagation considering 
	// correlated parameters:
	// sigma^2 = sum(i,j) der(i)*der(j)*covmat[i][j] -> TO CHECK!
	for(int j = 0; j < 150; j++)
	{
		Double_t xx[1] = {(j+1)*(1./150.)*60.};
		Double_t ci[1];
		// Double_t err0 = func->GradientPar(0,xx,0.01)*sqrt(covmat[0][0]);
		// Double_t err1 = func->GradientPar(1,xx,0.01)*sqrt(covmat[1][1]);
		// Double_t err2 = func->GradientPar(2,xx,0.01)*sqrt(covmat[2][2]);
		// Double_t err3 = func->GradientPar(3,xx,0.01)*sqrt(covmat[3][3]);
		// Double_t err4 = func->GradientPar(4,xx,0.01)*sqrt(covmat[4][4]);
		// Double_t err5 = func->GradientPar(5,xx,0.01)*sqrt(covmat[5][5]);
		// Double_t errf = sqrt(err0*err0 + err1*err1 + err2*err2 + err3*err3 + err4*err4 + err5*err5);
		gerr->SetPoint(j,xx[0],func->Eval(xx[0]));
		// gMinuit->GetConfidenceIntervals(j,1,xx,ci);
		Double_t errf = error_prop(6,gMinuit,func,xx);
		gerr->SetPointError(j,0,errf);
	}

	gPad->SetLogy();
	c->cd();
	g3->SetMarkerStyle(20);
	g3->SetMarkerSize(0.5);
	g3->Draw("AP");
	func->Draw("SAME");
	
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	c1->cd();
	g3->Draw("AP");
	gerr->SetFillStyle(1001);
	gerr->SetFillColor(kRed);
	gerr->Draw("e3 SAME");
	// crs->Draw("SAME");
	t.Stop();
	t.Print("u");

	cout << '\a'; // sound at the end

	c->Print("output.pdf");
	c->SaveAs("results/ca40_1798.C");
	c1->SaveAs("prova.root");

	return 0;
}