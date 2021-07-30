#include "get_stdout.h"
#include <TMinuit.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TStyle.h>

#include <vector>
#include <cstdlib>
#include <fstream>


using namespace std;

// Int_t ndata;
vector<Double_t> x;// = new Double_t[ndata];
vector<Double_t> data;// = new Double_t[ndata]; 
vector<Double_t> error;// = new Double_t[ndata];

Double_t csline(Double_t* x,Double_t *par) {
	// run an external command and get its output in a TString
	// TString s = get_stdout(Form("./p_scan.sh -i 599.91 40.078 20.0 40.5 111.0 1.1 1.1 1.25 0.63 %f ang %f",par[0]/*, par[1], par[2], par[3], par[4], par[5]*/,x[0]));
	TString s = get_stdout(Form("./bin/antip_scan 599.91 40.078 20.0 %f %f %f %f 1.25 %f %f ang %f",par[0],par[1],par[2],par[3],par[4],par[5],x[0]));
	// split the TString using space as delimiter
	// Example taken from here: https://root-forum.cern.ch/t/split-tstring-by-delimeter-in-root-c/18228/2
	TObjArray *ts = s.Tokenize(" ");
	// get the 1st element in a Double_t and returns to the program
	Double_t ret = atof( ((TObjString *) (ts->At(1)) )->String());
	return ret;
}

Double_t csline(Double_t x,Double_t *par) {
	// run an external command and get its output in a TString
	// TString s = get_stdout(Form("./p_scan.sh -i 599.91 40.078 20.0 40.5 111.0 1.1 1.1 1.25 0.63 %f ang %f",par[0]/*, par[1], par[2], par[3], par[4], par[5]*/,x[0]));
	TString s = get_stdout(Form("./bin/antip_scan 599.91 40.078 20.0 %f %f %f %f 1.25 %f %f ang %f",par[0],par[1],par[2],par[3],par[4],par[5],x));
	// split the TString using space as delimiter
	// Example taken from here: https://root-forum.cern.ch/t/split-tstring-by-delimeter-in-root-c/18228/2
	TObjArray *ts = s.Tokenize(" ");
	// get the 1st element in a Double_t and returns to the program
	Double_t ret = atof( ((TObjString *) (ts->At(1)) )->String());
	return ret;
}

Double_t derivative(Double_t x, Double_t *par, Int_t npar)
{
	Double_t h[6] = {0.};
	h[npar] = 1e-4;
	Double_t parplus2[6], parplus[6], parminus[6], parminus2[6];
	for(int j=0; j<6; j++)
	{
		parplus[j] = par[j]+h[j];
		parplus2[j] = par[j]+2*h[j];
		parminus[j] = par[j]-h[j];
		parminus2[j] = par[j]-2*h[j]; 
	}
	Double_t D = (-csline(x,parplus2)+8*csline(x,parplus)-8*csline(x,parminus)+csline(x,parminus2))/(12*h[npar]);
	return D;
}

Double_t derivative(Double_t x, Double_t *par)
{
	Double_t h = 1e-4;
	Double_t D = (-csline(x+2*h,par)+8*csline(x+h,par)-8*csline(x-h,par)+csline(x-2*h,par))/(12*h);
	return D;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   Int_t i;

//calculate chisquare
   Double_t chisq = 0;
   Double_t delta,den1,den2;
   Int_t ndata = (unsigned) x.size();
   for (i=0; i<ndata; i++) {
   	 den1 = error[i];
   	 den2 = derivative(x[i],par)*0.2;//(data[i]/csline(x[i],par))*0.2;
     delta  = (data[i]-csline(x[i],par));// /error[i];
     // chisq += delta*delta;
     chisq += delta*delta/(den1*den1+den2*den2);
   }
   f = chisq;
}

int provafit_minuit(){

	// this is just to record the execution time
	TStopwatch t;

	string datapath = "data/ca40_1798_err.dat";
	ifstream in(datapath.c_str(), ifstream::in );

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

	// for(i=0; i<44; i++) cout << x[i] << " " << data[i] << " " << error[i] << endl;

	// TMinuit *gMinuit = new TMinuit(6);  //initialize TMinuit with a maximum of 6 params
   	// gMinuit->SetFCN(fcn);
   	Double_t arglist[10];
   	// Int_t ierflg = 0;

   	// arglist[0] = 1;
   	// t.Start();
   	// gMinuit->mnexcm("SET ERR", arglist , 1, ierflg);

// Set starting values and step sizes for parameters
   	static Double_t vstart[6] = {100.0, 100.0, 1.0, 1.0, 0.25, 0.25}; //if free
   	// static Double_t vstart[6] = {40.5, 111.0, 1.1, 1.1, 0.5, 0.5}; // if fixed
   	static Double_t step[6] = {0.5 , 0.5 , 0.015 , 0.015, 0.005, 0.005};
   	
   	TF1* func = new TF1("func",csline,1.,60.,6);
   	func->SetParameters(vstart);
   	func->SetParNames("U0", "W0", "R0r", "R0i", "A0r", "A0i");

   	TVirtualFitter::SetDefaultFitter("Minuit");
   	TVirtualFitter *gMinuit = TVirtualFitter::Fitter(0,6);

   	for(int j = 0; j < 6; j++) gMinuit->SetParameter(j, func->GetParName(j), func->GetParameter(j), 0.01, 0,0);
   	gMinuit->SetFCN(fcn);
   	
   	// gMinuit->SetErrorDef(1.);
   	// gMinuit->mnparm(0, "U0", vstart[0], step[0], 0, 0,ierflg);
   	// gMinuit->mnparm(1, "W0", vstart[1], step[1], 0, 0,ierflg);
   	// gMinuit->mnparm(2, "R0r", vstart[2], step[2], 0, 0,ierflg);
   	// gMinuit->mnparm(3, "R0i", vstart[3], step[3], 0, 0,ierflg);
   	// gMinuit->mnparm(4, "A0r", vstart[4], step[4], 0, 0,ierflg);
   	// gMinuit->mnparm(5, "A0i", vstart[5], step[5], 0, 0,ierflg);
   	// gMinuit->FixParameter(0);
   	// gMinuit->FixParameter(1);
   	// gMinuit->FixParameter(2);
   	// gMinuit->FixParameter(3);
   	arglist[0] = 2;
   	gMinuit->ExecuteCommand("SET PRINT",arglist,1);

// Now ready for minimization step
    arglist[0] = 50;
    arglist[1] = 1.;
    // gMinuit->mnexcm("SIM", arglist, 1, ierflg);
    // gMinuit->mnimpr();
    // gMinuit->mnexcm("MIG", arglist, 2, ierflg);
    // gMinuit->SetMaxIterations(100);
    // gMinuit->mnexcm("MINOS",arglist,1,ierflg);
    gMinuit->ExecuteCommand("MIGRAD", arglist,2);

// Print results
    // Double_t amin,edm,errdef;
    // Int_t nvpar,nparx,icstat;
    // gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    // gMinuit->mnprin(3,amin);
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
	
	// // TF1 *crs = new TF1("crs",csline,5,55,6);
	// Double_t u0,erru0,w0,errw0;
	// Double_t r0r,errr0r,r0i,errr0i;
	// Double_t a0r,erra0r,a0i,erra0i;
	// gMinuit->GetParameter(0,u0,erru0);
	// gMinuit->GetParameter(1,w0,errw0);
	// gMinuit->GetParameter(2,r0r,errr0r);
	// gMinuit->GetParameter(3,r0i,errr0i);
	// gMinuit->GetParameter(4,a0r,erra0r);
	// gMinuit->GetParameter(5,a0i,erra0i);
	// Double_t parm[6] = {u0,w0,r0r,r0i,a0r,a0i};
	// // cout << u0 << " "
	// // 	<< w0 << " "
	// // 	<< r0r << " "
	// // 	<< r0i << " "
	// // 	<< a0r << " "
	// // 	<< a0i << " " << endl;

	// TGraphErrors* grfit = new TGraphErrors(150);
	// // TGraphErrors* gerr = new TGraphErrors(100);

	// Double_t X = 0;
	// Double_t value;
	// for(int j = 0; j<150; j++)
	// {
	// 	X = (j+1)*(1./150.)*60.; 
	// 	value = csline(X,parm);
	// 	grfit->SetPoint(j,X,value);
	// 	Double_t err0 = derivative(X,parm,0)*erru0;
	// 	Double_t err1 = derivative(X,parm,1)*errw0;
	// 	Double_t err2 = derivative(X,parm,2)*errr0r;
	// 	Double_t err3 = derivative(X,parm,3)*errr0i;
	// 	Double_t err4 = derivative(X,parm,4)*erra0r;
	// 	Double_t err5 = derivative(X,parm,5)*erra0i;
	// 	Double_t errx = derivative(X,parm)*0.2;
	// 	Double_t errf = sqrt(err0*err0+err1*err1+err2*err2+err3*err3+err4*err4+err5*err5);
	// 	grfit->SetPointError(j,0,errf);
	// 	// gerr->SetPoint(j-1,X,0);
	// }
	// (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gerr,0.95);

	gPad->SetLogy();
	c->cd();
	g3->SetMarkerStyle(20);
	g3->SetMarkerSize(0.5);
	g3->Draw("AP");
	// grfit->SetFillColor(kRed);
	// grfit->SetFillStyle(3002);
	// grfit->Draw("e3 SAME");
	func->DrawCopy("cont1 same");
	// gerr->Draw("e3 SAME");
	// crs->Draw("SAME");
	t.Stop();
	t.Print("u");

	cout << '\a'; // sound at the end

	//c->Print("output.pdf");

	return 0;
}