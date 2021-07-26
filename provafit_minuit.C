#include "get_stdout.h"

using namespace std;

Double_t x[44], data[44], error[44];

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

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   const Int_t nbins = 44;
   Int_t i;

//calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   for (i=0;i<nbins; i++) {
     delta  = (data[i]-csline(x[i],par))/error[i];
     chisq += delta*delta;
   }
   f = chisq;
}


void provafit_minuit(){

	// this is just to record the execution time
	TStopwatch t;

	string datapath = "data/ca40_1798_err.dat";
	ifstream in(datapath.c_str(), ifstream::in );

	Int_t i = 0;

	while(!in.eof())
	{
		in >> x[i] >> data[i] >> error[i];
		i++;
	}

	in.close();

	// for(i=0; i<44; i++) cout << x[i] << " " << data[i] << " " << error[i] << endl;

	TMinuit *gMinuit = new TMinuit(6);  //initialize TMinuit with a maximum of 6 params
   	gMinuit->SetFCN(fcn);

   	Double_t arglist[10];
   	Int_t ierflg = 0;

   	arglist[0] = 1;
   	t.Start();
   	gMinuit->mnexcm("SET ERR", arglist , 1, ierflg);

// Set starting values and step sizes for parameters
   	static Double_t vstart[6] = {10.0, 50.0, 1.0, 1.0, 0.5, 0.5}; //if free
   	// static Double_t vstart[6] = {40.5, 111.0, 1.1, 1.1, 0.5, 0.5}; // if fixed
   	static Double_t step[6] = {0.01 , 0.01 , 0.005 , 0.005, 0.005, 0.005};
   	gMinuit->mnparm(0, "U0", vstart[0], step[0], 0.0, 500.,ierflg);
   	gMinuit->mnparm(1, "W0", vstart[1], step[1], 0.0, 500.,ierflg);
   	gMinuit->mnparm(2, "R0r", vstart[2], step[2], 0.8, 2.0,ierflg);
   	gMinuit->mnparm(3, "R0i", vstart[3], step[3], 0.8, 2.0,ierflg);
   	gMinuit->mnparm(4, "A0r", vstart[4], step[4], 0.3, 1.0,ierflg);
   	gMinuit->mnparm(5, "A0i", vstart[5], step[5], 0.3, 1.0,ierflg);
   	// gMinuit->FixParameter(0);
   	// gMinuit->FixParameter(1);
   	// gMinuit->FixParameter(2);
   	// gMinuit->FixParameter(3);

// Now ready for minimization step
    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

// Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    gMinuit->mnprin(3,amin);

	TGraphErrors * g3 = new TGraphErrors("data/ca40_1798_err.dat","%lg %lg %lg");
	g3->SetTitle("#bar{p} + ^{40}Ca @ 179.8 MeV (p=599.91 MeV/c)");

	gStyle->SetOptFit(1);
	TCanvas * c = new TCanvas("c", "c", 800,800);
	
	// TF1 *crs = new TF1("crs",csline,5,55,6);
	Double_t u0,erru0,w0,errw0;
	Double_t r0r,errr0r,r0i,errr0i;
	Double_t a0r,erra0r,a0i,erra0i;
	gMinuit->GetParameter(0,u0,erru0);
	gMinuit->GetParameter(1,w0,errw0);
	gMinuit->GetParameter(2,r0r,errr0r);
	gMinuit->GetParameter(3,r0i,errr0i);
	gMinuit->GetParameter(4,a0r,erra0r);
	gMinuit->GetParameter(5,a0i,erra0i);
	Double_t parm[6] = {u0,w0,r0r,r0i,a0r,a0i};
	TGraphErrors* grfit = new TGraphErrors(100);
	// TGraph* grfit2 = dynamic_cast<TGraph*>(gMinuit->GetPlot());
	Double_t X = 0;
	Double_t value;
	for(int j = 1; j<101; j++)
	{
		X = j*0.01*60.; 
		value = csline(X,parm);
		grfit->SetPoint(j,X,value);
	}
	
	gPad->SetLogy();
	c->cd();
	g3->SetMarkerStyle(20);
	g3->SetMarkerSize(0.5);
	g3->Draw("AP");

	grfit->Draw("SAME");
	// crs->Draw("SAME");
	t.Stop();
	t.Print("u");

	//c->Print("output.pdf");
}