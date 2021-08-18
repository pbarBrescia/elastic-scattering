#include "get_stdout.h"

Double_t mom, at, zt, charge;
string opt;

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

Double_t error_prop(Int_t npar, TF1 *func, const Double_t *xx)
{
	Double_t errf = 0.,sqerrf = 0.;
	for(int j = 0; j < npar; j++)
	{
		Double_t parerror = func->GetParError(j);
		errf += (func->GradientPar(j,xx))*parerror;
		sqerrf += errf*errf;
	}

	return sqrt(sqerrf);
}

void plot_band(string fitres, string data, Double_t Charge, Double_t Mom, Double_t At, Double_t Zt,string Opt)
{
	auto start = chrono::system_clock::now();
	time_t start_time = chrono::system_clock::to_time_t(start);
	cout << "Started at " << ctime(&start_time) << endl;

	charge = Charge;
	mom = Mom;
	at = At;
	zt = Zt;
	opt = Opt;

// open the file containing the TFitResultPtr
	TFile* ff = new TFile(fitres.c_str());
	TFitResult* frptr = (TFitResult*) ff->Get("fitresptr");
	Double_t xmin, xmax;
	if(opt=="ang") {xmin=1.; xmax=70.;}
	else if(opt=="mom") {xmin=mom; xmax=400.;}
	TF1* fitnew = new TF1("fitnew",csline,xmin,xmax,6);
	TGraphErrors* gdata = new TGraphErrors(data.c_str(),"%lg %lg %lg");
// Set pars and errors in TF1 from TFitResultsPtr
	fitnew->SetParameters(frptr->GetParams());
	fitnew->SetParErrors(frptr->GetErrors());

	Int_t npt = 200;
	TGraphErrors* gerr = new TGraphErrors(npt+1);
	Double_t q = 1.-0.025; // 1-alpha/2 for t-Student's distribution
	Double_t ndf = (gdata->GetN())-(frptr->NPar()); // degree of freedom (n. of data-par)
	Double_t t95 = ROOT::Math::tdistribution_quantile(q,ndf-1); // t-Student's distr

	// cout << xmin << " " << xmax << endl;
	for(int i=0; i<npt+1; i++)
	{
		Double_t errtot=0.;
		const Double_t xx[1] = {xmin+i*(1./npt)*(xmax-xmin)};
		// calculating the 95% Confidence Intervals
		Double_t err = t95*error_prop(fitnew->GetNpar(), fitnew, xx);
		gerr->SetPoint(i,xx[0],fitnew->Eval(xx[0]));
		if(i==0 || i%50 ==0) cout << xx[0] << " " << fitnew->Eval(xx[0]) << endl;
		gerr->SetPointError(i,0.,err);
	}

	TCanvas* c = new TCanvas("c","c",800,600);
	c->cd();
	gPad->SetLogy();
	gerr->SetTitle("Best fit + 95% CI");
	gerr->SetFillStyle(3003);
	gerr->SetFillColor(kBlue);
	gerr->SetLineColor(kBlack);
	gerr->Draw("AC E3");
	gdata->SetTitle("Exp. data");
	gdata->SetMarkerStyle(20);
	gdata->SetMarkerSize(0.5);
	gdata->Draw("PSAME");
	c->BuildLegend();

// Save root file and png of the plot in fig/ directory
	string proj;
	if (charge == -1.) proj = "antip";
	else if (charge == 0.) proj = "antin";
	gerr->SaveAs(Form("fig/plot_%.3f_%3.f_%s_%s.root",at,fitnew->GetParameter(0),opt.c_str(),proj.c_str()));
	c->SaveAs(Form("fig/plot_%.3f_%.3f_%s_%s.png",at,fitnew->GetParameter(0),opt.c_str(),proj.c_str()));

	auto end = chrono::system_clock::now();
	time_t end_time = chrono::system_clock::to_time_t(end);
	chrono::duration<double> delta_time = end-start;
	cout << "Ended at " << ctime(&end_time) << "after " << (int)(delta_time.count()/60.) << " min " << 60.*((delta_time.count()/60.)-(int)(delta_time.count()/60.))  << " s" << endl;

	cout << '\a'; // sound at the end
}