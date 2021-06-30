#include <iostream> //libreria per la stampa su schermo e lettura da tastiera
#include <fstream> // libreria per la stampa o lettura su/da file
#include <cmath>
#include<cstdlib>
#include <cstdio> //input/output in stile C
#include <cstring>
#include <ctime>
#include <list>
#include <queue> 
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <typeinfo>
#include <string>
#include <numeric>
#include <functional>
#include <sstream>
#include <complex>

using namespace std;

int main(int argc, char* argv[]){



  // Code added in 2021 to use command line arguments
    if (argc < 2) {
        std::cerr << "Usage:   " << argv[0] << " [input data file] [pbar lab momentum] [target mass] [target charge]" << std::endl;
        std::cerr << "Example: " << argv[0] << " input.dat 50.0 40.078 20.0" << std::endl;
        return 1;
    }

  string input_data_file(argv[1]);

  double pbar_lab_p    = atof(argv[2]);
  double target_mass   = atof(argv[3]);
  double target_charge = atof(argv[4]);

  //attenzione: l'unità di misura di 1/ak è la stessa dell'ampiezza finale
  //lambda è adimensionale

  //A calculated cross section is in the same units as (1/ak)^2. 
  //Presently dim(pmc) = MeV/c, and dim(ak) = 1/fm -> dim sigma = fm^2, 
  //so the final sigma is multiplied by 10 before printing to have 
  //sigma in mb (10 mb = 1 fm^2). Dsigma/DOmega is mb/sr
  ifstream in1(input_data_file);

  double pmc = pbar_lab_p ;// old value: 50; //p in mev/c
  double zeta = target_charge ; // old value: 20.;
  double mass = target_mass * 931.494; // old value: 40.078 , 931.78 = u.a.m.
  double ak = pmc/197.33; // k in 1/fm, old value: 197.32
  //cm effects removed for coherence with the other program
  //  double mcorrection = (mass*1.)/(mass+1.);
  //double lambda = -6.85*zeta*mcorrection/pmc;
  double mu = mass*938.27/(938.27+mass); // MeV
  double esq = 1.4399764; // MeV*fm
  double hxc = 197.327; // MeV*fm
  //double lambda = -6.85*zeta/pmc;
  double lambda = - esq*mu*zeta/(hxc*pmc);

  //cout << lambda << " lambda\n";

  double Rz = 1;
  double Iz = lambda;

  int N1 = 200, N2 = N1/5;
  double du = 20./N1;

  double Rsum = 0, Isum = 0;
  for (int i=-N1; i< N2; i++){
    double u = i*du;
    double ef = 0;
    if (u < 4) {
      double aarg = u*Rz - exp(u);
      ef = exp(aarg);
    }
    else {
      ef = 0;
    }
    double RRef = cos(Iz)*ef;
    double IIef = sin(Iz)*ef;
    Rsum += (RRef*du);
    Isum += (IIef*du);
    //cout << i << " " << u << " " << ef << " " << Rsum << " " << Isum << " "  
    // << exp(u) << "\n";
  }

  complex <double> eye  (0.,1.);
  complex <double> uno  (1.,0.);
  
  complex<double> amplitude1 (Rsum,Isum);  
  complex<double> amplitude2 = conj(amplitude1);
  complex<double> amplitude3 = amplitude1/amplitude2;
 
  cout << amplitude1 << " " << amplitude2 << " " << amplitude3 << "\n";

  int nstep = 200;
  double step = (double) 1./nstep;

  for (int itheta = 0; itheta<nstep; itheta++){

    double theta = (itheta)*step*M_PI;

    double sin12 = sin(theta/2);
    double sin122 = sin12*sin12;
    double lsin12 = log(sin12);

    double earg = -2.*lambda*lsin12;
    complex <double> cearg = eye*earg;
    complex <double> cesp = exp(cearg);

    //cout << earg << " " << eye << " " 
    //  << cearg << " " << abs(cesp)<< " " << abs(amplitude3) << "\n";

    //coulomb amplitude
    complex <double> amplitude = (amplitude3 * cesp * (0.5*lambda/ak) )/ sin122;

    double aaa = abs(amplitude);

    //coulomb differential sigma (fm^2/sr)
    double dsigmaomega = aaa*aaa;

    double nucleartheta=0, rnuclear=0, inuclear=0;
    in1 >> nucleartheta >> rnuclear >> inuclear;

    double dsigmanuclear = rnuclear*rnuclear + inuclear*inuclear;

    double atotr = real(amplitude) + rnuclear;
    double atoti = imag(amplitude) + inuclear;

    double dsigmatot = atotr*atotr + atoti*atoti;

    if (theta>0)
    cout << theta*(180./M_PI) << " " << nucleartheta*(180./M_PI) 
	 << " " << 10*dsigmanuclear << " " 
	 << 10*dsigmaomega << " " << 10*dsigmatot << "\n";
   // 10* -> conversion to millibarn
  }

  in1.close();

  return 0;
}
