#include <iostream> //libreria per la stampa su schermo e lettura da tastiera
#include <fstream> // libreria per la stampa o lettura su/da file
#include <cmath>
#include <cstdlib>
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

complex<double> drop_imag(complex<double> z)
{
  double epsilon = 1e-07;

  if (abs(imag(z)) <= epsilon) z = real(z);

  return z;
}


// gamma function with Lanczos approximation
// re-adapting the python implementation: 
// https://en.wikipedia.org/wiki/Lanczos_approximation
complex<double> gamma(complex<double> z)
{
   vector<double> p {676.5203681218851
    ,-1259.1392167224028
    ,771.32342877765313
    ,-176.61502916214059
    ,12.507343278686905
    ,-0.13857109526572012
    ,9.9843695780195716e-6
    ,1.5056327351493116e-7};

  complex <double> eye  (0.,1.);
  complex <double> uno  (1.,0.);
  complex<double> y,x,t;
  if (real(z) < 0.5)
  {
    y = M_PI/(sin(M_PI*z)*gamma(uno-z));
  }
  else
  {
    z -= 1;
    x = 0.99999999999980993;
    for (int i=0; i < p.size(); i++)
    {
      complex<double> j (i,0);
      x += p.at(i) / (z + j + uno);
    }
    t = z + (complex<double>)p.size() - 0.5;
    y = sqrt(2 * M_PI) * pow(t,(z + 0.5)) * exp(-t) * x;
    return drop_imag(y);
  }
}

int main(int argc, char* argv[]){



  // Code added in 2021 to use command line arguments
    if (argc < 2) {
        std::cerr << "Usage:   " << argv[0] << " [input data file] [pbar lab momentum] [target mass] [target charge] [proj charge]" << std::endl;
        std::cerr << "Example: " << argv[0] << " input.dat 50.0 40.078 20.0 -1.0" << std::endl;
        return 1;
    }

  string input_data_file(argv[1]);
  // cout << input_data_file.c_str() << endl;

  double pbar_lab_p    = atof(argv[2]);
  double target_mass   = atof(argv[3]);
  double target_charge = atof(argv[4]);
  double proj_charge = atof(argv[5]);

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
  double zp = proj_charge;
  // double zp = -1.0;
  double ak = pmc/197.327; // k in 1/fm, old value: 197.32
  //cm effects removed for coherence with the other program
  //  double mcorrection = (mass*1.)/(mass+1.);
  //double lambda = -6.85*zeta*mcorrection/pmc;
  double mproj;
  if(zp == -1.0) mproj = 938.27;
  else if (zp == 0.0) mproj = 939.565;
  double mu = mass*mproj/(mproj+mass); // MeV
  double esq = 1.4399764; // MeV*fm
  double hxc = 197.327; // MeV*fm
  //double lambda = -6.85*zeta/pmc;
  double lambda = zp*esq*mu*zeta/(hxc*pmc);

  //cout << lambda << " lambda\n";

  double Rz = 1;
  double Iz = lambda;

  complex<double> z (Rz,Iz);

  complex <double> eye  (0.,1.);
  complex <double> uno  (1.,0.);
  
  complex<double> amplitude1 = gamma(z);  
  complex<double> amplitude2 = gamma(conj(z));
  complex<double> amplitude3 = amplitude1/amplitude2;
 
  //cout << amplitude1 << " " << amplitude2 << " " << amplitude3 << "\n";

  int nstep=500;
  double step = (double) 1./nstep;
  // cout << step << endl;

  for (int itheta = 0; itheta<nstep; itheta++){

    double theta = (itheta)*step*M_PI;

    double sin12 = sin(theta/2.0);
    double sin122 = sin12*sin12;
    double lsin12 = log(sin12);

    double earg = -2.*lambda*lsin12;
    complex <double> cearg = eye*earg;
    complex <double> cesp = exp(cearg);

    //cout << earg << " " << eye << " " 
    //  << cearg << " " << abs(cesp)<< " " << abs(amplitude3) << "\n";

    //coulomb amplitude
    complex <double> amplitude = -(amplitude3 * cesp * (0.5*lambda/ak) )/ sin122;

    double aaa = abs(amplitude);

    //coulomb differential sigma (fm^2/sr)
    double dsigmaomega = aaa*aaa;

    double nucleartheta = 0, rnuclear = 0, inuclear = 0; //, rtot=0, itot=0;
    in1 >> nucleartheta >> rnuclear >> inuclear; //>> rtot >> itot;
    // cout << nucleartheta << rnuclear << inuclear << endl;

    double dsigmanuclear = rnuclear*rnuclear + inuclear*inuclear;

    double atotr = real(amplitude) + rnuclear;
    double atoti = imag(amplitude) + inuclear;

    double dsigmatot = atotr*atotr + atoti*atoti;

    if (theta>0)
    cout << theta*(180./M_PI) << " " << nucleartheta*(180./M_PI) 
	 << " " << 10.*dsigmanuclear << " " 
	 << 10.*dsigmaomega << " " << 10.*dsigmatot << "\n";
   // 10* -> conversion to millibarn
  }

  in1.close();

  return 0;
}
