void conv(string conversion = "kin_to_mom", double k = 50, double m = 938.27)
{
	if(conversion == "kin_to_mom")
	{
		// k = kinetic energy in MeV
		// m = mass in MeV/c^2
		double p = sqrt(k*(k+2*m));
		// p = momentum in MeV/c
		cout << "Kinetic energy: " << k << " MeV \n"
			 << "Particle mass: " << m << " MeV/c^2 \n"
			 << "Momentum: " << p << " MeV/c \n"; 
	}
	else if(conversion == "mom_to_kin")
	{
		// k = momentum in MeV/c
		// m = mass in MeV/c^2
		double kin = m*(sqrt(1+(k/m)*(k/m))-1);
		// p = momentum in MeV/c
		cout << "Kinetic energy: " << kin << " MeV \n"
			 << "Particle mass: " << m << " MeV/c^2 \n"
			 << "Momentum: " << k << " MeV/c \n";
	}
	else
	{
		cerr << "ERROR: please, provide 3 parameters: \n"
			 << "root -l 'conv.cc(<conversion>, <momentum/kinetic energy>, <mass>)' \n"
			 << "<conversion> : \"kin_to_mom\" or \"mom_to_kin\" \n";
		exit(1);
	}

}