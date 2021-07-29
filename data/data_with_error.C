using namespace std;

void data_with_error(string path_datafile)
{
	ifstream in(path_datafile.c_str(),ifstream::in);
	vector<double> theta, dsigma;
	double a, b;
	int size = 0;
	cerr << "Reading file " << path_datafile.c_str() << " \n";
	while(!in.eof())
	{
		in >> a >> b;
		theta.push_back(a);
		dsigma.push_back(b);
		// cout<<a<<" "<<b<<endl;
	}

	in.close();

	string outdata = path_datafile.c_str();
	string::iterator it = outdata.end();
	string add = "_err";
	outdata.insert(path_datafile.size()-4,add);

	ofstream fout(outdata.c_str(), ofstream::out);
	cerr << "Saving new data with error in " << outdata.c_str() << " \n";
	for(int i = 0; i<theta.size()-1; i++)
	{
		if((i+2)%2==0)
		{
			double err=abs(dsigma.at(i+1)-dsigma.at(i));
			fout << theta.at(i) << " " << dsigma.at(i) << " " << err << endl;
		}
	}
	fout.close();
}