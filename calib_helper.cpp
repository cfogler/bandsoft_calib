#include "calib_helper.h"


int getRunNumber( string filename ){
	string parsed = filename.substr( filename.find("noshift_") );
	string moreparse = parsed.substr(8,6);
	cout << "\t*Intepreted run number from file name: " << stoi(moreparse) << "\n";
        return stoi(moreparse);
}


void shiftsReader::LoadInitBar( string filename ){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitBar[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadInitBarFadc( string filename ){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitBarFadc[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadInitRun( string filename ){
	ifstream f;
	int runno;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> runno;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitRun[runno] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadInitRunFadc( string filename ){
	ifstream f;
	int runno;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> runno;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitRunFadc[runno] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}

double * shiftsReader::getInitBar(void){
	return InitBar;
}
double * shiftsReader::getInitBarFadc(void){
	return InitBarFadc;
}
double * shiftsReader::getInitRun(void){
	return InitRun;
}
double * shiftsReader::getInitRunFadc(void){
	return InitRunFadc;
}
