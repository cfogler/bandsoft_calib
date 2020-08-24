#ifndef __CALIB_HELPER_H__
#define __CALIB_HELPER_H__

#include <fstream>
#include <iostream>
#include <string>


using namespace std;

class shiftsReader {
	public:
		void LoadInitBar( string filename );
		void LoadInitBarFadc( string filename );
		void LoadInitRun( string filename );
		void LoadInitRunFadc( string filename );
		double * getInitBar(void);
		double * getInitBarFadc(void);
		double * getInitRun(void);
		double * getInitRunFadc(void);
	private:
		double InitBar[600] = {0.};
		double InitBarFadc[600] = {0.};
		double InitRun[100000] = {0.};
		double InitRunFadc[100000] = {0.};
};


const int maxNeutrons	= 15;
int getRunNumber( string filename );

#endif
