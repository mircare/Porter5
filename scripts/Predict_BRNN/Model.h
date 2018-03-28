
#ifndef Model_h
#define Model_h 1

#include <stdlib.h>
#include <math.h>
#include "Sequence.h"
#include "BRNN.h"


class Model {
 private:
  int NU;
  int NY;
  int NH;

  int context;
  int Moore; 

  int NF;
  int NB;
  int NH2;

  int CoF;
  int CoB;
  int Cseg;
  int Cwin;
  int shortcut;
  int Step;

  double* Thresholds;

  int cycles;
  double* dcycles;

  int modular;


  BRNN* Net;
  BRNN* NetF;

  int** Conf;

//  double temp_error;
  int temp_aas;
  
  int* counted;
  double squared_error;
  double error;
  int nerrors;
  int* nerrors_;

  double epsilon;

  void alloc();


 public:

  Model(int NU, int NY, int NH,  int context ,int Moore, int NF, int NB, int NH2, int CoF, 
        int CoB, int Cseg, int Cwin, int Step, int shortcut, double* Thresholds, int cycles=1);
  Model(istream& is);
  void read(istream& is);
  void write(ostream& os);

  void randomize(int seed);


  void extimation(Sequence* seq);
  void maximization();
  void maximizationL1();
  void maximizationClipped();

  void predict(Sequence* seq);
  void predict(Sequence* seq, int cy);
//  void predict(Sequence* seq, int W);
  double* out() {return NetF->out();}
  int** getConf() {return Conf;}

  int getNErrors() { return nerrors;};

  int getNErrors_(int i) { return nerrors_[i];};
  int getClasses() { return NY;};

  int* getCounted() {return counted;}


  double* getdcycles() {return dcycles;}

   void resetNErrors() { 
	error=0;
	nerrors=0;
	memset(nerrors_,0,NY*sizeof(int));
	memset(counted,0,NY*sizeof(int));
	for (int p=0;p<NY;p++)
	  for (int y=0;y<NY;y++)
		Conf[p][y]=0;
	  for (int c=0;c<cycles;c++) {
		  dcycles[c]=0;
	  }
	Net->resetError();
	NetF->resetError();
	};

  double get_error() { 
	return error;
	};
  double get_squared_error() { 
	return Net->getError();
	};
  double get_squared_errorf() { 
	return NetF->getError();
	};
  double get_squared_errorF() { 
	return Net->getErrorF();
	};
  double get_squared_errorB() { 
	return Net->getErrorB();
	};
  void reset_squared_error() { 
	Net->resetError();
	NetF->resetError();
	  for (int c=0;c<cycles;c++) {
		  dcycles[c]=0;
	  }
	};

  void setEpsilon(double eps) { 
	  Net->setEpsilon(eps);
	  NetF->setEpsilon(eps);
  };


};


#endif // Model_h
