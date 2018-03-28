
#ifndef MModel_h
#define MModel_h 1

#include <stdlib.h>
#include <math.h>
#include "Sequence.h"
#include "BRNN.h"


class MModel {
 private:
  int nModels;

  int* NU;
  int* NY;
  int* NH;

  int* context;
  int* Moore;

  int* NF;
  int* NB;
  int* NH2;

  int* CoF;
  int* CoB;
  int* Cseg;
  int* Cwin;
  int* shortcut;
  int* Step;

  double** Thresholds;

  int cycles;
  double* dcycles;

  int* modular;


  BRNN** Net;
  BRNN** NetF;

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

  MModel(istream& is);

  void predict(Sequence* seq);
  void predict(Sequence* seq, int cy);

  void predict_probs(Sequence* seq);
  void predict_probs(Sequence* seq, int cy);


  void predict(Sequence* seq, double threshold);



//  double* out() {return NetF->out();}
  int** getConf() {return Conf;}

  int getNErrors() { return nerrors;};

  int getNErrors_(int i) { return nerrors_[i];};
  int getClasses() { return NY[0];};

  int* getCounted() {return counted;}


  double* getdcycles() {return dcycles;}

   void resetNErrors() {
	error=0;
	nerrors=0;
	memset(nerrors_,0,NY[0]*sizeof(int));
	memset(counted,0,NY[0]*sizeof(int));
	for (int p=0;p<NY[0];p++)
	  for (int y=0;y<NY[0];y++)
		Conf[p][y]=0;
	  for (int c=0;c<cycles;c++) {
		  dcycles[c]=0;
	  }
	for (int n=0;n<nModels;n++) {
		Net[n]->resetError();
		NetF[n]->resetError();
	}
	};

  double get_error() {
	return error;
	};

  void reset_squared_error() {
	for (int n=0;n<nModels;n++) {
		Net[n]->resetError();
		NetF[n]->resetError();
	}
	  for (int c=0;c<cycles;c++) {
		  dcycles[c]=0;
	  }
	};



};


#endif // MModel_h
