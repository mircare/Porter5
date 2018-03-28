
#ifndef BRNN_h
#define BRNN_h 1

#include <stdlib.h>
#include <math.h>
//#include "NNr.h"
#include "NN.h"
#include "NNt.h"

//#define MAX 8192

#define NTH 1



// BRNNr ver. 3.01 (5/11/2003)
// Copyright (C) Gianluca Pollastri 2003
//
// BRNNr with linear output
// inputs real valued or categorical
//
// Version 3.0:
// -compatible with NN,NNt,NNr 3.0
// -can be used as module in larger architecture (see 'backthrough'
// and the 'backp' parameter)
// -full shortcuts now operative
// 



class BRNN {
 public:

  int MAX;

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
  int Step;
  int shortcut;
  int doubleo;
  int modular;


  NN* NetOut;
//  NN* NetOut2;
  NNt* NetF;
  NNt* NetB;

  double** FF;
  double** BB;
  double** FFbp;
  double** BBbp;

  double* P_F;
  double* P_B;

  double* Y;
  double* BP;

  double error;
  double errorF;
  double errorB;

  double epsilon;

  void alloc(int the_MAX);
  void dealloc();



  BRNN(int NU, int NY, int NH,  int context ,int Moore,
	  int NF, int NB, int NH2, int CoF, int CoB, int Step,
	  int shortcut, int doubleo=0);
  BRNN(istream& is);
  BRNN(BRNN* from, int seq_length);
  ~BRNN();
  void read(istream& is);
  void write(ostream& os);

  void resetGradient();
  void initWeights(int seed);


  void copy_dW(BRNN* from) {
	NetOut->copy_dW(from->NetOut);
	NetF->copy_dW(from->NetF);
	NetB->copy_dW(from->NetB);
	error += from->error;
	errorF += from->errorF;
	errorB += from->errorB;
  }

  void F1_F(double* seq, int t, int length);
  void F1_F(int* seq, int t, int length);
  void B1_B(double* seq, int t, int length);
  void B1_B(int* seq, int t, int length);
  void propagate(double* seq, int length);
  void propagate(int* seq, int length);

  void forward(double* seq, int t, int length);
  void forward(int* seq, int t, int length);

  void forward(double* seq, int t1, int t2, int length);
  void forward(int* seq, int t1, int t2, int length);

  void forward_backward(double* seq, int t, int length, double* target, int backp=0);
  void forward_backward(int* seq, int t, int length, double* target, int backp=0);

  void forward_backward(double* seq, int t1, int t2, int length, double** target, int backp=0);
  void forward_backward(int* seq, int t1, int t2, int length, double** target, int backp=0);



  void F1_Fbp(double* seq, int t, int length, int backp=0);
  void F1_Fbp(int* seq, int t, int length, int backp=0);
  void B1_Bbp(double* seq, int t, int length, int backp=0);
  void B1_Bbp(int* seq, int t, int length, int backp=0);
  void back_propagate(double* seq, int length, int backp=0);
  void back_propagate(int* seq, int length, int backp=0);


  void extimation(double* seq, int* y, int length, int backp=0);
  void extimation(int* seq, int* y, int length, int backp=0);

  void extimation(double* seq, double* y, int length, int backp=0);
  void extimation(int* seq, double* y, int length, int backp=0);

  void backthrough(double* seq, double* y, int length, int backp=0) {
    NetOut->set_output(0);
    extimation(seq,y,length,backp);
    NetOut->set_output(1);
  }
  void backthrough(int* seq, double* y, int length, int backp=0) {
    NetOut->set_output(0);
    extimation(seq,y,length,backp);
    NetOut->set_output(1);
  }


  void resetBP(int length);


  void maximization();
  void maximizationL1();
  void maximizationClipped();

  void Feed(double* seq, int length);
  void Feed(int* seq, int length);
  void predict(double* seq, int length);
  void predict(int* seq, int length);


  double* out() {return Y;}
  double* back_out() {return BP;}




  double getError() {
	return error;
	};
  double getErrorF() {
	return errorF;
	};
  double getErrorB() {
	return errorB;
	};
  void resetError() {
	error=0.0;
	errorF=0.0;
	errorB=0.0;
	};

  void setEpsilon(double eps) { epsilon=eps; };


};


#endif // BRNN_h
