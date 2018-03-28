
#ifndef Layer_h
#define Layer_h 1
#include <math.h>
#include <iostream>
#include <string.h>
using namespace std;




// Layer ver 3.03
// 12/12/2003
// Copyright (C) Gianluca Pollastri 2003
//
// ANN Layers
// Linear, tanh and softmax outputs
// Categorical (one-hot), real-valued and mixed inputs.
//
//
// In version 3.0:
// -fixed 'overflow problem' that output error=0
//  for saturated softmax units
// -all (but the last one..) compatibility issues fixed
// -added updateWeightsClipped
//
// In version 3.01
// -fixed all the versions of gradient();
// -NUtot made an attribute
//
// In version 3.02:
// -fixed Linux bug in initWeights
// -added updateWeightsL1
//
// In version 3.03:
// -gradient fixed for softmax: (y-t)x, no f1


class Layer
{
public:
int NY;
int NU;
int NUr;
int* NK;

double* Y;
double* A;
double* U;	 //NU*NK
double* delta;    //NY
double* backprop; //NU*NK

double* W;
double* dW;
double* d2W;

double* B;        //NY
double* dB;       //NY
double* d2B;       //NY

int output; 	//0=no,1=yes
int ninput;		//0=input layer,1=just real side backprop,2=full backprop

int NUtot,NUplain;


void alloc(int NY, int NU, int* NK);


//public:

void softmax();
void squash();

Layer(Layer* from) {
int y;
	NY = from->NY;
	NU = from->NU;
	NUr = from->NUr;

NK=new int[NU+NUr];
for (int i=0; i<NU; i++)
	NK[i]=from->NK[i];
for (int i=NU; i<NU+NUr; i++)
        NK[i]=1;
alloc(NY,NU+NUr,NK);
ninput=0;
output=0;

NUplain =0;
for (int u=0; u<NU; u++) {
  NUplain += NK[u];
}

for (y=0; y<NY; y++) {
	B[y] = from->B[y];
	dB[y] = 0;
	d2B[y] = 0;
        for (int u=0; u<NUtot; u++) {
			W[y*NUtot+u] = from->W[y*NUtot+u];
			dW[y*NUtot+u] = 0;
			d2W[y*NUtot+u] = 0;
                }
        }

}





void copy_dW(Layer* from) {
	for (int y=0; y<NY; y++) {
		dB[y] += from->dB[y];
//		cout << from->dB[y] << " ";
        	for (int u=0; u<NUtot; u++) {
				dW[y*NUtot+u] += from->dW[y*NUtot+u];
                }
        }
	
}

void dump_dW(ostream& os) {
	for (int y=0; y<NY; y++) {
		os << dB[y] << " ";
        	for (int u=0; u<NUtot; u++) {
				os << dW[y*NUtot+u] << " ";
                }
		os << "\n";
        }
}
void dump_W(ostream& os) {
	for (int y=0; y<NY; y++) {
		os << B[y] << " ";
        	for (int u=0; u<NUtot; u++) {
				os << W[y*NUtot+u] << " ";
                }
		os << "\n";
        }
}


~Layer() {
delete[] NK;
delete[] Y;
delete[] A;
delete[] U;

delete[] delta;
delete[] backprop;

delete[] B;
delete[] dB;
delete[] d2B;
delete[] W;
delete[] dW;
delete[] d2W;
}



// Constructor
// Categorical inputs

Layer(int t_NY, int* t_NK, int t_NU) :
	NY(t_NY), NU(t_NU)
{
NK=new int[NU];
for (int i=0; i<NU; i++)
	NK[i]=t_NK[i];
alloc(NY,NU,NK);
ninput=0;
output=0;

NUplain =0;
for (int u=0; u<NU; u++) {
  NUplain += NK[u];
}
}

// Constructor
// Real-valued inputs

Layer(int t_NY, int t_NU) :
	NY(t_NY), NU(t_NU)
{
NK=new int[NU];
for (int i=0; i<NU; i++)
	NK[i]=1;
NUr=0;
alloc(NY,NU,NK);
ninput=0;
output=0;

NUplain =0;
for (int u=0; u<NU; u++) {
  NUplain += NK[u];
}
//cout << NUplain << " " << NUtot << "a " << flush;
}

// Constructor
// Mixed inputs (NU categorical attributes, NUr real-valued)

Layer(int t_NY, int* t_NK, int t_NU, int t_NUr) :
	NY(t_NY), NU(t_NU), NUr(t_NUr)
{
int i;
NK=new int[NU+NUr];
for (i=0; i<NU; i++)
	NK[i]=t_NK[i];
for (i=NU; i<NU+NUr; i++)
	NK[i]=1;
alloc(NY,NU+NUr,NK);
ninput=0;
output=0;


NUplain =0;
for (int u=0; u<NU; u++) {
  NUplain += NK[u];
}
//cout << NY << " " << flush;
//cout << NUplain << " " << NUtot << "b " << flush;
}

Layer(istream& is);



void set_ninput(int vi) {
  ninput=vi;
};

void set_output(int vo) {
  output=vo;
};


void read(istream& is);
void write(ostream& os);

virtual void forward(int* I);
virtual void forward(double* I);
virtual void forward(int* I1, double* I2);
virtual void forward(double* I1, double* I2);
//virtual void forward(double* I, int nz_n, int* nz);

virtual double f1(int y);
virtual double f_cost(double* t);
double log_cost(double* t);
double sq_cost(double* t);

virtual double backward(double* t, double weight=1.0);

void gradient(int* I);
void gradient(double* I);
void gradient(int* I1, double* I2);
void gradient(double* I1, double* I2);
//void gradient(double* I, int nz_n,int*nz);
void gradient();


virtual void updateWeights(double epsilon);
virtual void updateWeightsL1(double epsilon);
virtual void updateWeightsClipped(double epsilon);
void resetGradient();
virtual void initWeights(int seed);

inline double* back_out() { return backprop; }
inline double* Aout() { return A; }
inline double* out() { return Y; }


inline int get_NY() { return NY; }
inline int get_NU() { return NU; }
inline int* get_NK() { return NK; }

inline double* get_dW() { return dW; }

double dlength();

void set_dW(double* newdW);


};


class Layer_tanh : public Layer
{

public:


Layer_tanh(int t_NY, int* t_NK, int t_NU) :
Layer(t_NY, t_NK, t_NU)
{
}

Layer_tanh(int t_NY, int t_NU) :
Layer(t_NY, t_NU)
{
}

Layer_tanh(int t_NY, int* t_NK, int t_NU, int t_NUr) :
Layer(t_NY, t_NK, t_NU, t_NUr)
{
}

Layer_tanh(istream& is) :
Layer(is)
{
}

Layer_tanh(Layer*from) :
Layer(from)
{
}




void forward(int* I)
{
Layer::forward(I);
squash();
}

void forward(double* I)
{
Layer::forward(I);
squash();
}

void forward(int* I1,double* I2)
{
Layer::forward(I1,I2);
squash();
}

void forward(double* I1,double* I2)
{
Layer::forward(I1,I2);
squash();
}

//void forward(double* I,int nz_n,int*nz)
//{
//Layer::forward(I,nz_n,nz);
//squash();
//}

double backward(double* t, double weight=1.0)
{
return Layer::backward(t,weight);
}



double f1(int y);

double f_cost(double* t);


void initWeights(int seed)
{
Layer::initWeights(seed);
}

void updateWeights(double epsilon)
{
Layer::updateWeights(epsilon);
}

};



class Layer_soft : public Layer
{

public :

Layer_soft(int t_NY, int* t_NK, int t_NU) :
Layer(t_NY, t_NK, t_NU)
{
}

Layer_soft(int t_NY, int t_NU) :
Layer(t_NY, t_NU)
{
}

Layer_soft(int t_NY, int* t_NK, int t_NU, int t_NUr) :
Layer(t_NY, t_NK, t_NU, t_NUr)
{
}

Layer_soft(istream& is) :
Layer(is)
{
}


Layer_soft(Layer*from) :
Layer(from)
{
}



void forward(int* I)
{
Layer::forward(I);
softmax();
}

void forward(double* I)
{
Layer::forward(I);
softmax();
}

void forward(int* I1,double* I2)
{
Layer::forward(I1,I2);
softmax();
}

void forward(double* I1,double* I2)
{
Layer::forward(I1,I2);
softmax();
}

//void forward(double* I,int nz_n,int*nz)
//{
//Layer::forward(I,nz_n,nz);
//softmax();
//}

double backward(double* t, double weight=1.0);

double f1(int y);

double f_cost(double* t);

void initWeights(int seed)
{
Layer::initWeights(seed);
}

void updateWeights(double epsilon)
{
Layer::updateWeights(epsilon);
}

};




/*

class Selector : public Layer
{
int* edges;

public:


Selector(int t_NY, int* t_NK, int t_NU) :
Layer(t_NY, t_NK, t_NU)
{
edges=new int[NY];
}

Selector(int t_NY, int t_NU) :
Layer(t_NY, t_NU)
{
edges=new int[NY];
}

Selector(int t_NY, int* t_NK, int t_NU, int t_NUr) :
Layer(t_NY, t_NK, t_NU, t_NUr)
{
edges=new int[NY];
}

Selector(istream& is) :
Layer(is)
{
edges=new int[NY];
}


void forward(int* I)
{}
void forward(double* I);
void forward(int* I1, double* I2)
{}
void forward(double* I1, double* I2)
{}

double f1(double a);

double f_cost(double* t);

double backward(double* t, double weight=1.0);

void setEdges();

void initWeights(int seed)
{
Layer::initWeights(seed);
setEdges();
}

void updateWeights(double epsilon)
{
Layer::updateWeights(epsilon);
setEdges();
}

};



*/

#endif
