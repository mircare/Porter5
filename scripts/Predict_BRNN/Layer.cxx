
// Layer ver 3.03
// 12/12/2003
// Copyright (C) Gianluca Pollastri 2003




#include "Layer.h"
#include <stdlib.h>
//#include <omp.h>

#define miny 1e-4

#define momentum 0.9
//#define momentum 0.0

//#define nthr 1

void 
Layer::softmax() {
  int y;

  int overflow=0;
  double max=A[0];
  int amax=0;
  double norm=0;
  for (y=0; y<NY; y++) {
    if (A[y]>85) {
      overflow=1;
      }
    else {
      norm += (double)exp(A[y]);
      }
    if (A[y]>max) {
      max = A[y];
      amax = y;
      }
    }

  if (overflow) {
    for (y=0; y<NY; y++) {
	Y[y] = miny;
    }
    Y[amax]=1.0-miny*(NY-1);
  } else {
    for (y=0; y<NY; y++) {
      Y[y] = (double)exp(A[y])/norm;
    }
  }
  for (y=0; y<NY; y++) {
    if (Y[y]<miny) {
	Y[y] = miny;
    }
  }
}



void
Layer::squash() {
for (int y=0; y<NY; y++)
	Y[y]=(double)tanh(A[y]);
}






void
Layer::alloc(int NY, int nu, int* NK)
{
int y,u;
NUtot=0;

for (u=0; u<nu; u++)
	NUtot += NK[u];

Y=new double[NY];
A=new double[NY];
U=new double[NUtot];

delta=new double[NY];
backprop=new double[NUtot];

W=new double[NY*NUtot];
dW=new double[NY*NUtot];
d2W=new double[NY*NUtot];
memset(d2W,0,NY*NUtot*sizeof(double));
B=new double[NY];
dB=new double[NY];
d2B=new double[NY];
memset(d2B,0,NY*sizeof(double));
}




Layer::Layer(istream& is)
{
int y,u,k;

is >> NY;
is >> NU;
is >> NUr;
NK=new int[NU+NUr];

for (u=0; u<NU+NUr; u++) is >> NK[u];
NUplain =0;
for (u=0; u<NU; u++) {
  NUplain += NK[u];
}

alloc(NY,NU+NUr,NK);

for (y=0; y<NY; y++) {
  for (u=0; u<NUtot; u++) {
      is >> W[y*NUtot+u];
    }
  is >> B[y];
  }
ninput=0;
output=0;
}



void
Layer::read(istream& is)
{
int y,u,k;

is >> NY;
is >> NU;
is >> NUr;

for (u=0; u<NU+NUr; u++) is >> NK[u];
NUplain =0;
for (u=0; u<NU; u++) {
  NUplain += NK[u];
}
NUtot = NUplain+NUr;

for (y=0; y<NY; y++) {
  for (u=0; u<NUtot; u++) {
      is >> W[y*NUtot+u];
    }
  is >> B[y];
  }
ninput=0;
output=0;

}





void Layer::write(ostream& os)
{
int y,u,k;

os << NY << "\n";
os << NU << "\n";
os << NUr << "\n";

for (u=0; u<NU+NUr; u++)
	os << NK[u] << " ";
os << "\n";

for (y=0; y<NY; y++) {
  for (u=0; u<NUtot; u++) {
      os << W[y*NUtot+u] << " ";
    }
  os << B[y] << "\n";
  }
}




void
Layer::forward(int* I)
{
int y,u,nur;
memset(U,0,NUtot*sizeof(double));

  for (y=0; y<NY; y++) {
    double a=B[y];
	nur=0;
    for (u=0; u<NU; u++) {
		if (I[u]>=0) {
			a+=W[y*NUtot+u*NK[0]+I[u]];
			U[nur+I[u]]=1.0;
		}
		nur += NK[u];
	}
	Y[y]=A[y]=a;
  }
}




void
Layer::forward(double* I)
{
int y,u,k,i;
memset(U,0,NUtot*sizeof(double));

  for (y=0; y<NY; y++) {
    i=0;
    double a=B[y];
    for (u=0; u<NUplain; u++) {
        U[i]=I[i];
        a += W[y*NUtot+u]*I[i++];
      }
    Y[y]=A[y]=a;
    }
}


void
Layer::forward(int* I1, double* I2)
{
int y,u,k,i,nur;
memset(U,0,NUtot*sizeof(double));


  for (y=0; y<NY; y++) {
    double a=B[y];
	nur=0;

    for (u=0; u<NU; u++) {
		if (I1[u]>=0) {
			a += W[y*NUtot+u*NK[0]+I1[u]];
			U[nur+I1[u]]=1.0;
		}
		nur += NK[u];
	}

	i=0;
    for (u=0; u<NUr; u++) {
        U[NUplain+i]=I2[i];
        a += W[y*NUtot+NUplain+u]*I2[i++];
      }
    Y[y]=A[y]=a;
    }
}

void
Layer::forward(double* I1, double* I2)
{
int y,u,k,i1,i2;

  for (y=0; y<NY; y++) {
    i1=0;
    i2=0;
    double a=B[y];
    for (u=0; u<NUplain; u++) {
			U[i1]=I1[i1];
			a += W[y*NUtot+u]*I1[i1];
			i1++;
      }
    for (u=0; u<NUr; u++) {
        U[NUplain+i2]=I2[i2];
        a += W[y*NUtot+NUplain+u]*I2[i2];
	   i2++;
      }
    Y[y]=A[y]=a;
    }
}
/*
void
Layer::forward(double* I, int nz_n,int*nz)
{
int y,u,k,i1,i2;

  for (y=0; y<NY; y++) {
    i1=0;
    i2=0;
    double a=B[y];
    for (u=0; u<nz_n; u++) {
			U[nz[u]]=I[nz[u]];
			a += W[y][nz[u]][0]*I[nz[u]];
      }
    Y[y]=A[y]=a;
    }
}
*/


double
Layer::f1(int y)
{
return 1.0;
}







double
Layer::f_cost(double* t)
{
double sum=0.0;
//	cout << "L"<<flush;

for (int y=0; y<NY; y++)
	sum += (t[y]-Y[y])*(t[y]-Y[y]);
return sum;
}






double
Layer::log_cost(double* t)
{
double sum=0.0;

for (int y=0; y<NY; y++) {
   if ((t[y]) && (Y[y]))
	sum -= t[y]*(double)log(Y[y]);
  }
return sum;
}




double
Layer::sq_cost(double* t)
{
double sum=0.0;

for (int y=0; y<NY; y++) {
	sum += (t[y]-Y[y])*(t[y]-Y[y]);
  }
return sum;
}





double
Layer::backward(double* rbackprop, double weight)
{
int y,u,k;

double BKD[1024];
  for (y=0; y<NY; y++)
	BKD[y]=rbackprop[y];

// If isn't an output layer
// rbackprop[] is a backprop contribution
// coming from upwards.
if (!output) {
  for (y=0; y<NY; y++) {
    BKD[y] *= f1(y);
    delta[y]=weight*BKD[y];
  }
}
// If this is an output layer
// rbackprop[] is the target vector.
else {
  for (y=0; y<NY; y++) {
    delta[y]=weight*(Y[y]-BKD[y])*f1(y);
    }
  }

double sum;
int i=0;

// If this isn't an input layer
// the backprop contribution
// must be computed, either for
// the real input part or fully.
if (ninput==1) {
  for (u=0; u<NUr; u++) {
      sum=0.0;
      for (y=0; y<NY; y++) {
		  sum += W[y*NUtot+NUplain+u]*delta[y];
        }
      backprop[NUplain+u]=sum;
    }
  }
else if (ninput==2) {
  for (u=0; u<NU+NUr; u++) {
      sum=0.0;
      for (y=0; y<NY; y++) {
		  sum += W[y*NUtot+u]*delta[y];
        }
      backprop[u]=sum;
    }
  }

double err=0.0;
if (output) {
	err=f_cost(rbackprop);
	}
else {
	for (int yyy=0;yyy<NY;yyy++)
		err+= delta[yyy]*delta[yyy];
	}
return err;
}





double
Layer_soft::backward(double* rbackprop, double weight)
{
int y,u,k;

double BKD[1024];
  for (y=0; y<NY; y++)
	BKD[y]=rbackprop[y];

// If isn't an output layer
// rbackprop[] is a backprop contribution
// coming from upwards.
if (!output) {
  for (y=0; y<NY; y++) {
    BKD[y] *= f1(y);
    delta[y]=weight*BKD[y];
  }
}
// If this is an output layer
// rbackprop[] is the target vector.
else {
  for (y=0; y<NY; y++) {
    delta[y]=weight*(Y[y]-BKD[y]);
    }
  }

double sum;
int i=0;

// If this isn't an input layer
// the backprop contribution
// must be computed, either for
// the real input part or fully.
if (ninput==1) {
  for (u=0; u<NUr; u++) {
      sum=0.0;
      for (y=0; y<NY; y++) {
		  sum += W[y*NUtot+NUplain+u]*delta[y];
        }
      backprop[NUplain+u]=sum;
    }
  }
//}
else if (ninput==2) {
  for (u=0; u<NU+NUr; u++) {
      sum=0.0;
      for (y=0; y<NY; y++) {
		  sum += W[y*NUtot+u]*delta[y];
        }
      backprop[u]=sum;
    }
  }

double err=0.0;
if (output) {
	err=f_cost(rbackprop);
	}
else {
	for (int yyy=0;yyy<NY;yyy++)
		err+= delta[yyy]*delta[yyy];
	}
return err;
}



void
Layer::gradient(int* I)
{
int y,u;

for (y=0; y<NY; y++) {
  for (u=0; u<NU; u++) {
    if (I[u]>=0) {
      dW[y*NUtot+u*NK[0]+I[u]] += delta[y];
      }
    }
  dB[y] += delta[y];
  }
}





void
Layer::gradient(double* I)
{
int y,u,k;
int i;

for (y=0; y<NY; y++) {
  i=0;
  for (u=0; u<NUplain; u++) {
      	dW[y*NUtot+u] += delta[y]*I[i++];
	}
  dB[y] += delta[y];
  }
}


void
Layer::gradient(int* I1,double* I2)
{
int y,u,k;
int i;

for (y=0; y<NY; y++) {
  for (u=0; u<NU; u++) {
    if (I1[u]>=0) {
      dW[y*NUtot+u*NK[0]+I1[u]] += delta[y];
      }
    }
  i=0;
  for (u=0; u<NUr; u++) {
      	dW[y*NUtot+NUplain+u] += delta[y]*I2[i++];
	}
  dB[y] += delta[y];
  }
}


void
Layer::gradient(double* I1,double* I2)
{
int y,u,k;
int i1,i2;

for (y=0; y<NY; y++) {
  i1=0;
    for (u=0; u<NUplain; u++) {
	dW[y*NUtot+u] += delta[y]*I1[i1];
	i1++;
    }
  i2=0;
  for (u=0; u<NUr; u++) {
	dW[y*NUtot+NUplain+u] += delta[y]*I2[i2];
	i2++;
	}
  dB[y] += delta[y];
  }
}

/*
void
Layer::gradient(double* I,int nz_n, int*nz)
{
int y,u,k;
int i1,i2;

for (y=0; y<NY; y++) {
  i1=0;
  for (u=0; u<nz_n; u++) {
      		dW[y][nz[u]][0] += delta[y]*I[nz[u]];
    }
  dB[y] += delta[y];
  }
}
*/

void
Layer::gradient()
{
gradient(U,&(U[NUplain]));
}


double sign(double a) {
if (a>0) return 1.0;
if (a<0) return -1.0;
return 0.0;
}

double clipped(double a) {
double b=sign(a)*a;
if (b>1) return sign(a)*1.0;
if (b<0.1) return sign(a)*0.1;
return a;
}


void
Layer::updateWeights(double epsilon)
{
int y,u,k;
for (y=0; y<NY; y++) {
  for (u=0; u<NUtot; u++) {
      d2W[y*NUtot+u] = momentum*d2W[y*NUtot+u] + dW[y*NUtot+u];
      W[y*NUtot+u] -= epsilon*d2W[y*NUtot+u];
    }
  d2B[y] = momentum*d2B[y] + dB[y];
  B[y] -= epsilon*d2B[y];
  }

}

void
Layer::updateWeightsL1(double epsilon)
{
int y,u,k;
double sum=0;

for (y=0; y<NY; y++) {
  for (u=0; u<NUtot; u++) {
      sum += dW[y*NUtot+u]*dW[y*NUtot+u];
    }
  sum += dB[y]*dB[y];
  }
sum = sqrt(sum);
for (y=0; y<NY; y++) {
  for (u=0; u<NUtot; u++) {
      W[y*NUtot+u] -= epsilon*dW[y*NUtot+u]/sum;
    }
  B[y] -= epsilon*dB[y]/sum;
  }
}


void
Layer::updateWeightsClipped(double epsilon)
{
int y,u,k;
for (y=0; y<NY; y++) {
  for (u=0; u<NUtot; u++) {
      W[y*NUtot+u] -= epsilon*clipped(dW[y*NUtot+u]);
    }
  B[y] -= epsilon*clipped(dB[y]);
  }

}




void
Layer::resetGradient()
{
int y,u;
memset(dW,0,NY*NUtot*sizeof(double));
memset(dB,0,NY*sizeof(double));
}



void
Layer::initWeights(int seed)
{
int y,u,k;
double D=(double)(NU+NUr);

//srand48(seed);
srand(seed);
for (y=0; y<NY; y++) {
  for (u=0; u<NUtot; u++) {
      W[y*NUtot+u] = (double)(0.5-(double)rand()/(double(RAND_MAX)))/D;
    }
  B[y] = (double)(0.5-(double)rand()/(double(RAND_MAX)))/D;
  }
}



void
Layer::set_dW(double* newdW) {
for (int y=0; y<NY; y++) {
  for (int u=0; u<NUtot; u++) {
	dW[y*NUtot+u]=newdW[y*NUtot+u];
    }
  }
}

double
Layer::dlength()
{
int y,u,k;
double sum=0.0;

for (y=0; y<NY; y++) {
  for (u=0; u<NUtot; u++) {
       sum += dW[y*NUtot+u]*dW[y*NUtot+u];
    }
  sum += dB[y]*dB[y];
  }
return sqrt(sum);
}






// Layer_soft


double
Layer_soft::f1(int y)
{
return (Y[y] - Y[y]*Y[y]);
}


double
Layer_soft::f_cost(double* t)
{
//	cout << "S"<<flush;
return Layer::log_cost(t);
}



// Layer_tanh


double
Layer_tanh::f1(int y)
{
return 1.0-(Y[y]*Y[y]);
}




double
Layer_tanh::f_cost(double* t)
{
//	cout << "T"<<flush;
return Layer::sq_cost(t);
}




