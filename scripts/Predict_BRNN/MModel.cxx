



#include "MModel.h"




void
MModel::alloc() {

	counted = new int[NY[0]];
	nerrors_ = new int[NY[0]];
	dcycles = new double[cycles];

	Conf=new int*[NY[0]];
	for (int y=0;y<NY[0];y++)
		Conf[y]=new int[NY[0]];
}





MModel::MModel(istream& is) {
char fname[256];
filebuf inbuf;

is >> nModels;
Net=new BRNN*[nModels];
NetF=new BRNN*[nModels];

  NU=new int[nModels];
  NY=new int[nModels];
  NH=new int[nModels];

  context=new int[nModels];
  Moore=new int[nModels];

  NF=new int[nModels];
  NB=new int[nModels];
  NH2=new int[nModels];

  CoF=new int[nModels];
  CoB=new int[nModels];
  Cseg=new int[nModels];
  Cwin=new int[nModels];
  shortcut=new int[nModels];
  Step=new int[nModels];

  Thresholds=new double*[nModels];

for (int n=0;n<nModels; n++) {
  is>>fname;
  inbuf.open(fname, ios::in);
  istream is2(&inbuf);

  is2 >> NU[n]>>NY[n]>>NH[n]>>context[n];
  is2>>NF[n]>>NB[n]>>NH2[n]>>CoF[n]>>CoB[n]>>Cseg[n]>>Cwin[n];
  is2>>Step[n]>>shortcut[n]>>Moore[n]>>cycles;
  Thresholds[n] = new double[NY[n]];
  for (int y=0;y<NY[n]-1;y++) {
        is2 >> Thresholds[n][y];
  }

  Net[n] = new BRNN(is2);
  Net[n]->resetGradient();
  NetF[n] = new BRNN(is2);
  NetF[n]->resetGradient();

  inbuf.close();
}

alloc();

}







void
MModel::predict(Sequence* seq) {

int t,y;
int a,c,cycle;//,m,maxm;
double sum=0;
double* If;
int* O;
double* app=new double[NY[0]*(seq->length+1)];
memset(app,0,NY[0]*(seq->length+1)*sizeof(double));

double* final=new double[NY[0]*(seq->length+1)];
memset(final,0,NY[0]*(seq->length+1)*sizeof(double));



	sum=0;

	O=new int[seq->length+1];
	for (t=1; t<=seq->length; t++) {
		int close = 0;
		for (y=0;y<NY[0]-1;y++) {
			if (seq->y[t]>Thresholds[0][y]) {
				close =y+1;
			}
		}
		O[t]= close;
		seq->yc[t] = close;
	}

for (int n=0;n<nModels;n++) {
	Net[n]->predict(seq->u,seq->length);
	sum = 0;
	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY[n]; c++) {
//			sum += (app[NY[n]*t+c]-Net[n]->out()[NY[n]*t+c])*(app[NY[n]*t+c]-Net[n]->out()[NY[n]*t+c]);
			app[NY[n]*t+c] = Net[n]->out()[NY[n]*t+c];
			if (cycles==1)
				final[NY[n]*t+c] += Net[n]->out()[NY[n]*t+c]/nModels;
		}
	}
//	dcycles[0] += sum;

	If=new double[NY[n]*(2*Cseg[n]+2)*(seq->length+1)];

for (cycle=1;cycle<cycles;cycle++) {
	sum=0;

	memset(If,0,NY[n]*(2*Cseg[n]+2)*sizeof(double)*(seq->length+1));

	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY[n]; c++) {
			If[(NY[n]*(2*Cseg[n]+2))*t+c]=app[NY[n]*t+c];
		}
		for (int cs=-Cseg[n];cs<=Cseg[n];cs++) {
			for (int tcs=t+cs*(2*Cwin[n]+1)-Cwin[n];tcs<=t+cs*(2*Cwin[n]+1)+Cwin[n];tcs++) {
				if (tcs>0 && tcs<=seq->length)
					for (c=0;c<NY[n];c++) {
						If[(NY[n]*(2*Cseg[n]+2))*t+NY[n]+NY[n]*(Cseg[n]+cs)+c] += app[NY[n]*tcs+c]/(2*Cwin[n]+1);
					}
				else
					for (c=0;c<NY[n];c++) {
						If[(NY[n]*(2*Cseg[n]+2))*t+NY[n]+NY[n]*(Cseg[n]+cs)+c] += 0;
					}
			}
		}
//		O[t]=seq->y[t];
	}
	NetF[n]->predict(If,seq->length);
	sum = 0;
	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY[n]; c++) {
//			sum += (app[NY[n]*t+c]-NetF[n]->out()[NY[n]*t+c])*
//				(app[NY[n]*t+c]-NetF[n]->out()[NY[n]*t+c]);
			final[NY[n]*t+c] += NetF[n]->out()[NY[n]*t+c]/nModels;
		}
	}
	dcycles[cycle] += sum;
}
	delete[] If;
}
	delete[] O;


for (t=1; t<=seq->length; t++) {
	  double pred=0.0;
	  int argp=-1;

	  for (int c=0; c<NY[0]; c++) {
			seq->y_pred_probs[NY[0]*t+c]=final[NY[0]*t+c];
		  if (final[NY[0]*t+c]>pred) {
			  pred = final[NY[0]*t+c];
			  argp = c;
		  }
	  }
	  seq->y_pred[t] = argp;
}


for (t=1; t<=seq->length; t++) {

	  if (seq->y_pred[t]!=seq->yc[t]) {
		    nerrors++;
		    nerrors_[seq->yc[t]]++;
	  }

	if (seq->yc[t] != -1 && seq->y_pred[t] != -1) {
	  Conf[seq->y_pred[t]][seq->yc[t]]++;
	  counted[seq->yc[t]]++;
	}
}

delete[] app;
delete[] final;
}



void
MModel::predict(Sequence* seq, int cy) {
int temp=cycles;
cycles=cy;
predict(seq);
cycles=temp;
}














void
MModel::predict_probs(Sequence* seq) {

int t,y;
int a,c,cycle;//,m,maxm;
double sum=0;
double* If;
int* O;
double* app=new double[NY[0]*(seq->length+1)];
memset(app,0,NY[0]*(seq->length+1)*sizeof(double));

double* final=new double[NY[0]*(seq->length+1)];
memset(final,0,NY[0]*(seq->length+1)*sizeof(double));



	sum=0;

	O=new int[seq->length+1];
	for (t=1; t<=seq->length; t++) {
		int close = 0;
		for (y=0;y<NY[0]-1;y++) {
			if (seq->y[t]>Thresholds[0][y]) {
				close =y+1;
			}
		}
		O[t]= close;
		seq->yc[t] = close;
	}

for (int n=0;n<nModels;n++) {
	Net[n]->predict(seq->u,seq->length);
	sum = 0;
	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY[n]; c++) {
//			sum += (app[NY[n]*t+c]-Net[n]->out()[NY[n]*t+c])*(app[NY[n]*t+c]-Net[n]->out()[NY[n]*t+c]);
			app[NY[n]*t+c] = Net[n]->out()[NY[n]*t+c];
			if (cycles==1)
				final[NY[n]*t+c] += Net[n]->out()[NY[n]*t+c]/nModels;
		}
	}
//	dcycles[0] += sum;

	If=new double[NY[n]*(2*Cseg[n]+2)*(seq->length+1)];

for (cycle=1;cycle<cycles;cycle++) {
	sum=0;

	memset(If,0,NY[n]*(2*Cseg[n]+2)*sizeof(double)*(seq->length+1));

	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY[n]; c++) {
			If[(NY[n]*(2*Cseg[n]+2))*t+c]=app[NY[n]*t+c];
		}
		for (int cs=-Cseg[n];cs<=Cseg[n];cs++) {
			for (int tcs=t+cs*(2*Cwin[n]+1)-Cwin[n];tcs<=t+cs*(2*Cwin[n]+1)+Cwin[n];tcs++) {
				if (tcs>0 && tcs<=seq->length)
					for (c=0;c<NY[n];c++) {
						If[(NY[n]*(2*Cseg[n]+2))*t+NY[n]+NY[n]*(Cseg[n]+cs)+c] += app[NY[n]*tcs+c]/(2*Cwin[n]+1);
					}
				else
					for (c=0;c<NY[n];c++) {
						If[(NY[n]*(2*Cseg[n]+2))*t+NY[n]+NY[n]*(Cseg[n]+cs)+c] += 0;
					}
			}
		}
//		O[t]=seq->y[t];
	}
	NetF[n]->predict(If,seq->length);
	sum = 0;
	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY[n]; c++) {
//			sum += (app[NY[n]*t+c]-NetF[n]->out()[NY[n]*t+c])*
//				(app[NY[n]*t+c]-NetF[n]->out()[NY[n]*t+c]);
			final[NY[n]*t+c] += NetF[n]->out()[NY[n]*t+c]/nModels;
		}
	}
	dcycles[cycle] += sum;
}
	delete[] If;
}
	delete[] O;


	for (t=1; t<=seq->length; t++) {
	  for (int c=0; c<NY[0]; c++)
			seq->y_pred_probs[NY[0]*t+c]=final[NY[0]*t+c];

	seq->y_pred[t] = seq->yc[t];
}

delete[] app;
delete[] final;
}



void
MModel::predict_probs(Sequence* seq, int cy) {
int temp=cycles;
cycles=cy;
predict_probs(seq);
cycles=temp;
}











void
MModel::predict(Sequence* seq, double threshold) {

int t,y;
int a,c,cycle;//,m,maxm;
double sum=0;
double* If;
int* O;
double* app=new double[NY[0]*(seq->length+1)];
memset(app,0,NY[0]*(seq->length+1)*sizeof(double));

double* final=new double[NY[0]*(seq->length+1)];
memset(final,0,NY[0]*(seq->length+1)*sizeof(double));



	sum=0;

	O=new int[seq->length+1];
	for (t=1; t<=seq->length; t++) {
		int close = 0;
		for (y=0;y<NY[0]-1;y++) {
			if (seq->y[t]>Thresholds[0][y]) {
				close =y+1;
			}
		}
		O[t]= close;
		seq->yc[t] = close;
	}

for (int n=0;n<nModels;n++) {
	Net[n]->predict(seq->u,seq->length);
	sum = 0;
	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY[n]; c++) {
//			sum += (app[NY[n]*t+c]-Net[n]->out()[NY[n]*t+c])*(app[NY[n]*t+c]-Net[n]->out()[NY[n]*t+c]);
			app[NY[n]*t+c] = Net[n]->out()[NY[n]*t+c];
			if (cycles==1)
				final[NY[n]*t+c] += Net[n]->out()[NY[n]*t+c]/nModels;
		}
	}
//	dcycles[0] += sum;

	If=new double[NY[n]*(2*Cseg[n]+2)*(seq->length+1)];

for (cycle=1;cycle<cycles;cycle++) {
	sum=0;

	memset(If,0,NY[n]*(2*Cseg[n]+2)*sizeof(double)*(seq->length+1));

	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY[n]; c++) {
			If[(NY[n]*(2*Cseg[n]+2))*t+c]=app[NY[n]*t+c];
		}
		for (int cs=-Cseg[n];cs<=Cseg[n];cs++) {
			for (int tcs=t+cs*(2*Cwin[n]+1)-Cwin[n];tcs<=t+cs*(2*Cwin[n]+1)+Cwin[n];tcs++) {
				if (tcs>0 && tcs<=seq->length)
					for (c=0;c<NY[n];c++) {
						If[(NY[n]*(2*Cseg[n]+2))*t+NY[n]+NY[n]*(Cseg[n]+cs)+c] += app[NY[n]*tcs+c]/(2*Cwin[n]+1);
					}
				else
					for (c=0;c<NY[n];c++) {
						If[(NY[n]*(2*Cseg[n]+2))*t+NY[n]+NY[n]*(Cseg[n]+cs)+c] += 0;
					}
			}
		}
//		O[t]=seq->y[t];
	}


	NetF[n]->predict(If,seq->length);
	sum = 0;
	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY[n]; c++) {
			sum += (app[NY[n]*t+c]-NetF[n]->out()[NY[n]*t+c])*
				(app[NY[n]*t+c]-NetF[n]->out()[NY[n]*t+c]);
			final[NY[n]*t+c] += NetF[n]->out()[NY[n]*t+c]/nModels;
		}
	}
	dcycles[cycle] += sum;
}
	delete[] If;
}
	delete[] O;


for (t=1; t<=seq->length; t++) {
	  double pred=0.0;
	  int argp=-1;

//cout << threshold << " " << final[NY[0]*t] <<  " " << final[NY[0]*t+1] << "\n" << flush;

	  if (final[NY[0]*t]> threshold) {
		argp = 0;
	  } else {
		argp = 1;
	  }
	  seq->y_pred[t] = argp;
}


for (t=1; t<=seq->length; t++) {

	  if (seq->y_pred[t]!=seq->yc[t]) {
		    nerrors++;
		    nerrors_[seq->yc[t]]++;
	  }

	if (seq->yc[t] != -1 && seq->y_pred[t] != -1) {
	  Conf[seq->y_pred[t]][seq->yc[t]]++;
	  counted[seq->yc[t]]++;
	}
}

delete[] app;
delete[] final;
}
