



#include "Model.h"




void
Model::alloc() {

	counted = new int[NY];
	nerrors_ = new int[NY];
	dcycles = new double[cycles];

	Conf=new int*[NY];
	for (int y=0;y<NY;y++)
		Conf[y]=new int[NY];
}



Model::Model(int the_NU, int the_NY, int the_NH, int the_context, int the_Moore, int the_NF,
	int the_NB, int the_NH2, int the_CoF, int the_CoB, int the_Cseg, int the_Cwin, int the_Step, 
	int the_shortcut, double* the_Thresholds, int the_cycles) :
NU(the_NU), NY(the_NY), NH(the_NH), context(the_context), Moore(the_Moore), NF(the_NF), NB(the_NB),
NH2(the_NH2), CoF(the_CoF), CoB(the_CoB), Cseg(the_Cseg), Cwin(the_Cwin), Step(the_Step), shortcut(the_shortcut), 
cycles(the_cycles)
{

Thresholds = new double[NY];
for (int y=0;y<NY-1;y++) 
	Thresholds[y] = the_Thresholds[y];

Net = new BRNN(NU,NY,NH,context,Moore,NF,NB,NH2,CoF,CoB,Step,shortcut,0);
Net->resetGradient();
NetF = new BRNN(NY*(2*Cseg+2),NY,(int)(0.5*NH),context,Moore,(int)(0.5*NF),(int)(0.5*NB),(int)(0.5*NH2),CoF,CoB,Step,shortcut,0);
NetF->resetGradient();

alloc();
}





Model::Model(istream& is) {
is >> NU >> NY >> NH >> context;
is >> NF >> NB >> NH2 >> CoF >> CoB >> Cseg >> Cwin >> Step >> shortcut >> Moore >> cycles;

Thresholds = new double[NY];
for (int y=0;y<NY-1;y++) 
	is >> Thresholds[y];

Net = new BRNN(is);
Net->resetGradient();
NetF = new BRNN(is);
NetF->resetGradient();



alloc();
}





void
Model::read(istream& is) {
is >> NU >> NY >> NH >> context;
is >> NF >> NB >> NH2 >> CoF >> CoB >> Cseg >> Cwin >> Step >> shortcut >> Moore >> cycles;
for (int y=0;y<NY-1;y++) 
	is >> Thresholds[y];

Net->read(is);
Net->resetGradient();
NetF->read(is);
NetF->resetGradient();


}




void
Model::write(ostream& os) {
os << NU << " " << NY << " " << NH << " " << context << "\n";
os <<NF<<" "<<NB<<" "<<NH2<<" "<<CoF<<" "<<CoB<<" "<<Cseg<<" "<<Cwin<<" "<<Step<<" "<<shortcut<<" "<<Moore<<" "<<cycles<<"\n";

for (int y=0;y<NY-1;y++) 
	os << Thresholds[y] << " ";
os << "\n";

Net->write(os);
NetF->write(os);
}



void
Model::randomize(int seed) {

Net->initWeights(seed);
NetF->initWeights(seed);
}



void
Model::extimation(Sequence *seq) {

int t,y;

int* O;
int a,c,cycle;//,m,maxm;
double sum=0;
double* If;

double* app=new double[NY*(seq->length+1)];


if (1) {

	sum=0;

	O=new int[seq->length+1];

	for (t=1; t<=seq->length; t++) {

		int close = 0;
		for (y=0;y<NY-1;y++) {
			if (seq->y[t]>Thresholds[y]) {
				close =y+1;
			}
		}
		O[t]= close;
		seq->yc[t] = close;
	}

	BRNN* tempo = new BRNN(Net, seq->length+1);
	tempo->resetError();

	tempo->extimation(seq->u,O,seq->length);
	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY; c++) {
			app[NY*t+c]=tempo->out()[NY*t+c];
		}
	}
	dcycles[0] += sum;

#pragma omp critical(getweights)
{
Net->copy_dW(tempo);
}
delete tempo;



	If=new double[NY*(2*Cseg+2)*(seq->length+1)];

for (cycle=1;cycle<cycles;cycle++) {
	sum=0;

	memset(If,0,NY*(2*Cseg+2)*sizeof(double)*(seq->length+1));

	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY; c++) {
			If[(NY*(2*Cseg+2))*t+c]=app[NY*t+c];
		}
		for (int cs=-Cseg;cs<=Cseg;cs++) {
			for (int tcs=t+cs*(2*Cwin+1)-Cwin;tcs<=t+cs*(2*Cwin+1)+Cwin;tcs++) {
				if (tcs>0 && tcs<=seq->length)
					for (c=0;c<NY;c++) {
						If[(NY*(2*Cseg+2))*t+NY+NY*(Cseg+cs)+c] += app[NY*tcs+c]/(2*Cwin+1);
					}
				else
					for (c=0;c<NY;c++) {
						If[(NY*(2*Cseg+2))*t+NY+NY*(Cseg+cs)+c] += 0;
					}
			}
		}
	}

	BRNN* tempo = new BRNN(NetF, seq->length+1);
	tempo->resetError();
	tempo->extimation(If,O,seq->length);

	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY; c++) {
			sum += (app[NY*t+c]-tempo->out()[NY*t+c])*
				(app[NY*t+c]-tempo->out()[NY*t+c]);
			app[NY*t+c]=tempo->out()[NY*t+c];
		}
	}

#pragma omp critical(getweightsF)
{
NetF->copy_dW(tempo);
}
delete tempo;


	dcycles[cycle] += sum;
}
	delete[] If;
	delete[] O;
} else {
}
//cout << "\n"<<flush;
delete[] app;
}





void
Model::maximization() {
Net->maximization();
NetF->maximization();
}
void
Model::maximizationL1() {
Net->maximizationL1();
NetF->maximizationL1();
}
void
Model::maximizationClipped() {
Net->maximizationClipped();
NetF->maximizationClipped();
}






void
Model::predict(Sequence* seq) {

int t,y;
int a,c,cycle;//,m,maxm;
double sum=0;
double* If;
int* O;
double* app=new double[NY*(seq->length+1)];



if (1) {

	sum=0;

	O=new int[seq->length+1];
	for (t=1; t<=seq->length; t++) {
		int close = 0;
		for (y=0;y<NY-1;y++) {
			if (seq->y[t]>Thresholds[y]) {
				close =y+1;
			}
		}
		O[t]= close;
		seq->yc[t] = close;
	}


	BRNN* tempo = new BRNN(Net, seq->length+1);
	tempo->resetError();
	tempo->predict(seq->u,seq->length);
	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY; c++) {
			sum += (app[NY*t+c]-tempo->out()[NY*t+c])*(app[NY*t+c]-tempo->out()[NY*t+c]);
			app[NY*t+c]=tempo->out()[NY*t+c];
		}
	}
	dcycles[0] += sum;
	delete tempo;

	If=new double[NY*(2*Cseg+2)*(seq->length+1)];

for (cycle=1;cycle<cycles;cycle++) {
	sum=0;

	memset(If,0,NY*(2*Cseg+2)*sizeof(double)*(seq->length+1));

	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY; c++) {
			If[(NY*(2*Cseg+2))*t+c]=app[NY*t+c];
		}
		for (int cs=-Cseg;cs<=Cseg;cs++) {
			for (int tcs=t+cs*(2*Cwin+1)-Cwin;tcs<=t+cs*(2*Cwin+1)+Cwin;tcs++) {
				if (tcs>0 && tcs<=seq->length)
					for (c=0;c<NY;c++) {
						If[(NY*(2*Cseg+2))*t+NY+NY*(Cseg+cs)+c] += app[NY*tcs+c]/(2*Cwin+1);
					}
				else
					for (c=0;c<NY;c++) {
						If[(NY*(2*Cseg+2))*t+NY+NY*(Cseg+cs)+c] += 0;
					}
			}
		}
//		O[t]=seq->y[t];
	}

	BRNN* tempo = new BRNN(NetF, seq->length+1);
	tempo->resetError();
	tempo->predict(If,seq->length);
	for (t=1; t<=seq->length; t++) {
		for (c=0; c<NY; c++) {
			sum += (app[NY*t+c]-tempo->out()[NY*t+c])*
				(app[NY*t+c]-tempo->out()[NY*t+c]);
			app[NY*t+c]=tempo->out()[NY*t+c];
		}
	}
	delete tempo;
	dcycles[cycle] += sum;
}
	delete[] If;
	delete[] O;
}

for (t=1; t<=seq->length; t++) {
	  double pred=0.0;
	  int argp=-1;

	  for (int c=0; c<NY; c++) {
		  if (app[NY*t+c]>pred) {
			  pred = app[NY*t+c];
			  argp = c;
		  }
	  }
	  seq->y_pred[t] = argp;
}


#pragma omp critical(errors)
{
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
}

delete[] app;
}



void
Model::predict(Sequence* seq, int cy) {
int temp=cycles;
cycles=cy;
predict(seq);
cycles=temp;
}












