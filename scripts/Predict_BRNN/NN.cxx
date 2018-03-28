
// NN ver. 3.02 (27/11/2003)
//
// Copyright (C) Gianluca Pollastri 2003


#include "NN.h"


NN::NN(istream& is)
{
  is >> NO >> NH >> NI >> NIr >> which >> outp >> inp;
  upper= new Layer_soft(is);
  if (outp)
    upper->set_output(1);
  upper->set_ninput(2);

  lower=new Layer_tanh(is);
  lower->set_ninput(inp);

  int i;
  NK=new int[NI];
  NItot=0;
  for (i=0; i<NI; i++) {
	NK[i]=lower->get_NK()[i];
	NItot += NK[i];
  }
  NK2=new int[NIr];
  for (i=0; i<NIr; i++)
	NK2[i]=lower->get_NK()[NI+i];
  backprop=new double[NItot+NIr];
}


void
NN::read(istream& is)
{
  is >> NO >> NH >> NI >> NIr >> which >> outp >> inp;
  upper->read(is);
  if (outp)
  	upper->set_output(1);
  upper->set_ninput(2);

  lower->read(is);
  lower->set_ninput(inp);

  int i;
  NItot =0;
  for (i=0; i<NI; i++) {
	NK[i]=lower->get_NK()[i];
	NItot += NK[i];
  }
  for (i=0; i<NIr; i++)
	NK2[i]=lower->get_NK()[NI+i];
}



void
NN::forward(int* I)
{
  lower->forward(I);
  upper->forward(lower->out(),lower->out());
}

void
NN::forward(double* I)
{
  lower->forward(I);
  upper->forward(lower->out(),lower->out());
}

void
NN::forward(int* I1,double* I2)
{
  lower->forward(I1,I2);
  upper->forward(lower->out(),lower->out());
}

void
NN::forward(double* I1,double* I2)
{
  lower->forward(I1,I2);
  upper->forward(lower->out(),lower->out());
}

double
NN::backward(double* t, double weight)
{
  double err=upper->backward(t,weight);
  double BKD[1024];
  for (int i=0;i<NH;i++)
	BKD[i]=upper->back_out()[i];
  lower->backward(BKD,weight);
  if (inp==1)
    for (int r=NItot;r<NItot+NIr;r++)
      backprop[r]=lower->back_out()[r];
  else if (inp==2)
    for (int r=0;r<NItot+NIr;r++)
      backprop[r]=lower->back_out()[r];
  return err;
}

void
NN::gradient(int* I, double* t)
{
  upper->gradient();
  lower->gradient(I);
}
void
NN::gradient(double* I, double* t)
{
  upper->gradient();
  lower->gradient(I);
}
void
NN::gradient(int* I1,double* I2, double* t)
{
  upper->gradient();
  lower->gradient(I1,I2);
}
void
NN::gradient(double* I1,double* I2, double* t)
{
  upper->gradient();
  lower->gradient(I1,I2);
}


void
NN::resetGradient()
{
  lower->resetGradient();
  upper->resetGradient();
}

void
NN::updateWeights(double epsilon)
{
  lower->updateWeights(epsilon);
  upper->updateWeights(epsilon);
}

void
NN::updateWeightsL1(double epsilon)
{
  lower->updateWeightsL1(epsilon);
  upper->updateWeightsL1(epsilon);
}

void
NN::updateWeightsClipped(double epsilon)
{
  lower->updateWeightsClipped(epsilon);
  upper->updateWeightsClipped(epsilon);
}

void
NN::initWeights(int seed)
{
  lower->initWeights(seed);
  upper->initWeights(seed);
}


void
NN::write(ostream& os)
{
  os << NO << " " << NH<< " " << NI<< " " << NIr <<" ";
  os << which << " " << outp << " " << inp << "\n";
  upper->write(os);
  lower->write(os);
}



