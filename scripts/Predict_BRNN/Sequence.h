#ifndef Sequence_h
#define Sequence_h 1

#include <iostream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>

#define MAX_T 8196

using namespace std;

class Sequence {
public:


char name[256];

double* u;
double* y;
int* y_pred;
double* y_pred_probs;
int* yc;

int length;

int attributes;
int classes;




Sequence(istream& is, int the_attributes, int the_classes, int quot=0) {

  int i;
  char c;


attributes = the_attributes;
classes = the_classes;

  char temp[MAX_T];

if (quot == 0)  is >> name;
  is >> length;

//  cout << name << " " << length << " - " << flush;
  u = new double[attributes*(length+1)];
  y = new double[length+1];
  yc = new int[length+1];
  y_pred = new int[length+1];
  y_pred_probs = new double[classes*(length+1)];
  memset(yc,0,(length+1)*sizeof(int));
  memset(y_pred,0,(length+1)*sizeof(int));
  memset(u,0,(length+1)*attributes*sizeof(double));

  for (i=0;i<length*attributes;i++) {
    is >> u[attributes+i];
  }

if (quot == 0)
  for (i=1;i<=length;i++) {
    is >> y[i];
  }
}




void write(ostream& os) {
  int i;

  os<<name<<"\n";
  os<<length<<"\n";


  for (i=0;i<length*attributes;i++) {
    os << u[i+attributes] << " ";
  }
  os << "\n";

  for (i=1;i<=length;i++) {
    os << y[i] << " ";
  }
  os << "\n";
  for (i=1;i<=length;i++) {
    os << yc[i];
  }
  os << "\n";

  for (i=1;i<=length;i++) {
    os << y_pred[i];
  }
  os << "\n\n";
}



void write_probs(ostream& os) {
  int i,t;

  os<<name<<"\n";
  os<<length<<"\n";

  /*for (i=0;i<length*attributes;i++) {

    os << u[attributes+i] << " ";
  }
  os << "\n";

  for (i=1;i<=length;i++) {
    os << y_pred[i] << " ";
  }
  os << "\n";*/

  for (t=1; t<=length; t++) {
	  for (i=0;i<classes;i++) {
		  if (y_pred[t]==-1) os << "0.0000\t";
		  else {
			  char num[16];
			  sprintf(num, "%.4f", y_pred_probs[classes*t+i]);
			  os << num<<" ";
		  }
	  }
//	  os<<"\n";
  }

  os << "\n";


  for (i=1;i<=length;i++) {
    os << y_pred[i] << " ";
  }
  os << "\n\n";
};


void write_predictions(ostream& os) {
  int i,t;

  for (i=0;i<length*attributes;i++) {
    os << u[attributes+i] << " ";
  }
  os << "\n";

  for (i=1;i<=length;i++) {
    os << y_pred[i] << " ";
  }
  os << "\n";

  for (i=0;i<classes;i++) {
	  for (t=1; t<=length; t++) {
		  if (y_pred[t]==-1) os << "0.0000\t";
		  else {
			  char num[16];
			  sprintf(num, "%.4f", y_pred_probs[classes*t+i]);
			  os << num<<"\t";
		  }
	  }
	  os<<"\n";
  }

      os << "\n\n";
};





};
















class DataSet {
public:
  int length;
  Sequence** seq;
  int totSize;

  int attributes;
  int classes;

  DataSet() {};

  DataSet(int the_length) {
	totSize=0;
	length=the_length;
	seq = new Sequence*[length];
  }

  DataSet(istream& is, int quot=0)
    {
	  totSize=0;
      is >> length;

	  is >> attributes >> classes;

	cout << length << " sequences, " << attributes << " attributes, " << classes << " classes\n" << flush;
      seq = new Sequence*[length];
      for (int p=0; p<length; p++) {
		seq[p] = new Sequence(is,attributes,classes,quot);
		totSize += seq[p]->length;
      }
    };

  void write(ostream& os)
    {
      os << length << "\n";
      for (int p=0; p<length; p++) {
	seq[p]->write(os);
      }
    };
  void write(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
	ostream os(&outbuf);
	this->write(os);
      } else {
//	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };


  void write_probs(ostream& os)
    {
      os << length << "\n";
      for (int p=0; p<length; p++) {
		seq[p]->write_probs(os);
      }
    };
  void write_probs(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
		ostream os(&outbuf);
		this->write_probs(os);
      } else {
//	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };


  void write_predictions(ostream& os)
    {
      os << length << "\n";
      for (int p=0; p<length; p++) {
		seq[p]->write_predictions(os);
      }
    };
  void write_predictions(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
		ostream os(&outbuf);
		this->write_predictions(os);
      } else {
//	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };

};



#endif // Sequence_h
