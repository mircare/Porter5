OBJECTST = Layer.o MModel.o BRNN.o NNt.o NN.o Predict_probs.o
HEADERST = Layer.h MModel.h BRNN.h NNt.h NN.h Sequence.h 

CC	= g++
CFLAGS  = -O3

Predict : $(OBJECTST)
	$(CC) $(CFLAGS) -o $@ $(OBJECTST) $(LLIBS)

Predict_probs.o : Predict_probs.cxx
	$(CC) -c $(CFLAGS) $*.cxx

Layer.o : Layer.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) $*.cxx

MModel.o : MModel.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) $*.cxx

NNt.o : NNt.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) $*.cxx

NN.o : NN.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) $*.cxx
	
BRNN.o : BRNN.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) $*.cxx


clean:
	rm *.o
