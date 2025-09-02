CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -D_GLIBCXX_DEBUG -I  /opt/homebrew/Cellar/gsl/2.8/include/ -I/opt/homebrew/opt/eigen/include/eigen3/
LD_FLAGS        = -L/opt/homebrew/Cellar/gsl/2.8/lib  -lgsl -lgslcblas -lm -lstdc++ 
REC_OBJECTS	= aligntools.o utilities.o io.o distance_matrix.o diversity.o editing.o timesplit.o rgen.o

align: $(REC_OBJECTS)
	$(CC) $(CC_FLAGS) $(REC_OBJECTS) -o run_align $(LD_FLAGS)
aligntools.o: aligntools.cpp
	$(CC) $(CC_FLAGS) -c aligntools.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
distance_matrix.o: distance_matrix.cpp
	$(CC) $(CC_FLAGS) -c distance_matrix.cpp
diversity.o: diversity.cpp
	$(CC) $(CC_FLAGS) -c diversity.cpp
editing.o: editing.cpp
	$(CC) $(CC_FLAGS) -c editing.cpp
timesplit.o: timesplit.cpp
	$(CC) $(CC_FLAGS) -c timesplit.cpp
rgen.o: rgen.cpp
	$(CC) $(CC_FLAGS) -c rgen.cpp

