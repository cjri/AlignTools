CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -D_GLIBCXX_DEBUG 
LD_FLAGS        = -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lm -lstdc++ 
REC_OBJECTS	= aligntools.o utilities.o

align: $(REC_OBJECTS)
	$(CC) $(CC_FLAGS) $(REC_OBJECTS) -o run_align $(LD_FLAGS)
aligntools.o: aligntools.cpp
	$(CC) $(CC_FLAGS) -c aligntools.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp

