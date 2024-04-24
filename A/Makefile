FLAGS=-O3
# CUDA_FLAGS = -arch=sm_35

# CXX=g++

RM=rm -f

EXEC=pi

all: $(EXEC)

$(EXEC):
	$(CXX) $(FLAGS) $(EXEC).cpp -c -o $(EXEC).o
	$(CXX) $(FLAGS) $(EXEC).o -o $(EXEC)

clean:
	$(RM) $(EXEC).o $(EXEC)
