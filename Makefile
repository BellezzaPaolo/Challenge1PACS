#define the compiler
CXX=g++ -std=c++20

#define the executablee
EXEC=main

#define the headers, the objects, the source files to get the executable
OBJECTS=main.o Point.o Grad.o Parameters.o
SRC=main.cpp Point.cpp Grad.cpp Parameters.cpp
HEADERS=Point.hpp Grad.hpp Parameters.hpp

# Link the executable
$(EXEC) : $(OBJECTS)
	$(CXX) $^ -o $@

# Compile the source files into object files
main.o: main.cpp $(HEADERS)
	$(CXX) -c $< -o $@

utility_functions.o: utility_functions.cpp $(HEADERS)
	$(CXX) -c $< -o $@

# Clean rule to remove object files and the executable
clean:
	rm -f *.o $(EXEC)