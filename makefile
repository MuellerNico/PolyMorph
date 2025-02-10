# Compiler and flags
CXX = g++
CXXFLAGS = -Iinclude -O3 -fopenmp

# Source and executable
SRC = main.cpp
EXEC = polymorph

# Default target
all: $(EXEC)

# Build rule
$(EXEC): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

# Clean rule
clean:
	rm -f $(EXEC)
