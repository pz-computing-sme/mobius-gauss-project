# ============================================================================
# Gaussian Möbius Sieve - Build Configuration
# ============================================================================

CXX = g++
CXXFLAGS = -O3 -fopenmp -std=c++17 -Wall -Wextra
TARGET = sieve.exe

# Default target
all: $(TARGET)

# Build executable
$(TARGET): gaussian_mobius_sieve.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) gaussian_mobius_sieve.cpp

# Clean build artifacts
clean:
	rm -f $(TARGET) *.csv

# Run with 10^8
run: $(TARGET)
	.\$(TARGET) 100000000 mobius_gauss_1e8.csv

# Run small test (10^6)
test: $(TARGET)
	.\$(TARGET) 1000000 mobius_gauss_test.csv

.PHONY: all clean run test