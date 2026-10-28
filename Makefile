# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2 -fPIC
LDFLAGS = -lmuparser -lfftw3 -lm
CXXFLAGS += -g -O0
# Targets
TARGET = test
SOURCES = main.cpp get_input_var.cpp variables.cpp time_integral.cpp dealiaser.cpp DNO_calc.cpp
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = variables.h get_input_var.h time_integral.h dealiaser.h DNO_calc.h

# Default target
all: $(TARGET)

# Link object files into final executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile source files into object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)

# Phony targets
.PHONY: all clean
