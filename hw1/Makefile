CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11
TARGET = hw1
SRC = hw1.cpp
OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)
	rm -f $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

compile: $(TARGET)

run:

.PHONY: all clean
