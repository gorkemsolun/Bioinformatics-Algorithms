CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11
TARGET = basic_matrix_operations
SRC = basic_matrix_operations.cpp
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
	./$(TARGET) -n 3 -s 5
	./$(TARGET) -n 4 -s -3
	./$(TARGET) -n 5 -s 10
	./$(TARGET) -n 6 -s 2

.PHONY: all clean
