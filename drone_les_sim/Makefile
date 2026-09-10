CXX = g++
CXXFLAGS = -std=c++20 -Wall -Isrc

# Trouve uniquement les fichiers .cpp qui existent vraiment
SRC = $(shell find src -type f -name "*.cpp")
OBJ = $(SRC:.cpp=.o)
EXEC = sim_drone

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(OBJ) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)