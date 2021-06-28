MAIN = BPG_TEST
CXX = g++
CXX_FLAGS = -std=c++17 -O3 -Wall -Wextra -g
CPP_FILES = $(wildcard src/*.cpp)
OBJ_FILES = $(addprefix obj/, $(notdir $(CPP_FILES:.cpp=.o)))
DEP_FILES = $(OBJ_FILES:.o=.d) 
INCLUDE = -I./inc 
LINK = -L /home/sdawson/projects/def-shi/sdawson/FFTW/fftw-3.3.8/lib -lfftw3

# PHONY - means rules are executed even if files exist
.PHONY: all clean

# linking executables from object files
# all: $(MAIN)

$(MAIN): $(OBJ_FILES)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -o $@ $^ $(LINK)

-include $(DEP_FILES)

clean:
	$(RM) $(OBJ_FILES) $(DEP_FILES) $(MAIN)

obj/%.o: src/%.cpp Makefile
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -MMD -MP -c -o $@ $<


