# All source files
SRC = $(wildcard *.c)
HEADER = $(wildcard *.h)
OBJECT = $(SRC:.c=.o)
TARGET = image_seg

# Compiler
CXX = mpicc
CXX_FLAGS = -std=gnu99 -O

#Build all
all:$(TARGET)

$(TARGET): $(OBJECT)
	$(CXX) $(CXX_FLAGS) $(OBJECT) -o $(TARGET)

$(OBJECT): $(SRC)
	$(CXX) $(CXX_FLAGS) -c $(SRC)

$(SRC): $(HEADER)

.PHONY: clean
clean: 
	rm -f -r $(OBJECT) $(TARGET)

