# All source files
SRC = $(wildcard *.c)
OBJECT = $(SRC:.c=.o)
TARGET = $(SRC:.c=)

# Compiler
CXX = mpicc
CXX_FLAGS = -std=gnu99

#Build all
all:$(TARGET)

#Build all
%: %.c util.h
	$(CXX) $(CXX_FLAGS)  $< -o $@

.PHONY: clean cleantemp
clean: 
	rm -f -r $(OBJECT) $(TARGET)
cleantemp:
	rm -f -r *~
