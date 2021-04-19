# the compilers
CC = gcc
FC = gfortran

# C compiler flags:
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -Wall

# the build target executable:
TARGET = main

all: $(TARGET)C $(TARGET)F

$(TARGET)C: $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET)C.out $(TARGET).c

$(TARGET)F: $(TARGET).f95
	$(FC) -o $(TARGET)F.out $(TARGET).f95

clean:
	$(RM) $(TARGET)
