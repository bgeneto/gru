# project name (generate executable with this name)
TARGET   = gru

# compiler and options
CC       = gcc
CFLAGS   = -std=c99 -Wall
OMPFLAGS = -O2 -fopenmp
LDFLAGS  = -lm

# change these to set the proper directories where each files shoould be
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

OBJECTS  = $(OBJDIR)/mytools.o $(OBJDIR)/mc-ies.o $(OBJDIR)/mc.o $(OBJDIR)/mc-calc.o $(OBJDIR)/gru-es.o

all: $(OBJECTS)
	$(CC) $(CFLAGS) $(TARGET).c $(OBJECTS) -o $(TARGET) $(LDFLAGS)

omp: $(OBJECTS)
	$(CC) $(CFLAGS) $(OMPFLAGS) $(TARGET).c $(OBJECTS) -o $(TARGET) $(LDFLAGS)

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJDIR)/*.o $(TARGET)