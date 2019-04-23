# DO NOT MERGE THIS PART WITH MASTER

TARGET=GenDisCal
LOCALSRC=GenDisCal.c

# The following code should handle all the headers and source files
# with the exception of the one containing "main" (specified above)

IDIR=include
ODIR=obj
BDIR=bin
SDIR=src

MKDIR_P = mkdir -p
CC=gcc
CFLAGS= -lm -I$(IDIR) -fopenmp -O3

CFILES := $(wildcard $(SDIR)/*.c)
OFILES := $(patsubst $(SDIR)/%.c,$(ODIR)/%.o,$(CFILES))
IFILES := $(wildcard $(IDIR)/*.h)

$(ODIR):
	$(MKDIR_P) $(ODIR)

$(BDIR):
	$(MKDIR_P) $(BDIR)

$(ODIR)/%.o: $(SDIR)/%.c $(IFILES)
	$(info "Compiling" $<)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BDIR)/$(TARGET): $(OFILES) $(LOCALSRC)
	$(CC) -o $@ $^ $(CFLAGS)

