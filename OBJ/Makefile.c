# Makefile for GECKO-A C version
# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -O2 -fPIC -I../LIB
LDFLAGS = -lm

# Directories
SRCDIR = ../LIB
OBJDIR = .
BINDIR = ../RUN

# Target executable
TARGET = cm

# Source files (starting with converted modules)
SRCS = $(SRCDIR)/keyparameter.c \
       $(SRCDIR)/keyflag.c \
       $(SRCDIR)/minidict.c \
       $(SRCDIR)/tempoci.c \
       $(SRCDIR)/references.c \
       $(SRCDIR)/tempflag.c \
       $(SRCDIR)/main.c

# Object files
OBJS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))

# Header files
HEADERS = $(SRCDIR)/keyparameter.h \
          $(SRCDIR)/keyflag.h \
          $(SRCDIR)/minidict.h \
          $(SRCDIR)/tempoci.h \
          $(SRCDIR)/references.h \
          $(SRCDIR)/tempflag.h

# Default target
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS)

# Compile C source files to object files
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJDIR)/*.o $(TARGET)

# Phony targets
.PHONY: all clean

# Test compilation of modules (no linking)
test-compile: $(OBJS)
	@echo "All modules compiled successfully"
