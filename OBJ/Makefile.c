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
       $(SRCDIR)/sortstring.c \
       $(SRCDIR)/toolbox.c \
       $(SRCDIR)/dictstackdb.c \
       $(SRCDIR)/database.c \
       $(SRCDIR)/searching.c \
       $(SRCDIR)/primetool.c \
       $(SRCDIR)/rjtool.c \
       $(SRCDIR)/main.c

# Object files
OBJS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))

# Header files
HEADERS = $(SRCDIR)/keyparameter.h \
          $(SRCDIR)/keyflag.h \
          $(SRCDIR)/minidict.h \
          $(SRCDIR)/tempoci.h \
          $(SRCDIR)/references.h \
          $(SRCDIR)/tempflag.h \
          $(SRCDIR)/sortstring.h \
          $(SRCDIR)/toolbox.h \
          $(SRCDIR)/dictstackdb.h \
          $(SRCDIR)/database.h \
          $(SRCDIR)/searching.h \
          $(SRCDIR)/primetool.h \
          $(SRCDIR)/rjtool.h

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
.PHONY: all clean python clean-python

# Build Python extension
python:
	@echo "Building Python extension module..."
	cd ../python && python3 setup.py build_ext --inplace
	@echo "Python module built successfully"
	@echo "Test with: cd ../python && python3 -c 'import geckoa; print(geckoa.get_config())'"

# Clean Python build artifacts
clean-python:
	cd ../python && rm -rf build *.so *.pyc __pycache__
