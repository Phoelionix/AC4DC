
CPP := /usr/local/opt/llvm/bin/clang++

SRCDIR := src
BUILDDIR := build
TARGET := bin/ac4dc



SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
# OBJECTS := $(RAW_OBJECTS:$(BUILDDIR)/Wigner/%=$(BUILDDIR)/%)
CPPFLAGS := -std=c++11 -fopenmp # -O3
LIB := -L/usr/local/opt/llvm/lib -L$(HOME)/Programming/lib -fopenmp
INC := -I/usr/local/opt/llvm/include -I$(HOME)/Programming/include

$(TARGET): $(OBJECTS)
	@echo " Linking... "
	@echo " $(CPP) $^ -o $(TARGET) $(LIB) "; $(CPP) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) bin $(BUILDDIR)/Wigner
	@echo " $(CPP) $(CPPFLAGS) $(INC) -c -o $@ $<"; $(CPP) $(CPPFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

debug: $(OBJECTS)
	@echo " Linking with lldb flags... "

# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester

.PHONY: clean
