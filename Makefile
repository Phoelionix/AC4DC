
CPP := /usr/local/opt/llvm/bin/clang++

SRCDIR := src
BUILDDIR := build
TARGET := bin/plasMD



SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
RAW_OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(RAW_OBJECTS:$(BUILDDIR)/Wigner/%=$(BUILDDIR)/%)
CPPFLAGS := -O3 -std=c++11 -fopenmp -stdlib=libc++
LIB := -L/usr/local/opt/llvm/lib -L$(HOME)/Programming/lib
INC := -I/usr/local/opt/llvm/include -I$(HOME)/Programming/include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CPP) $^ -v -o $(TARGET) $(LIB) "; $(CC) $^ -v -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) bin
	@echo " $(CPP) $(CPPFLAGS) $(INC) -c -o $@ $<"; $(CPP) $(CPPFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR)/%.o: $(SRCDIR)/Wigner/%.$(SRCEXT)
	@echo " $(CPP) $(CPPFLAGS) $(INC) -c -o $@ $<"; $(CPP) $(CPPFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester

.PHONY: clean
