# Compiler and flags
CXX = g++
CXXFLAGS = -O2 -std=c++17 -Iinclude
LIBS = -larmadillo

# Directories
OBJDIR = build/obj
EXEDIR = build/exe

# Library sources (no main())
LIBSRC = src/analytic_solver.cpp src/jacobi_solver.cpp src/matrix_utils.cpp
LIBOBJ = $(patsubst src/%.cpp,$(OBJDIR)/%.o,$(LIBSRC))

# Programs with main()
PROGS = jacobi_test problem2 test_max_offdiag problem5 problem6

# Default target
all: $(PROGS) run_all 

# Rule for executables
$(PROGS): %: src/%.cpp $(LIBOBJ) | $(EXEDIR)
	$(CXX) $(CXXFLAGS) $^ -o $(EXEDIR)/$@ $(LIBS)

# Rule for objects
$(OBJDIR)/%.o: src/%.cpp include/%.hpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Ensure directories exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(EXEDIR):
	mkdir -p $(EXEDIR)

# Run all executables
run_all: $(PROGS)
	@echo "Running all executables..."
	@for prog in $(PROGS); do \
		echo ">>> Running $$prog"; \
		./$(EXEDIR)/$$prog; \
	done

# Clean up
clean:
	rm -rf $(OBJDIR) $(EXEDIR)
