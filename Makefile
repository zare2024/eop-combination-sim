# eop-combination-sim
# SOFA is built separately; copy its libsofa.a into lib/ (see README).

FC     = gfortran
FFLAGS = -O2 -g -Wall -fcheck=bounds -fbacktrace
LIBS   = -Llib -lsofa -llapack -lblas

# modules in dependency order
MODULES = src/constants.f90 src/sofa_interface.f90 src/data_io.f90 src/frames.f90 \
          src/vlbi_model.f90 src/scheduler.f90 src/lsq.f90

all: check_setup test_model vlbi_session

check_setup: src/check_setup.f90
	$(FC) $(FFLAGS) -o $@ $< $(LIBS)

test_model: $(MODULES) src/test_model.f90
	@mkdir -p build
	$(FC) $(FFLAGS) -Jbuild -o $@ $^ $(LIBS)

vlbi_session: $(MODULES) src/main_vlbi_session.f90
	@mkdir -p build
	$(FC) $(FFLAGS) -Jbuild -o $@ $^ $(LIBS)

.PHONY: all check test run clean

# toolchain check: LAPACK and SOFA
check: check_setup
	./check_setup

# tests of the model building blocks
test: test_model
	./test_model

# Milestone 1 session simulation; the report is also saved to results/summary.txt
run: vlbi_session
	@mkdir -p results
	./vlbi_session config/session.nml | tee results/summary.txt

clean:
	rm -rf build check_setup test_model vlbi_session
