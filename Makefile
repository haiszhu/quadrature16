# quadrature16 MEX bindings for near/self and mid quadrature and target projection.
# make        -> build/quadSD{,_helm,_stok}_w16_mex.m, build/quadSD{,_helm,_stok}_mid_w16_mex.m,
#                build/project_t0_h_w16_mex.m, build/prepare_quadSD_targets_w16_mex.m and their gateway
#                build/quadrature16.<mexext>
# make clean  -> remove build/
# Outside this folder: mwrap, MATLAB (headers, libmx/libmex/libmat, ILP64 libmwblas/libmwlapack), gfortran, cc.
# Override MW, MATLAB_ROOT, FC or CC on the command line.
.DEFAULT_GOAL := all
SHELL := /bin/bash
ROOT := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
SRC := $(ROOT)/src
BUILD := $(ROOT)/build

FC = gfortran
CC = cc
MW ?= $(or $(firstword $(wildcard $(HOME)/git/mwrap/mwrap)),$(shell command -v mwrap 2>/dev/null),mwrap)
ifeq ($(shell uname -s),Darwin)
  MATLAB_ROOT ?= $(lastword $(sort $(wildcard /Applications/MATLAB_R*.app)))
  MATLAB_ARCH ?= $(if $(filter arm64,$(shell uname -m)),maca64,maci64)
  MEX_LDFLAGS ?= -bundle -Wl,-undefined,dynamic_lookup
  USE_VFORCE ?= 1
else
  MATLAB_ROOT ?= $(lastword $(sort $(wildcard /usr/local/MATLAB/R* /opt/MATLAB/R* $(HOME)/MATLAB/R*)))
  MATLAB_ARCH ?= glnxa64
  MEX_LDFLAGS ?= -shared
  USE_VFORCE ?= 0
endif
MATLIB := $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)
MEX_EXT := $(patsubst mac%,mexmac%,$(patsubst glnxa64,mexa64,$(MATLAB_ARCH)))

# Near/self preserves floating-point operation order. Mid enables FMA contraction and,
# when supported by the compiler, MATMUL through external ILP64 BLAS.
FFLAGS := -O3 -fPIC -fdefault-integer-8 -frecursive -ffp-contract=off -fno-unsafe-math-optimizations -fprotect-parens
BLAS_FFLAGS := $(shell $(FC) -Q --help=fortran 2>/dev/null | awk '$$1 == "-fexternal-blas64" { print $$1 }')
MID_FFLAGS := $(filter-out -ffp-contract=%,$(FFLAGS)) -ffp-contract=fast $(BLAS_FFLAGS)
MID_CPPFLAGS := -cpp
ifeq ($(USE_VFORCE),1)
  MID_CPPFLAGS += -DSURFQUADR_USE_VFORCE
  VFORCE_LIBS := -framework Accelerate
endif
MATH_CONFIG := $(BUILD)/vforce-$(USE_VFORCE).stamp
CFLAGS := -O2 -fPIC -std=c11
MEX_CPPFLAGS := -I"$(MATLAB_ROOT)/extern/include" -DMATLAB_MEX_FILE -DMATLAB_DEFAULT_RELEASE=R2018a -DMWF77_UNDERSCORE1
MWFLAGS := -c99complex -i8 -mex quadrature16

DRIVERS := quadS_on_w16 quadD_on_w16 quadSn_on_w16 quadDn_on_w16 quadS_off_w16 quadD_off_w16 quadSn_off_w16 quadDn_off_w16 \
  quadS_helm_on_w16 quadD_helm_on_w16 quadSn_helm_on_w16 quadDn_helm_on_w16 \
  quadS_helm_off_w16 quadD_helm_off_w16 quadSn_helm_off_w16 quadDn_helm_off_w16 \
  quadS_stok_on_w16 quadD_stok_on_w16 quadSn_stok_on_w16 quadDn_stok_on_w16 quadS_stok_off_w16 quadD_stok_off_w16 \
  quadSn_stok_off_w16 quadDn_stok_off_w16
OBJECTS := $(BUILD)/quadrature16_mid.o $(BUILD)/quadrature16_utils.o $(DRIVERS:%=$(BUILD)/%.o) \
  $(BUILD)/quadrature16.o $(BUILD)/quadrature16_gateway.o
MEX := $(BUILD)/quadrature16.$(MEX_EXT)
WRAPPER := $(BUILD)/quadSD_w16_mex.m
# The gateway and Fortran objects use distinct names.
GATEWAY := $(BUILD)/quadrature16_gateway.c

.PHONY: all clean
all: $(MEX) $(WRAPPER)

$(BUILD):
	mkdir -p "$@"

# Switching the math backend rebuilds the mid object and relinks the MEX.
$(MATH_CONFIG): | $(BUILD)
	rm -f "$(BUILD)/vforce-0.stamp" "$(BUILD)/vforce-1.stamp"
	touch "$@"

$(WRAPPER): $(SRC)/quadrature16.mw | $(BUILD)
	cd "$(BUILD)" && "$(MW)" $(MWFLAGS) -mb "$<"

$(GATEWAY): $(SRC)/quadrature16.mw | $(BUILD)
	"$(MW)" $(MWFLAGS) -c "$@" "$<"

$(BUILD)/quadrature16_gateway.o: $(GATEWAY)
	$(CC) $(CFLAGS) $(MEX_CPPFLAGS) -c "$<" -o "$@"

# Modules: surfquadr_mid (quadrature16_mid.f90), surfquadr (quadrature16_utils.f90); .mod files go to build/.
$(BUILD)/quadrature16_mid.o: $(SRC)/quadrature16_mid.f90 $(MATH_CONFIG) | $(BUILD)
	$(FC) $(MID_FFLAGS) $(MID_CPPFLAGS) -J"$(BUILD)" -c "$<" -o "$@"

$(BUILD)/quadrature16_utils.o: $(SRC)/quadrature16_utils.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J"$(BUILD)" -c "$<" -o "$@"

$(BUILD)/quad%_w16.o: $(SRC)/quad%_w16.f90 $(BUILD)/quadrature16_utils.o
	$(FC) $(FFLAGS) -I"$(BUILD)" -c "$<" -o "$@"

$(BUILD)/quadrature16.o: $(SRC)/quadrature16.f90 $(BUILD)/quadrature16_utils.o $(BUILD)/quadrature16_mid.o
	$(FC) $(FFLAGS) -I"$(BUILD)" -c "$<" -o "$@"

$(MEX): $(OBJECTS)
	$(FC) $(MEX_LDFLAGS) $^ -L"$(MATLIB)" -Wl,-rpath,"$(MATLIB)" -lmx -lmex -lmat -lmwblas -lmwlapack $(VFORCE_LIBS) -o "$@"

clean:
	rm -rf "$(BUILD)"
