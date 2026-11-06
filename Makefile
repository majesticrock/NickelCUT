# ------------------------------------------------------------------
#   USER-CONFIGURABLE SETTINGS
# ------------------------------------------------------------------
L ?= 10
# lattice size (must be even)
CMAKE := cmake
CTEST := ctest
NPROC ?= $(shell nproc)
# how many parallel jobs make may use

# ------------------------------------------------------------------
#   BUILD-TREE NAMES
# ------------------------------------------------------------------
BUILD_DEFAULT := build/default/L$(L)
BUILD_DEBUG := build/debug/L$(L)
BUILD_BETATEST := build/BetaTest/L$(L)
BUILD_LW := build/lw-generator/L$(L)

# New "cluster" build trees (one per CPU family)
BUILD_CL := build/cascadelake/L$(L)
BUILD_ICE := build/icelake/L$(L)
ICE_ARCH := icelake-server

# ------------------------------------------------------------------
#   PUBLIC (phony) TARGETS
# ------------------------------------------------------------------
.PHONY: all clean debug BetaTest test commute mean_field \
		lw cascadelake icelake build/main.pdf cluster preview-cluster

# ------------------------------------------------------------------
#   ORIGINAL TARGETS – left untouched (kept for backward compatibility)
# ------------------------------------------------------------------
all:
	$(CMAKE) -S . -B $(BUILD_DEFAULT) -DNICKELCUT_L=$(L)
	$(CMAKE) --build $(BUILD_DEFAULT) --parallel $(NPROC)

debug:
	$(CMAKE) -S . -B $(BUILD_DEBUG) -DNICKELCUT_L=$(L) -DCMAKE_BUILD_TYPE=Debug
	$(CMAKE) --build $(BUILD_DEBUG) --parallel $(NPROC)

BetaTest:
	$(CMAKE) -S . -B $(BUILD_BETATEST) -DNICKELCUT_L=$(L) -DCMAKE_BUILD_TYPE=BetaTest
	$(CMAKE) --build $(BUILD_BETATEST) --parallel $(NPROC)

test:
	$(CMAKE) -S . -B $(BUILD_DEFAULT) -DNICKELCUT_L=$(L) -DBUILD_TESTING=ON
	$(CMAKE) --build $(BUILD_DEFAULT) --parallel $(NPROC)
	$(CTEST) --test-dir $(BUILD_DEFAULT) --output-on-failure

commute:
	$(CMAKE) -S . -B $(BUILD_DEFAULT) -DNICKELCUT_L=$(L)
	$(CMAKE) --build $(BUILD_DEFAULT) --target run_commute --parallel $(NPROC)

mean_field:
	$(CMAKE) -S . -B $(BUILD_DEFAULT) -DNICKELCUT_L=$(L)
	$(CMAKE) --build $(BUILD_DEFAULT) --target mean_field --parallel $(NPROC)

lw:
	$(CMAKE) -S . -B $(BUILD_LW) -DNICKELCUT_L=$(L) -DUSE_LW_GENERATOR=ON
	$(CMAKE) --build $(BUILD_LW) --parallel $(NPROC)

cascadelake:
	$(CMAKE) -S . -B $(BUILD_CL) -DNICKELCUT_L=$(L) -DTARGET_ARCH=cascadelake
	$(CMAKE) --build $(BUILD_CL) --parallel $(NPROC)

icelake:
	$(CMAKE) -S . -B $(BUILD_ICE) -DNICKELCUT_L=$(L) -DTARGET_ARCH=$(ICE_ARCH)
	$(CMAKE) --build $(BUILD_ICE) --parallel $(NPROC)

build/main.pdf: main.tex
	latexmk -lualatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex 1> build/log || cat build/log

clean:
	rm -rf build
	rm -rf auto_generated*

# ------------------------------------------------------------------
#   NEW “cluster” TARGET – builds both families, only flow & resume_flow
# ------------------------------------------------------------------
cluster: $(BUILD_CL)/flow $(BUILD_CL)/resume_flow $(BUILD_ICE)/flow $(BUILD_ICE)/resume_flow

# ------------------------------------------------------------------
#   Helper target – just shows what will be built (optional)
# ------------------------------------------------------------------
preview-cluster:
	@echo "The following four binaries will be built:"
	@echo "  $(BUILD_CL)/flow		  (cascadelake)"
	@echo "  $(BUILD_CL)/resume_flow   (cascadelake)"
	@echo "  $(BUILD_ICE)/flow		 (icelake)"
	@echo "  $(BUILD_ICE)/resume_flow  (icelake)"
	@echo "Run:  make cluster  (or  make cluster L=32  to change L)"

# ------------------------------------------------------------------
#   Binary targets: create the enclosing directory inside each recipe.
#   Avoiding separate directory targets prevents make from treating the
#   build tree root as a duplicate/recursive target.
# ------------------------------------------------------------------
$(BUILD_CL)/flow:
	mkdir -p $(dir $@)
	$(CMAKE) -S . -B $(BUILD_CL) -DNICKELCUT_L=$(L) -DTARGET_ARCH=cascadelake
	$(CMAKE) --build $(BUILD_CL) --target flow --parallel $(NPROC)

$(BUILD_CL)/resume_flow:
	mkdir -p $(dir $@)
	$(CMAKE) -S . -B $(BUILD_CL) -DNICKELCUT_L=$(L) -DTARGET_ARCH=cascadelake
	$(CMAKE) --build $(BUILD_CL) --target resume_flow --parallel $(NPROC)

$(BUILD_ICE)/flow: Makefile
	mkdir -p $(dir $@)
	$(CMAKE) -S . -B $(BUILD_ICE) -DNICKELCUT_L=$(L) -DTARGET_ARCH=$(ICE_ARCH)
	$(CMAKE) --build $(BUILD_ICE) --target flow --parallel $(NPROC)

$(BUILD_ICE)/resume_flow: Makefile
	mkdir -p $(dir $@)
	$(CMAKE) -S . -B $(BUILD_ICE) -DNICKELCUT_L=$(L) -DTARGET_ARCH=$(ICE_ARCH)
	$(CMAKE) --build $(BUILD_ICE) --target resume_flow --parallel $(NPROC)