L ?= 10
BUILD_DIR ?= build/default/L$(L)

all:
	cmake -S . -B $(BUILD_DIR) -DNICKELCUT_L=$(L)
	cmake --build $(BUILD_DIR) --parallel

debug:
	cmake -S . -B build/debug/L$(L) -DNICKELCUT_L=$(L) -DCMAKE_BUILD_TYPE=Debug
	cmake --build build/debug/L$(L) --parallel

BetaTest:
	cmake -S . -B build/BetaTest/L$(L) -DNICKELCUT_L=$(L) -DCMAKE_BUILD_TYPE=BetaTest
	cmake --build build/BetaTest/L$(L) --parallel

test:
	cmake -S . -B $(BUILD_DIR) -DNICKELCUT_L=$(L) -DBUILD_TESTING=ON
	cmake --build $(BUILD_DIR) --parallel
	ctest --test-dir $(BUILD_DIR) --output-on-failure

commute:
	cmake -S . -B $(BUILD_DIR) -DNICKELCUT_L=$(L)
	cmake --build $(BUILD_DIR) --target run_commute --parallel

mean_field:
	cmake -S . -B $(BUILD_DIR) -DNICKELCUT_L=$(L)
	cmake --build $(BUILD_DIR) --target mean_field --parallel

lw:
	cmake -S . -B build/lw-generator/L$(L) -DNICKELCUT_L=$(L) -DUSE_LW_GENERATOR=ON
	cmake --build build/lw-generator/L$(L) --parallel

cascadelake:
	cmake -S . -B build/cascadelake/L$(L) -DNICKELCUT_L=$(L) -DTARGET_ARCH=cascadelake
	cmake --build build/cascadelake/L$(L) --parallel

build/main.pdf: main.tex
	latexmk -lualatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex 1> build/log || cat build/log

clean:
	rm -rf build
	rm -rf auto_generated*

.PHONY: all clean commute debug BetaTest cascadelake lw