BUILD_DIR ?= build/default

all:
	cmake -S . -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR) --parallel

debug:
	cmake --preset debug
	cmake --build --preset debug --parallel

test:
	cmake --preset test
	cmake --build --preset test --parallel
	ctest --preset test

commute:
	cmake -S . -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR) --target run_commute --parallel

mean_field:
	cmake -S . -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR) --target mean_field --parallel

build/main.pdf: main.tex
	latexmk -lualatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex 1> build/log || cat build/log

clean:
	rm -rf build
	rm -rf auto_generated*

.PHONY: all clean commute debug