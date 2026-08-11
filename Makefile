BUILD_DIR ?= build/default

all:
	cmake -S . -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR) --parallel

commute:
	cmake -S . -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR) --target run_commute --parallel

clean:
	rm -rf build
	rm -rf auto_generated*

.PHONY: all clean commute