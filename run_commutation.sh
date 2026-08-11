#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

cmake -S . -B build/default
cmake --build build/default --target run_commute --parallel

if [ -f main.tex ]; then
    mkdir -p build
    latexmk -lualatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex 1> build/log || cat build/log
fi