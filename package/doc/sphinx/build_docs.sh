#!/usr/bin/env bash

#./build_docs.sh 2>&1 | tee html_output

eval "$(mamba shell hook --shell bash)"
mamba activate mdanalysis-dev

rm -rfv ../html/* && make html

cd ..

#python -m http.server 9191
