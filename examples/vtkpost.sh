#!/usr/bin/env bash

inputFile=./vtk_examples/ELEMENTS.vtu
resFile=./vtk_examples/ss_09.h5
outFile=./vtk_examples/ELEMENTS_WITH_RESULTS.vtu

python3 test_vtkpost.py --vtuFile $inputFile --resultsFile $resFile --outputFile $outFile