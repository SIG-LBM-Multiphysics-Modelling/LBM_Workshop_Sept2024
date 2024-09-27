#!/bin/sh

gnuplot< plot.gnu
pdflatex SEIRD_Validation

rm *.eps
rm *.tex
rm *.log
rm *-inc-eps-converted-to.pdf
rm *.aux

pdfcrop SEIRD_Validation.pdf
rm SEIRD_Validation.pdf
mv SEIRD_Validation-crop.pdf SEIRD_Validation.pdf
