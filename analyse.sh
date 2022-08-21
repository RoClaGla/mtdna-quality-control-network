#!/bin/bash

head -n1 output-4-1.050-100-0.1.csv > output-1.050-100.csv

cat output-4-1.050-100-0.1.csv output-16-1.050-100-0.1.csv output-64-1.050-100-0.1.csv | grep -v "h" >> output-1.050-100.csv
cat output-4-1.050-100-0.5.csv output-16-1.050-100-0.5.csv output-64-1.050-100-0.5.csv | grep -v "h" >> output-1.050-100.csv


Rscript plotcell.R

Rscript plot.R output-1.050-100.csv 100
