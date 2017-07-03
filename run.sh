#!/bin/bash

# clean cache
rm cache/*

matlab -nodisplay -nosplash -nodesktop -r "getImpliedDistribution; exit" | tail -n +11
soffice loadData.csv > /dev/null 2>&1 &&
firefox docs.google.com/spreadsheets/d/18GAQ_Rk8IFyiwPd2wjh_4VAdpe0bmxBvUxWXxIIhoZ4

# create graphs
python graphs/generate_graphs.py

# merge and cleanup eps files
gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=validation/validation.pdf validation/*.eps
rm validation/*.eps
