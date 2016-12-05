#!/bin/bs

# clean cache
rm cache/*

matlab -nodisplay -nosplash -nodesktop -r "getImpliedDistribution; exit" | tail -n +12
soffice loadData.csv > /dev/null 2>&1 &&
firefox docs.google.com/spreadsheets/d/18GAQ_Rk8IFyiwPd2wjh_4VAdpe0bmxBvUxWXxIIhoZ4

# merge and cleanup eps files
gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=validation/validation.pdf validation/*.eps
rm validation/*.eps
