mkdir -p dat_files

octave matlabcode.m 200 
octave matlabcode.m 400
octave matlabcode.m 800 

ipython a3.py

pdflatex report.tex
pdflatex report.tex

rm *.png
rm *.log
rm *.aux

evince report.pdf
