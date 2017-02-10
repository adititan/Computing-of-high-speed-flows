mkdir -p dat_files

#varying alpha
octave matlabcode.m 10.0 100 0.9 0.2
octave matlabcode.m 10.0 100 0.6 0.2
octave matlabcode.m 10.0 100 0.4 0.2
octave matlabcode.m 10.0 100 0.3 0.2
octave matlabcode.m 10.0 100 0.2 0.2
octave matlabcode.m 10.0 100 0.15 0.2

#varying grid size
octave matlabcode.m 10.0 100 0.5 0.2
octave matlabcode.m 10.0 200 0.5 0.2
octave matlabcode.m 10.0 400 0.5 0.2
octave matlabcode.m 10.0 600 0.5 0.2
octave matlabcode.m 10.0 800 0.5 0.2
octave matlabcode.m 10.0 1000 0.5 0.2

#varying CFL
octave matlabcode.m 10.0 100 0.6 0.9
octave matlabcode.m 10.0 100 0.6 0.7
octave matlabcode.m 10.0 100 0.6 0.5
octave matlabcode.m 10.0 100 0.6 0.3
octave matlabcode.m 10.0 100 0.6 0.2
octave matlabcode.m 10.0 100 0.6 0.1

ipython a2.py

pdflatex report.tex
pdflatex report.tex

rm *.png

evince report.pdf
