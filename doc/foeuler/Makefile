TEX=latex
# set latexfile to the name of the main file without the .tex
latexfile=guide
# put the names of figure files here.  include the .eps
figures=
# support subdirectories
VPATH=

pdf:$(latexfile).pdf

$(latexfile).dvi:$(figures) $(latexfile).tex
	$(TEX) $(latexfile)
	$(TEX) $(latexfile)

$(latexfile).pdf:$(latexfile).dvi
	dvips -t letter -Ppdf $(latexfile).dvi
	ps2pdf $(latexfile).ps

clean:
	rm *.aux *.dvi *.ps *.log *.toc *.pdf *.out
