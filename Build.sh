#!/bin/sh

# Tex this mofo -
latex $1.tex
bibtex $1
latex $1.tex
latex $1.tex

# Make the pdf -
dvipdfm -p "letter" $1.dvi
