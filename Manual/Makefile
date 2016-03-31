#This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
########################################################################
# MACROS
########################################################################
NAME := manual
figures := $(wildcard figures/*.pdf) $(wildcard figures/*.png)

LATEX := pdflatex
BIBTEX := bibtex

TMPSUFFIX := .aux .log .bbl .blg .out .toc Notes.bib
TMPFILES = $(foreach sfx,$(TMPSUFFIX),$(addsuffix $(sfx), $(NAME)))
########################################################################
# IMPLICIT RULES
########################################################################
.SUFFIXES: .tex .pdf

%.pdf: %.tex
	$(LATEX) -draftmode -halt-on-error $<
	$(BIBTEX) $* || true
	$(LATEX) -draftmode $< >/dev/null
	$(LATEX) -draftmode $< >/dev/null
	$(LATEX) $<
	@grep -i "undefined" $*.log; true

########################################################################
# TARGETS
########################################################################
.PHONY: tidy

all: $(NAME).pdf

tidy:
	rm -f $(TMPFILES)

clean: tidy
	rm -f $(NAME).pdf diff.pdf

########################################################################
# RULES
########################################################################
$(NAME).pdf: $(NAME).tex $(figures)

diff.tex: $(NAME).tex $(NAME)_old.tex
	latexdiff ms_old.tex ms.tex >diff.tex
