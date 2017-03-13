SHELL := /usr/bin/env bash

# configure engine
## LaTeX engine
### LaTeX workflow: pdf; xelatex; lualatex
latexmkEngine := xelatex
### pandoc workflow: pdflatex; xelatex; lualatex
pandocEngine := xelatex
## HTML
HTMLVersion := html5

MD := $(wildcard *.md)
# MD2TeX := $(patsubst %.md,%.tex,$(MD))
MD2PDF := $(patsubst %.md,%.pdf,$(MD))
MD2HTML := $(patsubst %.md,%.html,$(MD))

CSSURL:=https://ickc.github.io/markdown-latex-css
# command line arguments
pandocArgCommon := -f markdown+autolink_bare_uris-fancy_lists --toc --normalize -S -V linkcolorblue -V citecolor=blue -V urlcolor=blue -V toccolor=blue --latex-engine=$(pandocEngine) -M date="`date "+%B %e, %Y"`" --filter=pandoc-citeproc --bibliography=citation.bib --csl=ieee-with-url.csl
## MD
# pandocArgMD := -f markdown+abbreviations+autolink_bare_uris+markdown_attribute+mmd_header_identifiers+mmd_link_attributes+mmd_title_block+tex_math_double_backslash-latex_macros-auto_identifiers -t markdown+raw_tex-native_spans-simple_tables-multiline_tables-grid_tables-latex_macros --normalize -s --wrap=none --column=999 --atx-headers --reference-location=block --file-scope
## TeX/PDF
### LaTeX workflow
latexmkArg := -$(latexmkEngine)
pandocArgFragment := $(pandocArgCommon)
### pandoc workflow
pandocArgStandalone := $(pandocArgFragment) --toc-depth=1 -s -N
## HTML/ePub
pandocArgHTML := $(pandocArgFragment) -t $(HTMLVersion) --toc-depth=2 -s -N --mathjax -c $(CSSURL)/css/common.css -c $(CSSURL)/fonts/fonts.css
# pandocArgePub := $(pandocArgFragment) --toc-depth=2 -s -N --mathjax -c $(CSSURL)/css/common.css -c $(CSSURL)/fonts/fonts.css -t $(ePubVersion) --epub-chapter-level=2 --self-contained
# GitHub README
pandocArgReadmeGitHub := $(pandocArgCommon) --toc-depth=2 -s -t markdown_github --reference-location=block

####################################################################################################################################

all: $(MD2PDF) $(MD2HTML)
	echo $(MD)

clean:
	rm -f  $(MD2PDF) $(MD2HTML)

####################################################################################################################################

### md to PDF: $(pandocPDF)
%.pdf: %.md
	pandoc $(pandocArgStandalone) -o $@ $<

### md to HTML: $(pandocHTML)
%.html: %.md
	pandoc $(pandocArgHTML) -o $@ $<
