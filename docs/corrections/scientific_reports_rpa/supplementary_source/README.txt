Supplementary Information LaTeX source.

  supplementary.tex         clean Supplementary Information
  supplementary_lineno.tex  identical content with line numbers, for review
  figures/                  the eight supplementary figures S1-S8 (vector PDF)
  sn-jnl.cls                Springer Nature journal class
  sn-nature.bst             Springer Nature bibliography style

Compile with pdflatex; run it twice so cross-references settle:

  pdflatex supplementary && pdflatex supplementary
