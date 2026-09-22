Manuscript LaTeX source.

  main.tex          clean publication-style article
  main_lineno.tex   identical content with line numbers, for review
  figures/          the five main display figures (vector PDF)
  sn-jnl.cls        Springer Nature journal class
  sn-nature.bst     Springer Nature bibliography style

Compile with pdflatex; run it three times so cross-references settle:

  pdflatex main && pdflatex main && pdflatex main

No bibtex pass is required: the bibliography is embedded as a
thebibliography environment inside the .tex file.
