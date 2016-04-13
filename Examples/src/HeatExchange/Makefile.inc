heat_exchange.pdf:
	(cd Figure; make)
	pdflatex heat_exchange.tex
LDLIBS+=$(LIBD) $(addsuffix -mt-d, -lboost_iostreams -lboost_filesystem)
