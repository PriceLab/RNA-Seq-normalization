all:  docs install

docs:
	R -e "devtools::document()"
vig:
	R -e "devtools::build_vignettes()"

build:
	(cd ..; R CMD build --no-build-vignettes rnaSeqNormalization)

install:
	(cd ..; R CMD INSTALL --no-test-load rnaSeqNormalization)

check:
	(cd ..; R CMD check `ls -t rnaSeqNormalization_* | head -1`)

biocCheck:
	(cd ..; R CMD BiocCheck `ls -t rnaSeqNormalization_* | head -1`)


test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

