all:  roxy install

roxy:
	R -e "devtools::document()"

build:
	(cd ..; R CMD build --no-build-vignettes rnaSeqNormalizer)

install:
	(cd ..; R CMD INSTALL --no-test-load rnaSeqNormalizer)

check:
	(cd ..; R CMD check `ls -t rnaSeqNormalizer_* | head -1`)

biocCheck:
	(cd ..; R CMD BiocCheck `ls -t rnaSeqNormalizer_* | head -1`)

vig:
	R -e "devtools::build_vignettes()"

site:
	R -e "devtools::build_site()"

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

