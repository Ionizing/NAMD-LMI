#include make.inc

.PHONY: install

all: main

install: namd_lumi.x
namd_lumi.x: main
	cp src/$@ .

main:
	$(MAKE) -C src
cleansrc:
	$(MAKE) clean -C src

t: test
test: main
	$(MAKE) -C tests
cleantests:
	$(MAKE) clean -C tests

d: docs
docs: doc/Doxygenfile
	$(MAKE) -C doc
cleandocs:
	$(MAKE) clean -C doc

c: clean
clean: cleansrc cleantests cleandocs
	@rm -f *.x
