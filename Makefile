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

clean: cleansrc cleantests
	@rm -f *.x
