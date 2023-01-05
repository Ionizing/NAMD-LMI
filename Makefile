.PHONY: install test tests

all: main

install: namd_lumi.x
namd_lumi.x: main
	cp src/$@ .

main: libexternal
	$(MAKE) -C src
cleansrc:
	$(MAKE) clean -C src

t: test
test: main
	$(MAKE) -C tests
cleantests:
	$(MAKE) clean -C tests

d: docs
docs: README.md
	ford $<
cleandocs:
	rm -rf doc

e: libexternal
libexternal:
	$(MAKE) -C external
cleanexternal:
	$(MAKE) clean -C external

c: clean
clean: cleansrc cleantests
	@rm -f *.x

cleanall: clean cleanexternal cleandocs
