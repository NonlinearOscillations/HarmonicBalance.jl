JULIA:=julia

default: help

setup:
	${JULIA} -e 'import Pkg; Pkg.add(["JuliaFormatter", "Changelog"])'

format:
	${JULIA} -e 'using JuliaFormatter; format(".")'

changelog:
	${JULIA} -e 'using Changelog; Changelog.generate(Changelog.CommonMark(), "CHANGELOG.md"; repo = "NonlinearOscillations/HarmonicBalance.jl")'

test:
	${JULIA} --project=test -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=test test/runtests.jl

docs:
	${JULIA} --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=docs docs/make.jl

vitepress:
	npm --prefix docs i
	npm --prefix docs run docs:dev

all: setup format changelog test docs vitepress

help:
	@echo "The following make commands are available:"
	@echo " - make setup: install the dependencies for make command"
	@echo " - make format: format codes with JuliaFormatter"
	@echo " - make changelog: generate changelog"
	@echo " - make test: run the tests"
	@echo " - make docs: instantiate and build the documentation"
	@echo " - make vitepress: start Vitepress site of documentation"
	@echo " - make all: run every commands in the above order"

.PHONY: default setup format changelog test docs vitepress all help