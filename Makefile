PYTHON ?= python
PYTEST ?= pytest
CTAGS ?= ctags
PIP ?= pip

# .PHONY defines parts of the makefile that are not dependant on any specific file
# This is most often used to store functions
.PHONY = help test tags clean

# Defines the default target that `make` will to try to make, or in the case of a phony target, execute the specified commands
# This target is executed whenever we just type `make`
.DEFAULT_GOAL = help

# The @ makes sure that the command itself isn't echoed in the terminal
help:
	@echo "---------------HELP-----------------"
	@echo "make setup to setup the project"
	@echo "make test to run the tests"
	@echo "make clean to remove build/dev files"
	@echo "make install to install"
	@echo "------------------------------------"

test:
	$(PIP) install --upgrade pytest pytest-cov
	$(PYTEST) --cov-config=.coveragerc --cov-report term-missing --cov . . 

codecov:
	bash <(curl -s https://codecov.io/bash) -t c8bcf054-f0cb-4bf1-8866-1b862905ef89

build:
	@echo "---------------Build variant_annotator-----------------"
	$(PYTHON) -m pip install --user --upgrade setuptools wheel
	$(PYTHON) setup.py sdist bdist_wheel

install:
	@echo "---------------Install variant_annotator-----------------"
	git submodule update --init
	$(PIP) install --upgrade pip setuptools wheel
	$(PIP) install -r reqs/requirements.txt
	$(PIP) install .

update:
	@echo "---------------Update vcf2maf-----------------"
	cd tools/vcf2maf
	git pull
	cd ../..
	@echo "---------------Update ensembl-vep-----------------"
	cd tools/ensembl-vep
	git pull
	cd ../..
tags:
	$(CTAGS) --python-kinds=-i --exclude=*/tests/* -R .

clean:
	rm -f tags
	rm -rf build
	rm -rf dist
	rm -rf *.egg-info
