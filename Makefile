PYTHON ?= python
PYTEST ?= pytest
CTAGS ?= ctags

# .PHONY defines parts of the makefile that are not dependant on any specific file
# This is most often used to store functions
.PHONY = clean develop install tags test

develop:
	@echo "---------------Install varannot in develop mode-----------------"
	$(PYTHON) -m pip install -e .

install:
	@echo "---------------Install varannot-----------------"
	$(PYTHON) setup.py install

tags:
	$(CTAGS) --python-kinds=-i --exclude=*/tests/* -R .

test:
	@echo "---------------Run varannot tests-----------------"
	$(PYTHON) -m pip install pytest pytest-cov
	$(PYTEST) --cov-config=.coveragerc --cov-report term-missing --cov . . 

codecov:
	bash <(curl -s https://codecov.io/bash) -t c8bcf054-f0cb-4bf1-8866-1b862905ef89

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
	rm -rf ./*.pyc
	rm -rf ./*.tgz

