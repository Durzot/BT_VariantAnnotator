PYTHON ?= python
PYTEST ?= pytest
CTAGS ?= ctags

ifeq (,$(CONDA_PREFIX))
    CONDA_ACTIVATED=True
else
    CONDA_ACTIVATED=False
endif

#environment:
#ifeq (True, $(CONDA_ACTIVATED))
#    @echo ">>> Found $(CONDA_PREFIX)  environment activated. Install the requirements in this environment."
#    conda env update --file reqs/requirements.yml
#else
#    pip install -r requirements.txt
#endif

clean:
	$(PYTHON) setup.py clean

# installation instructions are not clear to me:
#   setup.py install -> the package cannot be imported from outside the repository even though he is 
#   visible in the conda list as <develop> in the same way as when 
#   setup.py develop -> everything is ok
#install: clean-ctags
#	$(PYTHON) setup.py install

install: ctags
	$(PYTHON) setup.py develop

develop: ctags
	$(PYTHON) setup.py develop
test:
	$(PYTEST) --cov-config=.coveragerc --cov-report term-missing --cov variant_annotator variant_annotator

ctags:
	$(CTAGS) --python-kinds=-i --exclude=*/tests/* -R variant_annotator

clean-ctags:
	rm -f tags
