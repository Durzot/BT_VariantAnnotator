PYTHON ?= python
PYTEST ?= pytest
CTAGS ?= ctags

init:
	pip install -r requirements.txt

test:
	$(PYTEST) --cov-config=.coveragerc --cov-report term-missing --cov variant_annotator variant_annotator

ctags:
	$(CTAGS) --python-kinds=-i --exclude=*/tests/* variant_annotator

clean:
	rm -f tags
