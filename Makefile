PYTHON ?= python3
VENV ?= .venv
PIP := $(VENV)/bin/pip

DATA_DIR ?= ./data
QUERY_DIR ?= ./test/VGs
SMOKE_OUT ?= output_smoke
SMOKE_TMP ?= tmp_smoke
BLAST_BIN ?=

.PHONY: setup build-ext test smoke clean

setup:
	$(PYTHON) -m venv $(VENV)
	$(PIP) install --upgrade pip setuptools wheel
	$(PIP) install -r requirements.in

build-ext:
	. $(VENV)/bin/activate && CC=g++ python setup.py install --install-platlib=./src/

test:
	. $(VENV)/bin/activate && python -m unittest discover -s tests -v

smoke:
	. $(VENV)/bin/activate && \
	if [ -n "$(BLAST_BIN)" ]; then export PATH="$(BLAST_BIN):$$PATH"; fi; \
	rm -rf $(SMOKE_OUT) $(SMOKE_TMP) && mkdir -p $(SMOKE_OUT) && \
	python VirHostMatcher-Net.py -q $(QUERY_DIR) -o $(SMOKE_OUT) -i $(SMOKE_TMP) -n 1 -t 1 -d $(DATA_DIR)

clean:
	rm -rf build/ $(SMOKE_OUT) $(SMOKE_TMP)
