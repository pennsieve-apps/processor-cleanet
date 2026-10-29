.PHONY: help image test test-local fixture run clean

SERVICE_NAME ?= processor-cleanet
IMAGE ?= $(SERVICE_NAME):dev
PYTHON ?= python3
DEMO ?= ../workflow-immune-cytometry/data/zenodo-15388817
EVENTS ?= 20000
export DOCKER_DEFAULT_PLATFORM = linux/amd64

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test        - unit tests (tests/testthat, synthetic fixture built in R) inside the image"
	@echo "make test-local  - the same tests with a local R that has Cleanet, flowCore, jsonlite, testthat"
	@echo "make fixture     - fcs-prep stand-in: $(EVENTS) events/file of the demo data -> data/fixture/{fcs,metadata}"
	@echo "                   (EVENTS=0 for every event; DEMO=<dir with the Zenodo 15388817 files>)"
	@echo "make run         - build the image and run it over data/fixture/, writing to data/output/"
	@echo "make clean       - remove data/output/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (needs requirements-dev.txt: flowio, numpy)"

image:
	docker build --platform linux/amd64 -t $(IMAGE) .

test: image
	docker run --rm --platform linux/amd64 --memory 3g --entrypoint Rscript $(IMAGE) /app/tests/testthat.R

test-local:
	Rscript tests/testthat.R

fixture:
	$(PYTHON) tests/make_fixture.py $(DEMO) data/fixture --events $(EVENTS)

run: clean
	@test -d data/fixture/fcs || { echo "no data/fixture: run make fixture first"; exit 1; }
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort | head -40

clean:
	rm -rf data/output
