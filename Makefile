
run:
	python -m predator

t: test
test:
	python -m pytest -vv test --doctest-module

black:
	black predator


.PHONY: t test black run
