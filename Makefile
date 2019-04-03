
run:
	python -m predator networks/toy10_filled.xml

t: test
test:
	python -m pytest -vv test --doctest-module

black:
	black predator


.PHONY: t test black run
