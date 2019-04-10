
run:
	python -m predator networks/toy_1.sbml -v viz.png -vr viz-nor.png

t: test
test:
	python -m pytest -vv test --doctest-module

black:
	black predator


.PHONY: t test black run


# release cycle recipes
fullrelease:
	fullrelease
install_deps:
	python -c "import configparser; c = configparser.ConfigParser(); c.read('setup.cfg'); print(c['options']['install_requires'])" | xargs pip install -U
install:
	python setup.py install
