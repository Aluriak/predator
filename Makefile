NET=toy_2_small2.sbml

run:
	python -m predator networks/$(NET) -v viz.png -vr viz-nor.png -vd viz-dag.png

t: test
tf: testff
test:
	python -m pytest -vv test predator --doctest-module
testff:  # --failed-first and --exitfirst argument
	python -m pytest -vv test predator --failed-first --exitfirst --doctest-module

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
