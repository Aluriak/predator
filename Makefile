
run:
	python -m predator networks/toy10_filled.xml

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
