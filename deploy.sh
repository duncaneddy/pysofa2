# Script to deploy package
pip install twine
python setup.py sdist
twine upload --username=$(PYPI_USERNAME) --password=$(PYPI_PASSWORD)  --verbose dist/*