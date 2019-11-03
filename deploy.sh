# Script to deploy package
pip install twine
python setup.py sdist
twine upload --verbose dist/*