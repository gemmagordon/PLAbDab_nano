from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='PLAbDab-nano',
    version='0.0.1',
    description='Set of functions to scrape Genbank for immune receptor proteins',
    license='BSD 3-clause license',
    maintainer='Gemma Gordon',
    long_description=long_description,
    long_description_content_type='text/markdown',
    maintainer_email='gemma.gordon@wolfson.ox.ac.uk',
    include_package_data=True,
    packages=find_packages(include=('PLAbDab_nano', 'PLAbDab_nano.*')),
    install_requires=[
        'numpy',
        'biopython',
        'numba',
        'beautifulsoup4',
        'kasearch',
        'ablang',
        'pandas',
        'pyhmmer',
        'immunebuilder'
    ],
)
