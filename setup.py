import os

from setuptools import setup


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setupDir, 'genometreetk', 'VERSION')) as versionFile:
        return versionFile.readline().strip()


setup(
    name='GenomeTreeTk',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['genometreetk'],
    scripts=['bin/genometreetk'],
    package_data={'genometreetk': ['VERSION', 'cm_files/*']},
    url='http://pypi.python.org/pypi/genometreetk/',
    license='GPL3',
    description='A toolbox for working with genome trees.',
    python_requires='>=3.7',
    install_requires=[
        "numpy >= 1.17.0",
        "gtdblib >= 1.12.0, < 2",
        "dendropy >= 4.4.0",
        "tqdm"],
)
