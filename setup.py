r"""
`gmxbenchmark` helps evaluating performance of the GROMACS molecular 
simulation engine (www.gromacs.org).

https://github.com/ptmerz/gmxbenchmark
"""
from setuptools import setup
from os import path

#####################################
VERSION = "0.1a"
__version__ = VERSION

#####################################

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Get the requirements
with open(path.join(here, 'requirements.txt'), encoding='utf-8') as f:
    requirements = [line.strip() for line in f]

setup(
    name='gmxbenchmark',
    version=__version__,
    description='GROMACS benchmarking',
    long_description='\n' + long_description,
    url='https://gmxbenchmark.readthedocs.io',
    author='Pascal T. Merz',
    author_email='pascal.merz@colorado.edu',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development'
    ],
    keywords='benchmarking, molecular-simulation, molecular-dynamics, molecular-mechanics',
    packages=['gmxbenchmark'],
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'gmxbenchmark-prepare=gmxbenchmark.prepare:main',
            'gmxbenchmark-analyze=gmxbenchmark.analyze:main'
        ],
    },
    include_package_data=True,
    project_urls={
        'Bug Reports': 'https://github.com/ptmerz/gmxbenchmark/issues',
        'Documentation': 'https://github.com/ptmerz/gmxbenchmark',
        'Source': 'https://github.com/ptmerz/gmxbenchmark'
    },
    license="LGPLv2.1",
)
