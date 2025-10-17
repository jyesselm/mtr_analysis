#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()


with open("README.md", "r", encoding="UTF-8") as f:
    readme = f.read()

with open("requirements.txt", "r", encoding="UTF-8") as f:
    requirements = f.read().splitlines()

setup(
    name='mtr_analysis',
    version='0.1.0',
    description='some simple analysis scripts for processing mtr ribozyme data',
    long_description=readme,
    long_description_content_type="test/markdown",
    author='Joe Yesselman',
    author_email='jyesselm@unl.edu',
    url='https://github.com/jyesselm/mtr_analysis',
    packages=[
        'mtr_analysis',
    ],
    package_dir={'mtr_analysis': 'mtr_analysis'},
    py_modules=[
        'mtr_analysis/cli'
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='mtr_analysis',
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    entry_points = {
        'console_scripts' : [
        ]
    }
)
