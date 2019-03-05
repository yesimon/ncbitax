# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

setup(
    name='ncbitax',
    version='0.1.1',
    description='Python API for working with NCBI taxonomy db',
    entry_points={
        'console_scripts': [
            'ncbitax = ncbitax.main:main',
            ]
    },
    long_description=readme,
    author='Simon Ye',
    author_email='mail@yesimon.com',
    url='https://github.com/yesimon/ncbitax',
    packages=find_packages(exclude=('tests', 'docs')),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
