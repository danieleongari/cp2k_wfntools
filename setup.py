#!/usr/bin/env python
from __future__ import absolute_import
from setuptools import setup, find_packages

if __name__ == '__main__':
    setup(
        packages=find_packages(),
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',
        name="cp2k-wfntools",
        author="Daniele Ongari",
        author_email="daniele.ongari@epfl.ch",
        description="Tools for manipulating cp2k wave functions",
        url="https://github.com/danieleongari/cp2k_wfntools",
        license="Creative Commons",
        classifiers=["Programming Language :: Python"],
        version="0.1.0",
        install_requires=["numpy", "scipy"],
        scripts=["cp2k_wfntools/wfntools"],
        extras_require={
            "testing": [
                "mock==2.0.0", "pgtest==1.1.0", "sqlalchemy-diff==0.1.3",
                "wheel>=0.31", "coverage"
            ],
            "pre-commit": [
                "pre-commit==1.11.0",
                "yapf==0.24.0",
                "prospector==1.1.5",
                "pylint",
            ]
        })
