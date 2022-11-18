#!/usr/bin/env python3.9

from setuptools import setup



setup(
    name='Not-Alike',
    version='0.0.1',
    author='Javier Montalvo',
    author_email='buitrejma@gmail.com',
    py_modules=['not_alike.nal', 'not_alike.utils'],
    packages=['not_alike'],
    python_requires='>=3.9',
    description='Pipeline that finds not alike regions of query genome compared to a hugh list of different genomes.',
    long_description = open('README.md', 'r').read(),
    install_requires =['click'],
    license = 'Python Software Foundation License',
    url='https://www.github.com/exseivier/not-alike',
    entry_points = {
        'console_scripts' : [
                'not-alike=not_alike.nal:main'
                ]
            },
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Python Software Foundation License',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        )
