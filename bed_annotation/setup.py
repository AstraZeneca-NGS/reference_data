#!/usr/bin/env python
import sys
from os.path import join, isfile, abspath, dirname

import pip
from setuptools import setup, find_packages


name = 'bed_annotation'
package_name = 'ensembl'
version = '0.1'


print('Upgrading pip and setuptools...')
try:
    pip.main(['install', '--upgrade', 'setuptools', 'pip'])
except StandardError:
    sys.stderr.write('Cannot update pip and setuptools, that might cause errors '
                     'during the following intallation\n')

try:
    import ngs_utils
except ImportError:
    print('Installing NGS_Utils...')
    pip.main(['install', 'git+git://github.com/vladsaveliev/NGS_Utils.git'])
    import ngs_utils


setup(
    name=name,
    version=version,
    author='Vlad Saveliev',
    author_email='vladislav.sav@gmail.com',
    description='Annotation of BED files',
    keywords='bioinformatics',
    license='GPLv3',
    packages=[
        'ensembl',
        'bed_annotation',
    ],
    package_data={
        package_name: [
            'hg19/ensembl.bed.gz',
            'hg19/ensembl.bed.gz.tbi',
            'hg38/ensembl.bed.gz',
            'hg38/ensembl.bed.gz.tbi',
            'mm10/ensembl.bed.gz',
            'mm10/ensembl.bed.gz.tbi',
        ],
    },
    include_package_data=True,
    zip_safe=False,
    scripts=[
        'annotate_bed.py',
    ],
    install_requires=[
        'pybedtools',
    ],
    setup_requires=[
        'numpy'
    ],
    classifiers=[
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    test_suite='nose.collector',
    tests_require=['nose'],
)
