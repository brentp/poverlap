from setuptools import setup, find_packages
import os

here = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(here, 'README.md')).read()

version = '0.1'

install_requires = [
    'nose',
]

setup(name='poverlap',
    version=version,
    description="Significance testing for genomic interval overlap",
    long_description=README,
    classifiers=[
      # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      'Topic :: Scientific/Engineering',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Topic :: Utilities'
    ],
    keywords='bioinformatics',
    author='Brent Pedersen',
    author_email='bpederse@gmail.com',
    url='https://github.com/brentp/poverlap/',
    license='BSD (2-clause)',
    #packages=['.'],
    include_package_data=True,
    test_suite='nose.collector',
    zip_safe=False,
    install_requires=install_requires,
    requires=['commandr'],
    #entry_points={
    #    'console_scripts':
    #        ['poverlap=poverlap:main']
    #}
)
