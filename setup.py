from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
classifiers = [
'Development Status :: 5 - Production/Stable',
'Intended Audience :: Science/Research',
'Operating System :: Unix',
'License :: OSI Approved :: MIT License',
'Programming Language :: Python :: 3',
]
setup(
name = 'iRRBS',
version = '0.0.1',
description = 'RRBS tool for deleting artificial cytosins',
long_description = (this_directory / "README.md").read_text(),
url = 'https://github.com/fothia/iRRBS',
download_url = 'https://github.com/fothia/iRRBS/archive/refs/tags/0.0.1.tar.gz'
author = 'Abel Fothi',
author_email = 'fothi.abel@gmail.com',
license = 'MIT',
classifiers = classifiers,
keywords = '',
install_requires = [
          'pysam',
          'pybedtools',
      ],
)
