from setuptools import setup, find_packages
classifiers = [
“Development Status :: 5 — Production/Stable”,
“Intended Audience :: Science/Research”,
“Operating System :: Unix”,
“License :: OSI Approved :: MIT License”,
“Programming Language :: Python :: 3”,
]
setup(
name=”iRRBS”,
version=”0.0.1",
description=”RRBS tool for deleting artificial cytosins”,
long_description=open(“README.txt”).read() + “\n\n” + open(“CHANGELOG.txt”).read(),
url=””,
author=”Abel Fothi”,
author_email=”fothi.abel@gmail.com”,
license=”MIT”,
classifiers=classifiers,
keywords=””,
packages=find_packages(),
)
