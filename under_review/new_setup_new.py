#!/usr/bin/env python

"""Setup script for the PolyPhys package."""

from setuptools import setup, find_packages

# Read long description from README and HISTORY files
with open("README.rst", encoding="utf-8") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst", encoding="utf-8") as history_file:
    history = history_file.read()

# Runtime and test dependencies
requirements = ["Click>=8.2"]
test_requirements = ["pytest>=8.3"]

setup(
    name="polyphys",
    version="0.3.0",
    author="Amir Sadeghi",
    author_email="amirh.sadeghi@outlook.com",
    description=(
        "A Python package for managing polymer physics and molecular dynamics "
        "projects, focusing on computational modeling of bacterial chromosome "
        "organization."
    ),
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/x-rst",
    url="https://github.com/amirhs1/PolyPhys",
    packages=find_packages(include=["polyphys", "polyphys.*"]),
    entry_points={
        "console_scripts": [
            "polyphys=polyphys.cli:main",
        ],
    },
    include_package_data=True,
    install_requires=requirements,
    python_requires=">=3.12",
    license="MIT license",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
    ],
    keywords="polyphys polymer-physics molecular-dynamics",
    test_suite="tests",
    tests_require=test_requirements,
    zip_safe=False,
)
