import setuptools
import re
import sys

with open("README.md", "r") as fh:
    long_description = fh.read()

REGEX_COMMENT = re.compile(r'[\s^]#(.*)')


def parse_requirements(filename):
    with open(filename, 'rt') as filehandle:
        return tuple(filter(None, (
            REGEX_COMMENT.sub('', line).strip()
            for line in filehandle
        )))


setuptools.setup(
    name="reference_parser-MattWellie",
    version="0.1.0",
    author="MatthewWelland",
    author_email="matthewjwelland@gmail.com",
    description="A reference-sequence translator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MattWellie/Final-Referencer",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=parse_requirements("requirements/requirements.txt"),
)