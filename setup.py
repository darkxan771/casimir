import pathlib

from setuptools import find_packages, setup

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / "README.md").read_text(encoding="utf-8")


setup(
    name="hypocasimir",
    version="0.1",
    description="A package to manipulate Lie groups and their representations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/darkxan771/casimir",
    author="Pierre-Loïc Méliot",
    author_email="pierreloic.meliot@gmail.com",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Mathematicians",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.10",
    install_requires=["matplotlib"],
)
