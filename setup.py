from setuptools import find_packages, setup

version = {}
with open("coldsnap/__version.py") as f:
    exec(f.read(), version)

setup(
    name="coldsnap",
    author="Weronika Osmolska",
    version=version["__version__"],
    description="Package for tracking cold spells.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Zun-1-1/ColdSnapTracking",
    license="MIT",
    packages=find_packages(where=".", exclude=["tests"]),
)