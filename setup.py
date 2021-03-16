from setuptools import setup, find_packages
setup(
    name="bambu",
    version='0.0.1',
    packages=find_packages(),
    author="Isadora Leitzke Guidotti",
    author_email="leitzke.gi@gmail.com",
    description="bambu",
    keywords="bioinformatics machine-learning data science",
    entry_points = {'cli': ['bambu-preprocess = bambu.preprocessing:main']}
    )