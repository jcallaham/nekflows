from setuptools import setup, find_packages

setup(
    name='nekflows',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        "numpy>=1.10.0",
        "scipy",
        "matplotlib",
        "modred",
        "pymech",
        "mpi4py",
        "dmsuite"
    ]
)