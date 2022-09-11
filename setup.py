from setuptools import setup, find_packages, Extension

setup(
    name="spectralkmeans",
    author="Daniel Palmor and May Siva",
    packages=find_packages(where='.'),
    install_requires=['invoke'],
    classifiers=["Programming Language :: Python:: Implementation :: CPython"],
    ext_modules=[Extension("spectralkmeans", ["spkmeansmodule.c", "spkmeans.c",
                                              "kmeans.c", "matrixOperations.c",
                                              "vectorOperations.c"])],
    headers=["spkmeans.h", "matrixOperations.h", "vectorOperations.h"],
)
