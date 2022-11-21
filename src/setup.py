from setuptools import setup

setup(
    name="pyclusterseql",
    version="0.0.1",
    description="batch and sequential agglomerative clustering",
    py_modules=["seql_clustering"],
    package_dir={'': 'src'},
)

