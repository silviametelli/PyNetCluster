from setuptools import setup

setup(
    name="pyclusterseql",
    version="0.0.1",
    packages=["net_cluster"],
    #description="batch and sequential agglomerative clustering",
    #py_modules=["seql_clustering"],
    install_requires = [
        'numpy',
        'sys',
        'math',
        'operator',
    ],
    package_dir={'': 'src'},
)

