from setuptools import setup

setup(
    name="PyNetCluster",
    version="0.0.1",
    packages=["pynetclustering"],
    #description="batch and sequential agglomerative clustering",
    #py_modules=["seql_clustering"],
    install_requires = [
        'numpy',
        'sys',
        'os'
        'math',
        'operator',
    ],
    package_dir={'': 'pynetclustering'},
)

