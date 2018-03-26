from setuptools import setup, find_packages

setup(
    name='esmacs',
    version='0.1',
    packages=find_packages(),
    install_requires=['numpy'],
    entry_points={
        'console_scripts': ['esmacs = esmacs.cli:run_esmacs']
    }
)
