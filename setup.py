from setuptools import setup, find_packages

setup(
    name='esmacs',
    version='0.1',
    packages=find_packages(),
    install_requires=['numpy'],
    scripts=['CLI/run_esmacs.py', 'CLI/run_titan.py']
)
