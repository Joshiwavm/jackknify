from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name="jackknify", 
    version="0.3.0",
    author="Joshiwa van Marrewijk",
    author_email="joshiwa01@gmail.com",
    description="Jackknifing interferometric datasets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Joshiwavm/jackknify", 
    packages=find_packages(),
    install_requires = [
        'astropy',
        'casatasks',
        'casatools',
        'casadata',
        'jupyterlab',
        'matplotlib',
        'numpy',
        'scipy',
        'tqdm'
    ]
    )
