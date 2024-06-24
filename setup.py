from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name="jacked", 
    version="0.1.1",
    author="Joshiwa van Marrewijk",
    author_email="joshiwa01@gmail.com",
    description="Jackknifing interferometric datasets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Joshiwavm/jacked", 
    packages=find_packages(),
    install_requires = [
        'astropy',
        'casatasks',
        'casatools',
        'jupyterlab',
        'matplotlib',
        'numpy',
        'scipy',
        'tqdm'
    ]
    )
