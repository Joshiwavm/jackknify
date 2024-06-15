from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="jacked", 
    version="0.1.0",
    author="Joshiwa van Marrewijk",
    author_email="joshiwa01@gmail.com",
    description="A brief description of your project",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Joshiwavm/jacked", 
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)