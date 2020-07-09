import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="opfi", 
    version="0.1.0",
    author="Alexis M Hill, James Rybarski",
    author_email="alexismhill3@gmail.com",
    description="A suite of tools for finding and analyzing genomic systems of interest.",
    # long_description=long_description,
    # long_description_content_type="text/markdown",
    url="https://github.com/alexismhill3/CRISPR-Transposons",
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    setup_requires=['wheel'],
)
