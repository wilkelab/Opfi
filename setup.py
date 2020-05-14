import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="crisposon", 
    version="0.0.3",
    author="Alexis M Hill",
    author_email="alexismhill3@gmail.com",
    description="A suite of tools for finding and annotating putative CRISPR-Transposon elements",
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
