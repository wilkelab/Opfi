import setuptools

setuptools.setup(
    name="opfi", 
    version="0.1.1",
    author="Alexis M Hill, James Rybarski",
    author_email="alexismhill3@gmail.com",
    description="A package for discovery, annotation, and analysis of gene clusters in genomics or metagenomics datasets.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/wilkelab/Opfi",
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "biopython==1.76", 
        "pytest==5.3.2",
        "hypothesis==5.1.1", 
        "matplotlib==3.2.1", 
        "PyYAML==5.4", 
        "dna-features-viewer==3.0.1", 
        "more-itertools==8.4.0", 
        "parasail==1.2", 
        ]
)
