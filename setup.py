import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="igseq",
    version="0.0.1",
    description="IgSeq commands",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ShawHahnLab/igseq",
    install_requires=[
        "biopython>=1.79"
        ],
    python_requires='>=3.9',
    packages=setuptools.find_packages(exclude=["test_*"]),
    package_data={'igseq': ['data/*']},
    include_package_data=True,
    entry_points={'console_scripts': [
        'igseq=igseq.__main__:main',
    ]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: POSIX :: Linux",
    ],
)
