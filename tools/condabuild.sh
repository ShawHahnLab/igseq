#!/usr/bin/env bash

# As helpfully explained in:
# https://medium.com/analytics-vidhya/publish-a-python-package-to-conda-b352eb0bcb2e
conda build -c conda-forge -c bioconda -c defaults conda --python=3.9 --output-folder conda-out
# Then just:
# anaconda upload -u ShawHahnLab conda-out/linux-64/igseq-0.0.1-py39_0.tar.bz2
