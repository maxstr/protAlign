#!/usr/bin/env bash
# Select current version of virtualenv:
VERSION=13.0.3
# Name your first "bootstrap" environment:
INITIAL_ENV=.
# Set to whatever python interpreter you want for your first environment:
PYTHON=$(which python)
URL_BASE=https://pypi.python.org/packages/source/v/virtualenv

# --- Real work starts here ---
curl -O $URL_BASE/virtualenv-$VERSION.tar.gz
tar xzf virtualenv-$VERSION.tar.gz
# Create the first "bootstrap" environment.
$PYTHON virtualenv-$VERSION/virtualenv.py $INITIAL_ENV
# Don't need this anymore.
rm -rf virtualenv-$VERSION
# Install virtualenv into the environment.
$INITIAL_ENV/bin/pip install virtualenv-$VERSION.tar.gz
source $INITIAL_ENV/bin/activate
pip install -r requirements.txt
wget http://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.f
gfortran -static -O3 -ffast-math -lm -o TMalign TMalign.f
