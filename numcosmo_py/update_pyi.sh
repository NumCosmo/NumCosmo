#!/bin/sh

python ./generate_stubs.py NumCosmoMath 1.0 > ncm.pyi
python ./generate_stubs.py NumCosmo 1.0 > nc.pyi

black .
