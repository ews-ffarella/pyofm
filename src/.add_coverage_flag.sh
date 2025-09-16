#!/usr/bin/env bash

sed -i 's/finiteVolume\/lnInclude/finiteVolume\/lnInclude -fprofile-arcs -ftest-coverage/g' Make/options
sed -i 's/-lfiniteVolume/-lfiniteVolume -lgcov/g' Make/options

