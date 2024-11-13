#!/bin/bash

cp ../INPUT/gecko.nml ./
cp ../INPUT/cheminput.dat ./
cp ../OBJ/cm ./

./cm

rm cm gecko.nml cheminput.dat

