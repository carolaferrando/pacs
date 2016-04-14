#!/bin/bash
make distclean
make 
./main -p parametersL2.pot 
mv result.dat resultL2.dat
./main -p parametersH1.pot
mv result.dat resultH1.pot
