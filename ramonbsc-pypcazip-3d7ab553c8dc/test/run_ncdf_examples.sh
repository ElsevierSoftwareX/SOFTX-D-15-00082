#!/bin/bash
# This script will put the pyPCAZIP suite through its paces, exemplifying
# a lot of the most common ways it can be used.
echo " "
echo "Examples of the use of the pyPcazip toolkit with optional AMBER ncdf"
echo "format trajectory file support."
echo " "
mkdir -p output
rm -f output/*

echo " "
echo "pyPcazip: usage"
echo " "
pyPcazip -h
echo " "
echo "Example 1: Amber netcdf format trajectory, select Calpha atoms"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.ncdf -o output/2ozq_from_ncdf.pcz --selection "name CA" -vvv -p output/justCA.pdb

echo " "
echo "Example 2: Album including netcdf trajectory, select Calpha atoms"
echo " "
pyPcazip --topology 2ozq.pdb -a 2ozq+ncdf.alb -o output/2ozq_from_ncdf_alb.pcz --selection "name CA" -vvv -p output/justCA.pdb

echo " "
echo "Example 3: Unzipping to dcd format"
echo " "
pyPcaunzip --topology output/justCA.pdb --compressed output/2ozq_from_ncdf.pcz -o output/2ozq_from_ncdf.dcd -vvv
echo " "
echo " "
echo "Example 4: Print basic information from a compressed file:"
echo " "
pyPczdump --input output/2ozq_from_ncdf.pcz --info -vvv
echo " "
echo "Example 5: Print the average structure out of a compressed trajectory file:"
echo " "
pyPczdump --input output/2ozq_from_ncdf.pcz --avg --pdb output/justCA.pdb -vvv
echo " "
echo "Example 6: Print out the eigenvaluess in a compressed file."
echo " "
pyPczdump --input output/2ozq_from_ncdf.pcz --evals -vvv
echo " "
echo "Example 7: Print out a specific eigenvector."
echo " "
pyPczdump --input output/2ozq_from_ncdf.pcz --evec 1 -vvv
echo " "
echo "Example 8: Print out the projections of a specific eigenvector."
echo " "
pyPczdump --input output/2ozq_from_ncdf.pcz --proj 1 -vvv
echo " "
echo "Example 9: Basic example of pyPczcomp - usecase 2:"
echo " "
pyPcazip -i 2ozq.dcd -t 2ozq.pdb -s "name CA" -o output/2ozq.pcz
pyPczcomp -i output/2ozq.pcz output/2ozq_from_ncdf.pcz --nvecs 5
echo " "
echo " "
