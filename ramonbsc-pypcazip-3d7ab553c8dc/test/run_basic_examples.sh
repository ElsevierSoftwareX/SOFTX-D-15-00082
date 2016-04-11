#!/bin/bash
# This script will put the pyPCAZIP suite through its paces, exemplifying
# a lot of the most common ways it can be used.
mkdir -p output
rm -f output/*
echo "Basic examples of usage for the pyPcazip tools. Only AMBER mdcrd,"
echo "GROMACS xtc and CHARMM dcd formats are tested here."
echo " "
echo "pyPcazip: usage"
echo " "
pyPcazip -h
echo " "
echo "Example 1: DCD format trajectory, select Calpha atoms"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd -o output/2ozq.pcz --selection "name CA" -vvv -p output/justCA.pdb
echo " "
echo " "
echo "Example 2: XTC format trajectory, select Calpha atoms"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.xtc -o output/2ozq_from_xtc.pcz --selection "name CA" -vvv -p output/justCA.pdb
echo " "
echo " "
echo "Example 3: AMBeR mdcrd format trajectory, select Calpha atoms"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.mdcrd -o output/2ozq_from_xtc.pcz --selection "name CA" -vvv -p output/justCA.pdb
echo " "
echo "Example 4: with slicing of the trajectory file"
echo " "
pyPcazip --topology 2ozq.pdb -i '2ozq.dcd(1:10)' -o output/2ozq_sliced.pcz --selection "name CA" -vvv
echo " "
echo "Example 5: use --fast option and include all protein atoms"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd -o output/2ozq_all.pcz --fast -vvv
echo " "
echo "Example 6: with a mask file to select atoms to include"
echo " "
pyPcazip --topology 2ozq.pdb -i '2ozq.dcd(1:10)' -o output/2ozq_masked.pcz --mask pocket.pdb -vvv
echo " "
echo "Example 7:  with an album"
pyPcazip --topology 2ozq.pdb -a 2ozq.alb -o output/2ozq_album.pcz --mask pocket.pdb -vvv
echo " "
echo " "
echo "Example 8:  PZC4 format output"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd  -o output/2ozq_PCZ4.pcz --selection "name CA" -f PCZ4 -vvv
echo " "
echo " "
echo "Example 9:  enforce PZC6 format output"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd  -o output/2ozq_PCZ6.pcz --selection "name CA" -f PCZ6 -vvv
echo " "
echo " "
echo "Example 10:  write out a new trajectory file of the selected atoms/frames"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd  -o output/2ozq_not_made.pcz --selection "name CA" --trj_output output/CA.dcd --nopca  -vvv

echo "Example 11:  attempt to use PZC7 version"
echo " "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd -o output/2ozq_PCZ7.pcz --mask pocket.pdb -f PCZ7 -vvv

echo "pyPcaunzip - usage:"
echo " "
pyPcaunzip -h
echo " "
echo "Example 12: Unzipping to dcd format"
echo " "
pyPcaunzip --topology 2ozq.pdb --compressed output/2ozq_all.pcz -o output/2ozq_uncompressed.dcd -vvv
echo " "

echo "Example 13: Unzipping to xtc format"
echo " "
pyPcaunzip --topology 2ozq.pdb --compressed output/2ozq_all.pcz -o output/2ozq_uncompressed.xtc -vvv
echo " "

echo "pyPczdump - usage:"
echo " "
pyPczdump -h
echo " "
echo "Example 14: Print basic information from a compressed file:"
echo " "
pyPczdump --input output/2ozq.pcz --info -vvv
echo " "
echo "Example 15: Print the average structure out of a compressed trajectory file:"
echo " "
pyPczdump --input output/2ozq.pcz --avg --pdb output/justCA.pdb -vvv
echo " "
echo "Example 16: Print out the eigenvectors in a compressed file."
echo " "
pyPczdump --input output/2ozq.pcz --evals -vvv
echo " "
echo "Example 17: Print out a specific eigenvector."
echo " "
pyPczdump --input output/2ozq.pcz --evec 1 -vvv
echo " "
echo "Example 18: Print out the projections of a specific eigenvector."
echo " "
pyPczdump --input output/2ozq.pcz --proj 1 -vvv
echo " "
echo "Example 19: Print out the atomic fluctuations related to a specific eigenvector."
echo " "
pyPczdump --input output/2ozq.pcz --fluc 1 -vvv
echo " "
echo "Example 20: Produce an animation along eigenvector 1:"
echo " "
pyPczdump --input output/2ozq.pcz --anim 1 --pdb output/justCA.pdb -vvv
echo " "
echo "Example 21: Print out the rmsd of each frame from frame 1 (the second, as pyPcazip counts from 0)."
echo " "
pyPczdump --input output/2ozq.pcz --rms 1 -vvv
echo " "
echo "Example 22: Print out a collectivity metric for each eigenvector."
echo " "
pyPczdump --input output/2ozq.pcz --coll -vvv
echo " "

echo "pyPczcomp - usage:"
echo " "
pyPczcomp -h
echo " "
echo " "
echo "Example 23: Basic example of pyPczcomp:"
echo " "
pyPczcomp -i output/2ozq.pcz output/2ozq_sliced.pcz --nvecs 5
echo " "
echo " "
echo "Example 24: Basic example of pyPclust:"
echo " "
pyPczclust -i output/2ozq.pcz -o output/2ozq.clust -vv
echo " "
echo "End of basic examples."
