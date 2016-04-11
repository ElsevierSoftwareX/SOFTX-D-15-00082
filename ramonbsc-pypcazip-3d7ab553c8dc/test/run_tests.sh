#!/bin/bash
# This script will put the pyPCAZIP suite through its paces, exemplifying
# a lot of the most common ways it can be used.
function test {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "FAILED" >&2
    else
        echo "PASSED" >&2
    fi
    return $status
}

mkdir -p output
rm -f output/*

echo -n "Test 1a - DCD format compression/decompression: "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd -o output/test1.pcz --selection "name CA"  -p output/justCA.pdb
pyPcaunzip -c output/test1.pcz -t output/justCA.pdb -o output/test1.dcd
test ./traj_check.py output/justCA.pdb output/test1.dcd reference/test1.dcd 0.001 
echo -n "Test 1b - DCD format compression/decompression: "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd -o output/test1.pcz --selection "name CA"  -p output/justCA.pdb --trj_output output/ref1.dcd -q 95
pyPcaunzip -c output/test1.pcz -t output/justCA.pdb -o output/test1.dcd
test ./traj_check.py output/justCA.pdb output/test1.dcd output/ref1.dcd 0.2 
echo -n "Test 2a - XTC format compression/decompression: "
pyPcazip --topology 2ozq.pdb -i 2ozq.xtc -o output/test2.pcz --selection "name CA"  -p output/justCA.pdb
pyPcaunzip -c output/test2.pcz -t output/justCA.pdb -o output/test2.xtc 
test ./traj_check.py output/justCA.pdb output/test2.xtc reference/test2.xtc 0.002 
echo -n "Test 2b - XTC format compression/decompression: "
pyPcazip --topology 2ozq.pdb -i 2ozq.xtc -o output/test2.pcz --selection "name CA"  -p output/justCA.pdb --trj_output output/ref2.xtc -q 95
pyPcaunzip -c output/test2.pcz -t output/justCA.pdb -o output/test2.xtc
test ./traj_check.py output/justCA.pdb output/test2.xtc output/ref2.xtc 0.2 
echo -n "Test 3 - MDCRD format compression/DCD decompression: "
pyPcazip --topology 2ozq.pdb -i 2ozq.mdcrd -o output/test3.pcz --selection "name CA" -p output/justCA.pdb
pyPcaunzip -c output/test3.pcz -t output/justCA.pdb -o output/test3.dcd
test ./traj_check.py output/justCA.pdb output/test3.dcd reference/test3.dcd 0.001 
echo -n "Test 4 - Album compression/XTC decompression: "
pyPcazip --topology 2ozq.pdb -a 2ozq.alb -o output/test4.pcz --selection "name CA" -p output/justCA.pdb
pyPcaunzip -c output/test4.pcz -t output/justCA.pdb -o output/test4.xtc
test ./traj_check.py output/justCA.pdb output/test4.xtc reference/test4.xtc 0.002 
echo -n "Test 5 - NCDF compression/DCD decompression: "
python -c 'import netCDF4' >& /dev/null
if [ $? -ne 0 ]; then
echo "netCDF4 not available - SKIPPED"
else
pyPcazip --topology 2ozq.pdb -i 2ozq.ncdf -o output/test5.pcz --selection "name CA" -p output/justCA.pdb
pyPcaunzip -c output/test5.pcz -t output/justCA.pdb -o output/test5.dcd
test ./traj_check.py output/justCA.pdb output/test5.dcd reference/test5.dcd 0.001 
fi
echo -n "Test 5b - BINPOS compression/DCD decompression: "
python -c 'import mdtraj' >& /dev/null
if [ $? -ne 0 ]; then
echo "mdtraj not available - SKIPPED"
else
pyPcazip --topology 2ozq.pdb -i 2ozq.binpos -o output/test5b.pcz --selection "name CA" -p output/justCA.pdb
pyPcaunzip -c output/test5b.pcz -t output/justCA.pdb -o output/test5b.dcd
test ./traj_check.py output/justCA.pdb output/test5b.dcd reference/test5b.dcd 0.001 
fi
echo -n "Test 6 - use of --fast option: "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd -o output/test6.pcz --selection "name CA"  -p output/justCA.pdb --fast
pyPcaunzip -c output/test6.pcz -t output/justCA.pdb -o output/test6.dcd
test ./traj_check.py output/justCA.pdb output/test6.dcd reference/test1.dcd 0.001 
echo -n "Test 7 - use of --lowmem option: "
pyPcazip --topology 2ozq.pdb -i 2ozq.dcd -o output/test7.pcz --selection "name CA"  -p output/justCA.pdb --lowmem
pyPcaunzip -c output/test7.pcz -t output/justCA.pdb -o output/test7.dcd
test ./traj_check.py output/justCA.pdb output/test7.dcd reference/test1.dcd 0.001 
