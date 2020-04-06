#!/bin/sh

system=`uname -s`
nrandom=`date '+%d%H%M%S'`
mkdir $nrandom
cd $nrandom
case $system in
    Linux|LINUX|UNIX|Unix|Darwin|DARWIN)

        echo $nrandom > nseed_runtime
        ;;
    OSF*)
        echo $nrandom > nseed_runtime
        ;;
    *)
        echo "You are running operating system other than Linux/UNIX/OSF"
        echo "Modify this file first"
        exit 1;
        ;;
esac

# rm -f *.o
#make
mkdir ana
cp ../input.ampt ana/
cp ../input.ampt .
cp ../ampt .
cp ../makeAmptroot.C .
echo "#  AMPT started at " `date` > start.time
./ampt < nseed_runtime 
root -b -q makeAmptroot.C++
mv ana/test.root ../test_$nrandom.root
cd ..
rm -r $nrandom
