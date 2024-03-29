#!/bin/sh

IFS='_'
i=0
job_num=0
for x in $1; do
	if [ $i -eq 1 ]; then
		job_num=$x
	fi
	i=$((i+1))
done
system=`uname -s`
nrandom=`date '+%d%H%M%S'`$job_num"1"
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
echo $nrandom
echo "#  AMPT started at " `date` > start.time
echo "AMPT started at " `date`
./ampt < nseed_runtime 
echo "AMPT ended Root started at " `date`
root -b -q 'makeAmptroot.C++('\""$nrandom"\"')'
echo "Root ended at " `date`
mv ana/data_$nrandom.root ../test_$nrandom.root
cd ..
rm -r $nrandom
