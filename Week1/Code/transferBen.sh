#!/bin/sh

NDOUBLE="1 4 10 100 1000 10000 100000 1000000 2000000 5000000 10000000"
FILE="time.txt"
echo "Start"



/bin/rm -f $FILE

for nb in $NDOUBLE
do
	echo "Transfer $nb doubles"
	mpiexec -n 2 ./benmark $nb | grep -v CPU >> $FILE
done

exit 0
