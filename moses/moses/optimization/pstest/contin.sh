#!/bin/sh
make -j4 -C ../../../../build/

MOSES="../../../../build/moses/moses/main/moses"
FOLDER="-i ../../../../examples/example-data"
ALGO="-a hc -l DEBUG --contin-depth 2 -B 0 -E 1"
GDB="gdb -ex \"break hill-climbing.cc:108\" -ex \"run\" -args "

#eval $GDB \
#valgrind --tool=massif \
#valgrind --tool=callgrind \
#-B 0 -E 1

# Predicates all
#valgrind --tool=callgrind \
eval $GDB \
$MOSES -H it $FOLDER/predicates.csv -W1 -u pred -m1000 $ALGO

# Iris
#$MOSES $FOLDER/iris.data -u class -n sin -n log -n exp -n div -m2000 $ALGO

# Wine
#$MOSES $FOLDER/wine.data -u1 -n sin -n log -n exp -n div -m100 $ALGO

# Bank
#$MOSES -H it $FOLDER/bank.csv -W1 -u Q3 -n sin -n log -n exp -m200 $ALGO

# Polynomial factorization
#$MOSES -Hsr -k3 -m100 $ALGO
