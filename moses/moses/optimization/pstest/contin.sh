#!/bin/sh
make -C ../../../../build/ moses

MOSES="../../../../build/moses/moses/main/moses"
FOLDER="-i ../../../../examples/example-data"
ALGO="-a ps -l DEBUG"
GDB="gdb -ex \"break particle-swarm.cc:100\" -ex \"run\" -args "

# Disjunction


# Predicates all
#eval $GDB \
$MOSES -H it $FOLDER/predicates.csv -W1 -u pred -m2000 $ALGO

# Iris
#eval $GDB \
#$MOSES $FOLDER/iris.data -u class -n sin -n log -n exp -n div -m100 $ALGO

# Wine
#eval $GDB \
#$MOSES $FOLDER/wine.data -u1 -n sin -n log -n exp -n div -m100 $ALGO

# Bank
#eval $GDB \
#$MOSES -H it $FOLDER/bank.csv -W1 -u Q3 -n sin -n log -n exp -m100 $ALGO

# Polynomial factorization
#eval $GDB \
#$MOSES -Hsr -k3 -m100 $ALGO
