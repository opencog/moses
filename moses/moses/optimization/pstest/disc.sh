#!/bin/sh
make -C ../../../../build/ moses

MOSES="../../../../build/moses/moses/main/moses"
FOLDER="-i ../../../../examples/example-data"
ALGO="-a ps -l DEBUG"
GDB="gdb -ex \"break particle-swarm.cc:272\" -ex \"run\" -args "


# Disjunction bit+disc
#eval $GDB \
#$MOSES -H it $FOLDER/disjunction.csv $ALGO

# Column Labels bit+disc
#eval $GDB \
#$MOSES -H it $FOLDER/column-labels.csv -W1 -u or_them_all $ALGO

# Parity 3 bit+disc
#eval $GDB \
$MOSES -Hpa -k3 $ALGO

# Disjunction 3 bit+disc
#eval $GDB \
#$MOSES -Hdj -k3 $ALGO

# Multiplex 3 bit+disc
#eval $GDB \
#$MOSES -Hmux -k3 $ALGO

# Parity 5 bit+disc
#eval $GDB \
#$MOSES -Hpa -k5 $ALGO

