#!/usr/bin/python

import time
import os

MOSES = "../../../../build/moses/moses/main/moses"
FOLDER = "-i ../../../../examples/example-data"
ALGO = "-a ps -l DEBUG -B 0 -E 1 --contin-depth "
VALG = "valgrind --tool=callgrind"

script = MOSES + " -H it " + FOLDER + \
    "/predicates.csv -W1 -u pred -m1000 " + ALGO

for depth in range(2,33):
    mean = 0.0
    for seed in range(1,11):
        print "DEPTH: %d | SEED: %d" % (depth, seed)
        tempscript = script + `depth` + " -r " + `seed` + " >> time.out"

        start = time.time()
        os.system(tempscript)
        end = time.time()

        diff = end - start
        print "Time: %f" % diff
        mean = mean + diff

    print "Mean: %f" % (mean / 10)
