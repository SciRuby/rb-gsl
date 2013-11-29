#!/bin/sh
gsl-randist 0 10000 cauchy 30 | ruby cauchy.rb -100 100 200
