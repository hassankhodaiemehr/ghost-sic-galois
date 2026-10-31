#!/bin/bash
k=$1
timeout 1800 gp -q -s 2G <<GP
default(realprecision,850);
read("norms.gp"); read("in_$k.gp");
norms(L,32);
GP
