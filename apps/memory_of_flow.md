# Measured memory requirements of flow

A single FlowContainer object requires ~16*L^6 bytes of memory (+something of the order of L^2, which is negligible).
We notice that during the execution of the progam, around 92.5 times that number is allocated.
The numbers appear indepedent on the thread number (tested with 1 and 8).
Measuring yields

L=8  -> 370 MB
L=10 -> 1.4 GB
L=12 -> 4.1 GB

extrapolating yields

L=14 -> 10.4 GB
L=16 -> 23.1 GB
L=18 -> 46.9 GB
L=20 -> 88.2 GB

# Runtimes

It appears that the runtime scales as L^10.
The times themselves of course heavily depend on the available resources.