# Qubitserf
Distance finder for quantum codes.

Disclaimer: this code is still very much in early-development, but its main functionality is mostly considered to have been achieved.

## Build instructions
Call the script `compile.py`. You will find the generated files in the `build` folder, most importantly the executable `interface`.

## Use instructions for the interface
The program uses a meet in the middle algorithm by default which uses $\mathcal{O}(n^{\lfloor d/2 \rfloor})$ memory, where $n$ is the number of qubits and $d$ is the distance of the code. Add the stabilizers line by line as in the example below. The input is considered complete upon reading an empty line or `EOF`.
#### When the input is a CSS code
Run using the flag `--zx` for the program to output both the $Z$-distance and the $X$-distance.

If one wants to compute the distance using the Brouwer-Zimmerman algorithm (constant memory use, parallelizable) instead of the default, call the interface using the `--bz` flag; if you wish to parallelize the computation, use the `--threads` flag, followed by maximum number of threads to break the process into.
#### When the input is not a CSS code
There are no flags to use.

### Example interactions
In the following example, a large CSS code is provided as input and the program outputs `5 4`, meaning that the code has $Z$-distance $5$ and $X$-distance 4. It also lets the program break the task into 16 threads.
```
./interface --bz --zx --threads 16
ZZIIIZZIIIIIIIIIIIII
IXXIIIXXIIIIIIIIIIII
IIZZIIIZZIIIIIIIIIII
IIIXXIIIXXIIIIIIIIII
IIIIIXXIIIXXIIIIIIII
IIIIIIZZIIIZZIIIIIII
IIIIIIIXXIIIXXIIIIII
IIIIIIIIZZIIIZZIIIII
IIIIIIIIIIZZIIIZZIII
IIIIIIIIIIIXXIIIXXII
IIIIIIIIIIIIZZIIIZZI
IIIIIIIIIIIIIXXIIIXX
XXIIIIIIIIIIIIIIIIII
IIXXIIIIIIIIIIIIIIII
IIIIIIIIIIIIIIIXXIII
IIIIIIIIIIIIIIIIIXXI
IIIIIZIIIIZIIIIIIIII
IIIIZIIIIZIIIIIIIIII
IIIIIIIIIIIIIIZIIIIZ

5 4
```

In the following example, a large CSS code is provided as input and the program outputs `4`, meaning that the code has distance 4.
```
./interface
ZZIIIZZIIIIIIIIIIIII
IXXIIIXXIIIIIIIIIIII
IIZZIIIZZIIIIIIIIIII
IIIXXIIIXXIIIIIIIIII
IIIIIXXIIIXXIIIIIIII
IIIIIIZZIIIZZIIIIIII
IIIIIIIXXIIIXXIIIIII
IIIIIIIIZZIIIZZIIIII
IIIIIIIIIIZZIIIZZIII
IIIIIIIIIIIXXIIIXXII
IIIIIIIIIIIIZZIIIZZI
IIIIIIIIIIIIIXXIIIXX
XXIIIIIIIIIIIIIIIIII
IIXXIIIIIIIIIIIIIIII
IIIIIIIIIIIIIIIXXIII
IIIIIIIIIIIIIIIIIXXI
IIIIIZIIIIZIIIIIIIII
IIIIZIIIIZIIIIIIIIII
IIIIIIIIIIIIIIZIIIIZ

4
```

