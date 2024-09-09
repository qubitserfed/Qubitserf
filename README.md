# Qubitserf
Distance calculator for CSS codes based on the Brouwer-Zimmerman algorithm.

Disclaimer: this code is still very much in early-development, but its main functionality is mostly considered to have been achieved.

## Build instructions
Call the script `compile.py`. You will find the generated files in the `build` folder, most importantly the executable `interface`.

## Use instructions for the interface
Run using the flag `--zx` for the program to output both the $Z$-distance and the $X$-distance. Run using the `--multithreaded` flag to enable the use of multithreading. Add the stabilizers line by line as in the example below. The input is considered complete upon reading an empty line or `EOF`.

### Example interactions
In the following example, a large CSS code is provided as input and the program outputs `5 4`, meaning that the code has $Z$-distance $5$ and $X$-distance 4.
```
./interface --zx
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


## To do:
* Write a python wrapper
* Implement the dynamic Brouwer-Zimmerman algorithm
* Enable GPU acceleration for the exponential part
* \$\$\$\$\$\$
