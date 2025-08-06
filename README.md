# Qubitserf

A command‑line tool for computing the minimum distance of quantum codes, with support for CSS codes and parallelized algorithms.

> **Disclaimer:** Qubitserf is in early development. While main functionality is implemented, expect edge cases and performance tuning to evolve.

---

## Features

* **Distance computation** using a meet‑in‑the‑middle algorithm (default) with $\mathcal{O}(n^{\lfloor d/2 \rfloor})$ memory, where:

  * $n$ = number of qubits
  * $d$ = code distance
* **Brouwer‑Zimmerman algorithm** (`--bz`): constant memory, parallelizable, ideal for large codes.
* **CSS‑specific** multi‑distance output (`--zx`): returns both $Z$-distance and $X$-distance.
* **Threaded computation** (`--threads <count>`): specify number of worker threads for parallel runs.
* **Verbose logging** (`--verbose` or `-v`): enable detailed progress output.
* **Flexible input**: read stabilizer generators line by line; terminate input on an empty line or EOF.

---

## Installation

1. **Clone the repository**:

   ```bash
   git clone https://github.com/yourusername/qubitserf.git
   cd qubitserf
   ```
2. **Build**:

   ```bash
   make
   make interface
   ```
3. **Artifacts**: The `build/` directory will contain:

   * `interface` : the main executable
   * Supporting libraries and object files

---

## Usage

```bash
./interface [OPTIONS]
```

1. **Prepare input**: enter stabilizer generators (one per line) in Pauli notation (`I`, `X`, `Y`, `Z`), for example:

   ```text
   ZZIIXZ...  
   XIIZZ...   
   ...        
   <empty line or Ctrl‑D>
   ```
2. **Run with flags** (all flags are optional):

   * `--bz`
     : Use the Brouwer‑Zimmerman algorithm instead of the default.
   * `--zx`
     : For CSS codes only, output `Z`‑ and `X`‑distance as two integers.
   * `--threads <N>`
     : Parallelize using up to `N` threads. Must be a positive integer.
   * `--verbose`, `-v`
     : Print detailed progress and diagnostics to stderr.

### Notes

* `--bz` and `--zx` **require** CSS codes. The program will error if used on non‑CSS input.
* Input terminates on an empty line or EOF (`Ctrl‑D` on UNIX, `Ctrl‑Z` on Windows).

---

## Examples

### Default (meet‑in‑the‑middle) on a CSS code

```bash
$ ./interface --zx  
ZZIIIZZIIIII
IXXIIIXXIIII
... (more lines)
<empty line>
5 4
```

Outputs `5 4`, the $Z$- and $X$-distances.

### Meet in the middle, parallelized, verbose

```bash
$./interface -v --threads 16
ZIIIIIIIZXYYIIYYXIIIII
XIIIIIIIXYZZIIZZYIIIII
IZIIIIIIXXYXYIZXIIIIII
IXIIIIIIYYZYZIXYIIIIII
IIZIIIIIIXXYXYIZXIIIII
IIXIIIIIIYYZYZIXYIIIII
IIIZIIIIXYYYYXXZXIIIII
IIIXIIIIYZZZZYYXYIIIII
IIIIZIIIXZXXYYYYXIIIII
IIIIXIIIYXYYZZZZYIIIII
IIIIIZIIXZIYXYXXIIIIII
IIIIIXIIYXIZYZYYIIIIII
IIIIIIZIIXZIYXYXXIIIII
IIIIIIXIIYXIZYZYYIIIII
IIIIIIIZXYYIIYYXZIIIII
IIIIIIIXYZZIIZZYXIIIII
IIIIIIIIIIIIIIIIIZIIII
IIIIIIIIIIIIIIIIIIZIII
IIIIIIIIIIIIIIIIIIIZII
IIIIIIIIIIIIIIIIIIIIZI
IIIIIIIIIIIIIIIIIIIIIZ

Distance bound: >1
Elapsed:[15ms] 
Distance bound: >2
Elapsed:[0ms] 
Distance bound: >3
Elapsed:[15ms] 
Distance bound: >4
Elapsed:[9ms] 
Distance bound: >5
Elapsed:[119ms] 
Distance bound: >6
Elapsed:[101ms] 
Distance: =7
Elapsed:[1.087s] 
```

### Non-CSS code distance

```bash
$ ./interface
ZZIIIZIIIZ
IZIZIZIZIZ
... (more lines)
<empty line>
4
```

Prints the single distance `4` for a non-CSS code.

