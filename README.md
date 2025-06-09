# A Theoretical and Experimental Comparison of SAT‑Solving Algorithms

Minimal, header‑only C++20 implementation of three classical SAT algorithms  
(resolution, Davis–Putnam, Davis–Logemann–Loveland) together with three
branching heuristics (Random, MOMs, Jeroslow–Wang).

## Build (native)

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
