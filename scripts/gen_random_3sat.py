import random, sys
n, m, seed, path = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
random.seed(seed)
with open(path, "w") as f:
    f.write(f"p cnf {n} {m}\n")
    for _ in range(m):
        clause = set()
        while len(clause) < 3:
            lit = random.randrange(1, n + 1)
            if random.random() < 0.5:
                lit = -lit
            clause.add(lit)
        f.write(" ".join(map(str, clause)) + " 0\n")
