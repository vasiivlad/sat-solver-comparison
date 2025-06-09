import pandas as pd
import matplotlib.pyplot as plt, sys
df = pd.read_csv(sys.argv[1])
plt.plot(df['vars'], df['time_sec'])
plt.xlabel("variables")
plt.ylabel("seconds")
plt.title("DPLL (JW)")
plt.savefig(sys.argv[2])
