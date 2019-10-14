import pandas as pd
mv = pd.read_csv("meanVariances.csv", index_col=0)
o=pd.read_csv("O.dat", header=None)
mv.insert(3,'occurrence',o.values)
mv.to_csv("meanVariances.csv")
