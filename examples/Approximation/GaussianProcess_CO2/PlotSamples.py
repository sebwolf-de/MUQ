import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



fileName = 'Samps.csv'

df = pd.read_csv(fileName)


plt.plot(df['Points'], df['Concentration'])
plt.show()
