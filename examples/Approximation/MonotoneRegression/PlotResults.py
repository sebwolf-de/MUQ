import matplotlib.pyplot as plt
import h5py
import numpy as np

# Read the original data
fin = h5py.File('data/LTA01.h5')
strainData = np.array( fin['/Strain'] )
stressData = np.array( fin['/Stress'] )

# Read the predictions
fin = h5py.File('results/StressPredictions.h5')
strainPred = np.array( fin['/Strain'] )
stressPred = np.array( fin['/Stress'] )


plt.plot(strainData, stressData, '.', label='Observations')
plt.plot(strainPred, stressPred, label='Monotone Fit')
plt.show()
