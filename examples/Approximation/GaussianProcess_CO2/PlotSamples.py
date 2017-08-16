import h5py
import matplotlib.pyplot as plt

# Read origional data
f1 = h5py.File('data/MaunaLoaCO2.h5')
date1 = f1['/Weekly/Dates']
conc1 = f1['/Weekly/Concentrations']

# Read the prediction results
f2 = h5py.File('results/CO2_Prediction.h5')
date2 = f2['/Predict/Dates']
conc2 = f2['/Predict/Concentrations']

# Define HDF attributes
f2['/Predict'].attrs['Source'] = 'Results from GaussianProcess_CO2.cpp'
f2['/Predict/Concentrations'].attrs['Units'] = 'CO2 mole fraction as parts per million (ppm)'
f2['/Predict/Dates'].attrs['Units'] = '(year,month,day) as a decimal'

# Plot them together
plt.plot(date1, conc1)
plt.plot(date2[0,:], conc2[0,:])
plt.show()
