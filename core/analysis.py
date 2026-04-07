import numpy as np
from lmfitxps import models
import lmfit

class Data:
    def __init__(self, file, worksheet):
        '''
        Load data
        '''
        self.__load_data(file, worksheet)

        # Prepare fit
        self.models = {'bg': models.ShirleyBG(prefix='shirley_')}   # Add background
        self.pars = lmfit.Parameters()
        ## Background parameters
        self.pars.add('shirley_k', value=0.002, min=0, max=5)
        self.pars.add('shirley_const', value=self.counts[-1], min=0, max=2*np.max(self.counts))


    def __load_data(self, file, worksheet):
        import pandas as pd
        excel = file
        df = pd.read_excel(excel, worksheet)

        self.original_directory = df.iloc[2, 0]
        # Extract data from directory
        split = self.original_directory.split("\\")
        self.scan = split[-1][:-4]  # first select split, then remove the file extension
        self.xray = split[-3]

        # Clean and rename DataFrame
        df = df.drop(index=range(10)).reset_index(drop=True)
        df = df.drop(columns=[df.columns[1], df.columns[3]])
        df.columns = ['Binding Energy (eV)', 'Counts / s']
        
        # Convert to numpy
        self.energy = df['Binding Energy (eV)'].to_numpy(dtype=np.float64)
        self.counts = df['Counts / s'].to_numpy(dtype=np.float64)


    def add_singlet(self, key, pos):
        self.models[key] = models.ConvGaussianDoniachSinglett(prefix=f'{key}_')

        self.pars.add(f'{key}_center', value=pos, min=np.min(self.energy), max=np.max(self.energy))
        
        # Defaults
        self.pars.add(f'{key}_amplitude', value=np.max(self.counts), min=0.0)
        self.pars.add(f'{key}_sigma', value=0.2, min=0.05, max=3)
        self.pars.add(f'{key}_gaussian_sigma', value=0.1, min=1e-3)
        self.pars.add(f'{key}_gamma', value=0.0, min=1e-3, max=1.0)


    def add_doublet(self, key, pos, delta):
        self.models[key] = models.ConvGaussianDoniachDublett(prefix=f'{key}_')

        self.pars.add(f'{key}_center', value=pos, min=np.min(self.energy), max=np.max(self.energy))
        self.pars.add(f'{key}_soc', value=delta, vary=False)
        
        # Defaults
        self.pars.add(f'{key}_amplitude', value=np.max(self.counts), min=0.0)
        self.pars.add(f'{key}_sigma', value=0.2, min=0.05)
        self.pars.add(f'{key}_gaussian_sigma', value=0.1, min=1e-3)
        self.pars.add(f'{key}_gamma', value=1e-3, min=1e-3, max=1.0)
        self.pars.add(f'{key}_height_ratio', value=0.5, min=1e-3, max=1.0)
        self.pars.add(f'{key}_fct_coster_kronig', value=1.0, vary=False)
        

    def fit(self):
        return FitResult(self)
    

class FitResult:
    def __init__(self, data):
        self.data = data
        self.result = self.__run_fit()
        self.comps = self.result.eval_components()
        
    
    def __run_fit(self):
        data = self.data
        # Define model from dictionary
        for _, obj in data.models.items():
            try:
                mod += obj
            except UnboundLocalError:
                mod = obj
        
        # Run fit
        x = data.energy
        y = data.counts
        out = mod.fit(y, data.pars, y=y, x=x, weights=1 /(np.sqrt(y)))
        return out
    

    def plot(self):
        import matplotlib.pyplot as plt
        import toml

        settings = toml.load('core\\plot_settings.toml')

        data = self.data
        comps = self.comps
        bg = comps['shirley_']
        plt.plot(data.energy, bg, c='black', label='Shirley background')

        for (key, value), (_, colour) in zip(comps.items(), settings['colours']['core'].items()):
            if key == 'shirley_':
                continue
            else:
                plt.plot(data.energy, bg+value, label=key[:-1], color=colour)
                plt.fill_between(data.energy, bg, bg + value, alpha=0.3, color=colour)
