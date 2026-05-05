import numpy as np
from lmfitxps import models
import lmfit

class Data:
    '''
    Data object
    '''
    def __init__(self, file:str, worksheet:str):
        '''
        Load data from excel sheet

        # Arguments

        **file**: file name of the data (excel sheet) that must be loaded

        **worksheet**: worksheet to load

        # Defines

        **models**: dictionary containing the model elements (initially background; peaks can be added)

        **pars**: parameter object
        '''

        self.__load_data(file, worksheet)

        # Prepare fit
        self.models = {'bg': models.ShirleyBG(prefix='shirley_')}   # Add background
        self.pars = lmfit.Parameters()
        ## Background parameters
        self.pars.add('shirley_k', value=0.002, min=0, max=5)
        self.pars.add('shirley_const', value=self.counts[-1], min=0, max=2*np.max(self.counts))


    def __load_data(self, file, worksheet):
        '''
        Load the data

        # Arguments

        **file**: file name of the data (excel sheet) that must be loaded

        **worksheet**: worksheet to load

        # Defines

        **original_directory**: original data directory

        **scan**: scan name

        **xray**: name of the x-ray source

        **energy**: array containing the binding energies

        **counts**: array containting the counts ont he detector

        '''


        import pandas as pd

        # read file
        excel = file
        df = pd.read_excel(excel, worksheet)
        
        # Extract data from directory
        self.original_directory = df.iloc[2, 0]
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
        '''
        Adds singlet peak to models

        # Arguments
        **key**: key of the peak

        **pos**: expected position of the peak

        # Defines
        **models[key]**: adds peak to model

        **pars[{key}_center]**: position of the peak

        **pars[{key}_amplitude]**: amplitude of the peak

        **pars[{key}_sigma]**: broadening of the Doniach-Sunjic 

        **pars[{key}_gaussian_sigma]**: broadening of the gaussian kernel

        **pars[{key}_gamma]**: asymmetry of the Doniach-Sunjic 
        '''


        self.models[key] = models.ConvGaussianDoniachSinglett(prefix=f'{key}_')

        self.pars.add(f'{key}_center', value=pos, min=np.min(self.energy), max=np.max(self.energy))
        
        # Defaults
        self.pars.add(f'{key}_amplitude', value=np.max(self.counts), min=0.0)
        self.pars.add(f'{key}_sigma', value=0.2, min=1e-6)
        self.pars.add(f'{key}_gaussian_sigma', value=0.1, min=1e-6)
        self.pars.add(f'{key}_gamma', value=0.0, min=0, max=1.0)


    def add_doublet(self, key, pos, delta):
        '''
        Adds doublet peaks to models

        # Arguments
        **key**: key of the main peak

        **pos**: expected position of the main peak

        **delta**: separation of the soublet

        # Defines

        **models[key]**: adds peak to model

        **pars[{key}_center]**: position of the main peak

        **pars[{key}_soc]**: separation of the peaks

        **pars[{key}_amplitude]**: amplitude of the main peak

        **pars[{key}_sigma]**: broadening of the Doniach-Sunjic 

        **pars[{key}_gaussian_sigma]**: broadening of the gaussian kernel

        **pars[{key}_gamma]**: asymmetry of the Doniach-Sunjic 

        **pars[{key}_height_ratio]**: ratio of the amplitude of the two peaks

        **pars[{key}_fct_coster_kronig]**: ratio of the widths of the two peaks 
        '''


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
        '''
        Performs fit
        '''

        return FitResult(self)
    

class FitResult:
    '''
    Fit Result Class
    '''
    def __init__(self, data):
        '''
        Runs the fit on the procided data object

        # Arguments

        **data**: data object

        # Defines

        **data**: stores a copy of data within this class

        **result**: lmfit result object

        **comps**: dictionary of the different components of the fitting results, i.e., models
        '''


        self.data = data
        self.result = self.__run_fit()
        self.comps = self.result.eval_components()
        
    
    def __run_fit(self):
        '''
        Runs the fit

        # Returns

        **out**: lmfit result object
        '''


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
        weights = 1 /(np.sqrt(y))
        out = mod.fit(y, data.pars, y=y, x=x, weights=weights)
        return out


    def display(self):
        '''
        Print the fit report
        '''
        print(self.result.fit_report())


    def save(self, filename):
        '''
        Save the fit report to a text file

        # Argument
        **filename**: file name of the target file
        '''
        with open(filename, 'w') as f:
            f.write(self.result.fit_report())
