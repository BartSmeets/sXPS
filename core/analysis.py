import numpy as np

class data:
    def __init__(self, file, worksheet):
        '''
        Load data
        '''
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
        self.energy = self.df['Binding Energy (eV)'].to_numpy(dtype=np.float64)
        self.counts = self.df['Counts / s'].to_numpy(dtype=np.float64)

        return
    

    def baseline(self, lamb=1e6, ratio=1e-6, niter=10, mask=None):
                '''
                arPLS Baseline Correction 
                From Baek et al. Analyst (2014); DOI:10.1039/C4AN01061B
                '''
                # Initialisation
                from scipy import sparse
                from scipy.sparse.linalg import spsolve
                from numpy.linalg import norm

                ## Matrices
                y = self.counts
                N = len(y)
                diag = np.ones(N - 2)
                D = sparse.spdiags([diag, -2*diag, diag], [0, -1, -2], N, N-2)
                H = lamb * D.dot(D.T)
                w = np.ones(N)
                W = sparse.spdiags(w, 0, N, N)
                ## Loop
                count = 0 

                while True:
                    # Update weights
                    W.setdiag(w)
                    z = spsolve(W + H, W * y)

                    # Make d- for logistic function
                    d = y - z
                    dn = d[d < 0]
                    ## mean and std of d-
                    m = np.mean(dn)
                    s = np.std(dn)
                    ## Logistic function
                    k = 2   # logistic steepness
                    w_new = 1 / (1 + np.exp(k * (d - (2*s - m))/s))
                    
                    if mask is not None:
                        w_new[mask] = 1e-4

                    # Check conditions
                    count += 1 
                    if count > niter:
                        break
                    if norm(w - w_new)/norm(w) < ratio:
                        break
                    else:
                        w = w_new
                return z
    

class peak:
    def __init__(self):
        pass
     
    def __ulrik(x, E, F, a=0, b=0, m=0):
        '''
        a and b: asymmetry parameters
        F: FWHM
        E: peak position
        m: Lorentzian fraction
        '''


        def GL(x, F, E, m):
            return np.exp(-4*np.log(2)*(1-m)*(x-E)**2/(F**2))/(1 + 4*m*(x-E)**2/(F**2))
        def G(x, F, E):
            return np.exp(-4*np.log(2)*(x-E)**2/F**2)
        def w(a, b):
            return b*(0.7+0.3/(a+0.01))
        def AW(x, a, F, E):
            teller = 2*np.sqrt(np.log(2))*(x-E)
            noemer = F - a*2*np.sqrt(np.log(2))*(x-E)
            return np.exp(-(teller/noemer)**2)

        return np.where(x>E, GL(x, F, E, m), GL(x, F, E, m) + w(a,b)*(AW(x, a, F, E)-G(x, F, E)))
    

    def __call__(self, x, a=0, b=0, m=0):
        return self.amp*self.__ulrik(x, self.energy, self.FWHM, a, b, m)
    

class fit:
    def __init__(self, num):
        total = [peak()]*num
    
