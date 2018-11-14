from __future__ import division
import numpy as np
import cPickle as pickle
import yaml 
import pandas as pd
from matplotlib.patches import Ellipse
import matplotlib
import matplotlib.pyplot as plt


def target_values(q_data, u_data, delta_q_data, delta_u_data, ksi_q, ksi_u, delta_ksi_q, delta_ksi_u):
    """
    Finds the target values to put into a model given the target data values and the offsets of a particular model. 
    """
    target_q = q_data + ksi_q
    target_u = u_data + ksi_u 

    delta_target_q = np.sqrt(delta_q_data**2 + delta_ksi_q**2)
    delta_target_u = np.sqrt(delta_u_data**2 + delta_ksi_u**2)

    return target_q, target_u, delta_target_q, delta_target_u

    

def offsets( q_cont, u_cont, delta_q_cont, delta_u_cont, q_model, u_model):
    """
    Offset between the target data continuum q and u and a given model.    
    This offset is applied to target data values before they are modelled. 
    """
    ksi_q = q_model - q_cont
    ksi_u = u_model - u_cont

    delta_ksi_q = delta_q_cont
    delta_ksi_u = delta_u_cont


    return ksi_q, ksi_u, delta_ksi_q, delta_ksi_u



def lnlikelihood(model, data, error): # Tested
    """
    Find the Log likelihood given gaussian data

    Parameters
    ----------
    model : float
        Value found by the model
    data : float
        Value from the data tha we are trying to reproduce or TARGET value.
    error : float
        Error on the data value

    Returns
    -------
        The log likelihood as a float

    """
    return -0.5 * ( np.log(2*np.pi*error**2) + ((model-data)/error)**2 )


def read_yaml(input_yaml): # Tested
    """
    Reads in the parameters from a yaml file and creates a list of data frames, each containing the parameters for each run given in the yaml config file. 
    """
    with open(input_yaml, 'r') as stream:
        try:
            parameters = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    df_list=[]
    i = 0
    for run in parameters:
        param_df = pd.DataFrame(columns = parameters[run].keys())
        param_df.loc[0] = parameters[run].values()
        df_list.append(param_df)
        i+=1
    print "Created ",i," DataFrames containing the parameters in ", input_yaml, "\n"
        
    return df_list
        


def read_pkl(input_pkl): # Tested
    """
    Reads in and returns photosphere (or any) object from a pickle file
    """

    assert isinstance(input_pkl, str), "The input file name should be a string"

    try:
        with open(input_pkl, 'rb') as input:
            photosphere = pickle.load(input)
        print "Object imported from ", input_pkl
        return photosphere
    except IOError:
        print "File "+ input_pkl + " not found"


def list_params(lim_par1, lim_par2=None, lim_par3=None, step1=0.1, step2=0.1, step3=0.1): # Tested
    """
    Creates a list is parameters to try with the toy model

    Parameters
    ----------
    lim_par1 : tuple or list of 2 elements
        limits for the range of parameter 1 (upper limit non inclusive)
    lim_par2 : tuple or list of 2 elements
        limits for the range of parameter 2 (upper limit non inclusive)
    lim_par3 : tuple or list of 2 elements
        limits for the range of parameter 3 (upper limit non inclusive)
    step : float
        Increment for the ranges of parameters created. Same step is used for all parameter ranges.

    Returns
    -------
        List of lists containing the parameters.

    """
    params_ls = []

    if lim_par1 is not None and lim_par2 is None:
        par1ls = np.arange(lim_par1[0], lim_par1[1], step1)
        for x in par1ls:
            params_ls.append([x])

    elif lim_par1 is not None and lim_par2 is not None and lim_par3 is None:
        par1ls = np.arange(lim_par1[0], lim_par1[1], step1)
        par2ls = np.arange(lim_par2[0], lim_par2[1], step2)
        for x in par1ls:
            for y in par2ls:
                params_ls.append([x, y])

    else:
        par1ls = np.arange(lim_par1[0], lim_par1[1], step1)
        par2ls = np.arange(lim_par2[0], lim_par2[1], step2)
        par3ls = np.arange(lim_par3[0], lim_par3[1], step3)
        for x in par1ls:
            for y in par2ls:
                for z in par3ls:
                    params_ls.append([x, y, z])

    return params_ls


class ModelGrid(object):
    """
    Object containing the grid of models created by our sampler. 
    
    Parameters
    ----------
    grid : str or dataframe
        grid of models
        
    Notes
    -----
    Expected columns: [X	Y	R	Q	U	Fl/Fc	QCONT	UCONT]
    
    Attributes
    ----------
    grid : pandas.DataFrame
        The main grid attribute. It will undergo the tranformations by the methods
    grid_reset : pandas.DataFrame
        Attribute to save the initial grid into so it's easy to get back to. 
        Used by reset() to turn self.grid back to its original version
    grid_uncrop : None or pandas Data.Frame
        Uncropped grid (but if a continuum transformation has been applied this grid 
        will have it too, unlike grid_reset.
    targets : None or list or array (3 values)
        Target values for q, u and flux ratio. 
        None until crop_to_Xsigma() is called
    targets_r : None or list or array (3 values)
        Errors on the target values for q, u and flux ratio. 
        None until crop_to_Xsigma() is called
      
    Methods
    -------
    match_obs_cont()
    crop_to_Xsigma()
    plot()   
    reset()
    uncrop()
    
    """
    def __init__(self, grid):
        """
        Parameters
        ----------
        grid : str or dataframe
            grid of models
            
        Notes
        -----
        Expected columns: [X	Y	R	Q	U	Fl/Fc	QCONT	UCONT]
        
        
        """
        # imports the main grid dataframe
        if isinstance(grid, pd.core.frame.DataFrame): self.grid = grid
        elif isinstance(grid, str): self.grid = pd.read_csv(grid)
        
        # creates the "reset" grid. Used by the method "reset"
        self.grid_reset = self.grid
        
        # grid_uncrop will be set by "crop_to_Xsigma"
        # it allows to go back to uncropped grid without reseting entierly 
        self.grid_uncrop = None
        
        # will be used to store the target values of the STokes param and flux ratio
        # and the associated errors
        self.targets = None
        self.targets_r = None
 
    def match_obs_cont(self, cont, cont_r ):
        """
        Method used to transform the modelled Stokes parameters to match the observed continuum
        
        Parameters
        ----------
        cont : array, list or tuple - 2 values
            [ observed q_cont, observed u_cont]
        cont_r : array, list or tuple - 2 values
            [error on  osberved q_cont, error on observed u_cont]
        
        """
        ksi_q, ksi_u, delta_ksi_q, delta_ksi_u = offsets(cont[0], cont[1], cont_r[0], cont_r[1], self.grid_reset.QCONT[0], self.grid_reset.UCONT[0])

        self.grid.Q, self.grid.QCONT = self.grid_reset.Q - ksi_q, self.grid_reset.QCONT - ksi_q
        self.grid.U, self.grid.UCONT = self.grid_reset.U - ksi_u, self.grid_reset.UCONT - ksi_u
        
        self.grid_uncrop = self.grid

    def crop_to_Xsigma(self, targets, targets_r, Xsigma = 3):  
        """
        Method used to crop the grid to only contain models with Stokes parameters within X sigma of the target (observed) values
        
        Parameters
        ----------
        targets : array, list or tuple - 3 values
            [q_target, u_target, flux ratio target]
        targets_r : array, list or tuple - 3 values
            [error q_target, error u_target, error flux ratio target]
        Xsigma : int or float, optional
            Number of sigmas to use to crop. Default is 3.
            
        """
        if self.grid_uncrop is None: self.grid_uncrop = self.grid
        
        self.targets  = [float(tar) for tar in targets] # this just to make sure I have a list of floats
        self.targets_r = [float(tar_r) for tar_r in targets_r]
        
        qlim_low = self.targets[0] - Xsigma*self.targets_r[0]
        qlim_high = self.targets[0] + Xsigma*self.targets_r[0]

        ulim_low = self.targets[1] - Xsigma*self.targets_r[1]
        ulim_high = self.targets[1] + Xsigma*self.targets_r[1]

        flim_low = self.targets[2] - Xsigma*self.targets_r[2]
        flim_high = self.targets[2] + Xsigma*self.targets_r[2]
        
        self.grid = self.grid_uncrop
        self.grid = self.grid[(self.grid['Q'] < qlim_high) & (self.grid['Q'] > qlim_low) 
                    & (self.grid['U'] > ulim_low) & (self.grid['U'] < ulim_high) 
                    & (self.grid['Fl/Fc'] > flim_low)& (self.grid['Fl/Fc'] < flim_high)]  
        
        
    def plot(self, axis_ratio,transparency=0.15, phot_color = 'grey', cmap = 'Purples', outputfile = None, skip=None, lim=[-1,1]):    
        """
        Plotting method to visualise absorption regions on top of photosphere
        
        Parameters
        ----------
        axis_ratio : float
            axis ratio of the photosphere of the models in the grid
        transparency : float, optional
            transparency (alpha parameter) of individual absorption regions. Default = 0.15.
        phot_color : str (matplotlib color), optional
            color of the photosphere, Default is 'grey'
        cmap : matplotlib color map, optional
            colormap used to show the absorption regions, Default = 'Purples_r'
        outputfile : str, optional
             name of output file for the plot i fyou want to save the plot. Default is None (WILL NOT SAVE)
        
        """
        cmap = plt.cm.get_cmap(cmap)
        self.grid['Qres'] = self.targets[0] - self.grid.Q.values
        self.grid['Ures'] = self.targets[1] - self.grid.U.values
        self.grid['Fres'] = self.targets[2] - self.grid['Fl/Fc'].values
        
        self.grid['TOTres'] = abs(self.grid.Qres + self.grid.Ures + self.grid.Fres)        

        sorted_df = self.grid.sort_values('TOTres', ascending=False)
        
        if axis_ratio > 1 : axis_ratio = 1/axis_ratio
        
        photosphere = Ellipse(xy=(0,0), width = 2*0.99, height=2* axis_ratio*0.99)
        
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, aspect='equal')
        ax.add_artist(photosphere)
        photosphere.set_facecolor(phot_color)
        ax.set_ylim(lim)
        ax.set_xlim(lim)


        norm = matplotlib.colors.Normalize(vmin=sorted_df['TOTres'].min(), vmax=sorted_df['TOTres'].max())

        for X, Y, R, C, i in zip(sorted_df.X, sorted_df.Y, sorted_df.R, sorted_df['TOTres'], range(len(sorted_df.Y))):
            #print norm(C)
            if skip == None:
                circle = plt.Circle((X, Y), R, color=cmap(norm(C)), alpha=transparency)
                ax.add_artist(circle)
            elif skip != 0:
                if i % skip == 0:
                    circle = plt.Circle((X, Y), R, color=cmap(norm(C)), alpha=transparency)
                    ax.add_artist(circle)
                    
                if i == range(len(sorted_df.Y))[0] or  i == range(len(sorted_df.Y))[-1]:
                    circle = plt.Circle((X, Y), R, color=cmap(norm(C)), alpha=transparency)
                    ax.add_artist(circle)
                    

        ax.axvline(0, ls='--', c='k')
        ax.axhline(0, ls='--', c='k')

        circlebest = plt.Circle((sorted_df.X.values[-1], sorted_df.Y.values[-1]), sorted_df.R.values[-1], ls='--')

        ax.add_artist(circlebest)

        circlebest.set_facecolor('None')
        circlebest.set_edgecolor('k')
        
        if outputfile is not None:
            plt.savefig(outputfile)
        
        return fig, ax
        
    def uncrop(self):
        """
        If the grid of models was cropped, reverses to uncropped grid. 
        Will conserve any transformation of the model continuua through "mathc_obs_cont"
        """
        if self.grid_uncrop is not None: self.grid = self.grid_uncrop
        
    def reset(self):
        """
        Resets the grid to the version given on initialisation.
        """
        self.grid = self.grid_reset
        
