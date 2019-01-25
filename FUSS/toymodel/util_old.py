from __future__ import division
import numpy as np
import cPickle as pickle
import yaml 
import pandas as pd


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
