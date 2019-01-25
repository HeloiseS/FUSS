from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class GPfit2D(object): # Tested
    def __init__(self, input_file = None, param1 = None, param2 = None,
                 lnlik = None, fitlim1= None, fitlim2= None, res= None):

        if input_file is None:
            self.param1 = param1
            self.param2 = param2
            self.lnlik = lnlik

        else:
            self._read_input(input_file)

        self.fitlim1 = fitlim1
        self.fitlim2 = fitlim2
        self.res = res
        self.lnlik_pred = None
        self.fit_sigma = None

    def _read_input(self, input_file):
        """
        May need chaning as change output of model sampling function

        """
        data = pd.read_csv(input_file, names=['lnlik', 'x', 'y'])
        self.param1, self.param2 = data['x'].get_values(), data['y'].get_values()
        self.lnlik = data['lnlik'].get_values()

    def fit(self, plot = True, plot_save=False):
        X = np.atleast_2d([self.param1, self.param2]).T

        # need a mesh for surface/2D plot fit
        x1, x2 = np.meshgrid(np.linspace(self.fitlim1[0], self.fitlim1[1], self.res),
                             np.linspace(self.fitlim2[0], self.fitlim2[1], self.res))

        xx = np.vstack([x1.reshape(x1.size), x2.reshape(x2.size)]).T

        # Instanciate a Gaussian Process model
        kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=15)

        # Fit to data using Maximum Likelihood Estimation of the parameters
        gp.fit(X, self.lnlik.reshape(-1,1))

        # Make the prediction on the meshed x-axis (ask for MSE as well)
        self.lnlik_pred, self.fit_sigma = gp.predict(xx, return_std=True)

        ## Finding the maximum and location at max !!!
        lnlik_pred_reshaped = (self.lnlik_pred.T[0]).reshape(100,100)
        j,k = np.argwhere(lnlik_pred_reshaped == lnlik_pred_reshaped.max())[0]
        self.param1_best, self.param2_best, self.lnlik_best = x1[j,k], x2[j,k], lnlik_pred_reshaped[j,k]
        print "BEST param1, param2 and LNLIK"
        print self.param1_best, self.param2_best, self.lnlik_best

        if plot is True:
            fig = plt.figure(figsize=(10,10))
            #2D plot
            plt.scatter(x1, x2, c=lnlik_pred_reshaped, vmin=lnlik_pred_reshaped.min(), vmax=lnlik_pred_reshaped.max(), cmap = cm.magma_r)#gist_earth)
            plt.scatter(x1[j,k], x2[j,k], c='white')
            plt.show()

            #3D Figure
            fig = plt.figure(figsize=(10,10))
            ax = fig.gca(projection='3d')
            fitplot = ax.plot_surface(x1, x2, lnlik_pred_reshaped, alpha=0.5, cmap = cm.coolwarm, antialiased=True)
            fig.colorbar(fitplot,shrink=0.5 )
            ax.scatter(self.param1, self.param2, self.lnlik, c='k', zorder=1000)
            ax.scatter(self.param1_best, self.param2_best, self.lnlik_best, c='red')
            ax.set_xlabel('Center x')
            ax.set_ylabel('Center y')
            ax.set_zlabel('Center log likelihood')
            if plot_save is True:
                plt.savefig('lnlik_surface')
            plt.show()

        return x1, x2, self.lnlik_pred, self.fit_sigma