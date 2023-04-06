import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import pymech.neksuite as nek
from scipy.interpolate import griddata

class NekFlowConfig():
    def __init__(self, field_path, **kwargs):
        self.field_path = field_path
        self.load_mesh()
        self.interpolation_setup(**kwargs)
        
    def filename(self, t_idx):
        return f'{self.prefix}0.f{t_idx:05d}'

    def load_mesh(self):
        # Read base flow field for grid points
        field = nek.readnek(f'{self.field_path}/{self.filename(1)}')
        self.nel = len(field.elem) # Number of spectral elements
        self.nGLL = field.elem[0].pos.shape[-1] # Order of the spectral mesh
        self.n = self.nel*self.nGLL**2

        # Load cell centers
        self.Cx = np.array([field.elem[i].pos[0, 0, j, k]
                       for i in range(self.nel) for j in range(self.nGLL) for k in range(self.nGLL)])
        self.Cy = np.array([field.elem[i].pos[1, 0, j, k]
                       for i in range(self.nel) for j in range(self.nGLL) for k in range(self.nGLL)])

        dOmega = np.loadtxt(f'{self.field_path}/mass_matrix.dat')
        self.dOmega = np.concatenate((dOmega, dOmega))  # Double for both velocity components

        self.dot = lambda a, b: np.vdot(self.dOmega*a, b)

    def mask_domain(self, mask):
        """Set dOmega[mask] = 0"""
        dOmega = self.dOmega[:self.n]
        dOmega[mask] = 0
        print(f"Integration volume: {sum(dOmega)}")
        self.dOmega = np.concatenate((dOmega, dOmega))

    def get_coeffs(self, mode_path=None, fname='coeffs.dat'):
        if mode_path is None:
            mode_path = self.field_path
        data = np.loadtxt(f'{mode_path}/{fname}')
        t = data[:, 0]
        a = data[:, 1:]
        return t, a

    def load_modes(self, r, mode_path=None, vort=True):
        n = self.n
        if mode_path is None:
            mode_path = self.field_path

        # Load all modes and gradients  (first will usually be mean flow)
        self.U = np.zeros((2*n, r))
        for i in range(r):
            self.U[:, i] = self.get_velocity(f'{mode_path}/{self.prefix}0.f{i+1:05d}')

        if vort:
            self.vort = np.zeros((n, r))
            for i in range(r):
                self.vort[:, i] = self.get_vorticity(f'{mode_path}/{self.prefix}0.f{i+1:05d}')

    def load_modes_and_grads(self, r, mode_path=None, vort=False, pres_grad=False):
        if mode_path is None:
            mode_path = self.field_path
        self.load_modes(r, mode_path=mode_path, vort=vort)
        n = self.n

        self.gradUx = np.zeros(self.U.shape)  # [ux, uy]
        self.gradUy = np.zeros(self.U.shape)  # [vx, vy]
        lapUx = np.zeros(self.U.shape)   # [uxx, uyy]
        lapUy = np.zeros(self.U.shape)   # [vxx, vyy]

        for i in range(r):
            self.gradUx[:, i] = self.get_velocity(f'{mode_path}/du_{self.prefix}0.f{i+1:05d}')
            self.gradUy[:, i] = self.get_velocity(f'{mode_path}/dv_{self.prefix}0.f{i+1:05d}')
            lapUx[:, i] = self.get_velocity(f'{mode_path}/ddu{self.prefix}0.f{i+1:05d}')
            lapUy[:, i] = self.get_velocity(f'{mode_path}/ddv{self.prefix}0.f{i+1:05d}')
        
        # Combine Laplacian terms
        self.lapU = np.zeros(self.U.shape)
        self.lapU[:n, :] = lapUx[:n, :] + lapUx[n:, :]  # uxx + uyy
        self.lapU[n:, :] = lapUy[:n, :] + lapUy[n:, :]  # vxx + vyy
        
        if pres_grad:
            self.gradP = np.zeros(self.U.shape)
            for i in range(r):
                self.gradP[:, i] = self.get_velocity(f'{mode_path}/dp_{self.prefix}0.f{i+1:05d}')

    def load_base_flows_and_grads(self, nb, base_path=None):
        if base_path is None:
            assert(hasattr(self, 'base_path'))
            base_path = self.base_path
        n = self.n

        # Load base flow for each input directory
        self.UB = np.zeros((2*n, nb))
        self.gradUBx = np.zeros(self.UB.shape)
        self.gradUBy = np.zeros(self.UB.shape)
        for i in range(nb):
            self.UB[:, i] = self.get_velocity(f'{base_path}/{self.prefix}0.f{i+1:05d}')
            self.gradUBx[:, i] = self.get_velocity(f'{base_path}/du_{self.prefix}0.f{i+1:05d}')
            self.gradUBy[:, i] = self.get_velocity(f'{base_path}/dv_{self.prefix}0.f{i+1:05d}')

        pass

    def interpolation_setup(self, **kwargs):
        self.XX, self.YY = None, None
    
    def interp(self, field, method='cubic'):
        """
        field - 1D array of cell values
        Cx, Cy - cell x-y values
        X, Y - meshgrid x-y values
        grid - if exists, should be an ngrid-dim logical that will be set to zer
        """
        ngrid = len(self.XX.flatten())
        grid_field = np.squeeze( np.reshape(
            griddata((self.Cx, self.Cy), field, (self.XX, self.YY), method=method),
            (ngrid, 1)
        ))
        if hasattr(self, 'interp_mask') and self.interp_mask is not None:
            grid_field[self.interp_mask] = 0
        return grid_field
        
    def plot_field(self, field, clim=[-5, 5], levels=None):
        pass
        
    def get_velocity(self, file):
        field = nek.readnek(file)
        u = np.array([field.elem[i].vel[0, 0, j, k]
                   for i in range(self.nel) for j in range(self.nGLL) for k in range(self.nGLL)])
        v = np.array([field.elem[i].vel[1, 0, j, k]
                   for i in range(self.nel) for j in range(self.nGLL) for k in range(self.nGLL)])
        return np.concatenate((u, v))

    def get_vorticity(self, file):
        field = nek.readnek(file)
        vort = np.array([field.elem[i].temp[0, 0, j, k]
                   for i in range(self.nel) for j in range(self.nGLL) for k in range(self.nGLL)])
        return vort

class LidDrivenCavity(NekFlowConfig):
    def __init__(self, field_path, base_path=None, **kwargs):
        self.prefix = 'cav'
        self.base_path = base_path
        super().__init__(field_path, **kwargs)

    def interpolation_setup(self, nx=200, ny=200, **kwargs):
        self.nx, self.ny = nx, ny
        self.x = np.linspace(-1, 1, nx)
        self.y = np.linspace(-1, 1, ny)
        self.XX, self.YY = np.meshgrid(self.x, self.y)
        
    def plot_field(self, field, clim=[-5, 5], levels=None):
        if levels is None: levels=np.linspace(*clim, 10)
        plt.contourf(self.x, self.y, np.reshape(field, [self.nx, self.ny], order='F').T,
                     cmap='RdBu', levels=levels, vmin=clim[0], vmax=clim[1], extend='both')

class ShearDrivenCavity(NekFlowConfig):
    def __init__(self, field_path, base_path=None, **kwargs):
        self.prefix = 'cav'
        self.base_path = base_path
        super().__init__(field_path, **kwargs)
        
    def interpolation_setup(self, nx=400, ny=200, L0=-1, Lx=2.5, Ly=0.5, **kwargs):
        self.nx, self.ny = nx, ny
        self.x = np.linspace(L0, Lx, nx)
        self.y = np.linspace(-1, Ly, ny)
        self.XX, self.YY = np.meshgrid(self.x, self.y)
        
    def plot_field(self, field, clim=None, levels=None):
        if clim is None: clim = (-100, 100)
        if levels is None: levels=np.linspace(*clim, 10)
        plt.contourf(self.x, self.y, np.reshape(field, [self.nx, self.ny], order='F').T,
                     cmap='RdBu', levels=levels, vmin=clim[0], vmax=clim[1], extend='both')

        # Plot walls
        wall_1 = plt.Rectangle((0, 0), (-1.2), (-1), fill='true', edgecolor='k', facecolor='gray')
        wall_2 = plt.Rectangle((1, 0), (2), (-1), fill='true', edgecolor='k', facecolor='gray')
        plt.gcf().gca().add_artist(wall_1)
        plt.gcf().gca().add_artist(wall_2)

        plt.gcf().gca().set_xlim([-0.5, 2.5])
        plt.gcf().gca().set_position([0, 0, 1, 1])

class MixingLayer(NekFlowConfig):
    def __init__(self, field_path, **kwargs):
        self.prefix = 'mix'
        super().__init__(field_path, **kwargs)

    def interpolation_setup(self, nx=400, ny=120, Lx=250, Ly=10, L0=0, **kwargs):
        self.nx, self.ny = nx, ny
        self.x = np.linspace(L0, Lx, nx)
        self.y = np.linspace(-Ly, Ly, ny)
        self.XX, self.YY = np.meshgrid(self.x, self.y)

    def plot_field(self, field, clim=None, levels=None, cm=None, black_contours=True):
        if clim is None: clim = (-0.8, -0.05)
        if levels is None: levels=np.linspace(*clim, 5)
        if cm is None: cm = sns.color_palette("rocket", as_cmap=True)
        plt.contourf(self.x, self.y, np.reshape(field, [self.nx, self.ny], order='F').T,
                        cmap=cm, levels=levels, vmin=clim[0], vmax=clim[1], extend='both')

        if black_contours:
            plt.contour(self.x, self.y, np.reshape(field, [self.nx, self.ny], order='F').T,
                            colors='k', linestyles='-', linewidths=0.5, levels=levels, vmin=clim[0], vmax=clim[1], extend='both')
            

class CylinderWake(NekFlowConfig):
    def __init__(self, field_path, **kwargs):
        self.prefix = 'cyl'
        super().__init__(field_path, **kwargs)

    def interpolation_setup(self, xspan=(-3, 9), yspan=(-3, 3), nx=400, ny=200, **kwargs):
        self.nx, self.ny = nx, ny
        self.x = np.linspace(*xspan, nx)
        self.y = np.linspace(*yspan, ny)
        self.XX, self.YY = np.meshgrid(self.x, self.y)
        self.interp_mask = (np.sqrt(self.XX**2 + self.YY**2) < 0.5).flatten('C')  # Zero out locations inside cylinder
        
    def plot_field(self, field, clim=None, levels=None, black_contours=True, cm=None):
        if clim is None: clim = (-5, 5)
        if levels is None: levels=np.linspace(*clim, 10)
        if cm is None: cm = 'RdBu'
        plt.contourf(self.x, self.y, np.reshape(field, [self.nx, self.ny], order='F').T,
                     cmap=cm, levels=levels, vmin=clim[0], vmax=clim[1], extend='both')

        cyl = plt.Circle((0, 0), 0.5, edgecolor='k', facecolor='gray')
        plt.gcf().gca().add_artist(cyl)

        if black_contours:
            plt.contour(self.x, self.y, np.reshape(field, [self.nx, self.ny], order='F').T,
                        colors='k', linestyles='-', linewidths=0.5, levels=levels, vmin=clim[0], vmax=clim[1], extend='both')
            
