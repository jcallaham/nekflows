"""
Custom vector handles for working with Nek5000 data in modred
"""

import pymech.neksuite as nek
import modred as mr
import numpy as np
from copy import deepcopy

class NekVector(mr.Vector):
    def __init__(self, data_array):
        """
        data_array should be a numpy array here
            (or something else that support addition and scalar multiplication)
        """
        self.data_array = data_array

    def __add__(self, other):
        """Return a new object that is the sum of self and other"""
        sum_vec = deepcopy(self)
        sum_vec.data_array = self.data_array + other.data_array
        return sum_vec

    def __mul__(self, scalar):
        """Return a new object that is ``self * scalar`` """
        mult_vec = deepcopy(self)
        mult_vec.data_array = mult_vec.data_array * scalar
        return mult_vec


class NekHandle(mr.VecHandle):
    def __init__(self,
                 vec_path,
                 base_handle=None,
                 scale=None,
                 save_template=None):
        """
        vec_path - where to load the vector from
        template - optional path to load metadata from somewhere besides the vec_path on saving
        """
        mr.VecHandle.__init__(self, base_handle, scale)
        if save_template is None:  save_template = vec_path  # Default to using own path as template
        self.vec_path = vec_path
        self.save_template = save_template
            
        
    def _size(self, field):
        nel = len(field.elem) # Number of spectral elements
        nGLL = field.elem[0].pos.shape[-1] # Order of the spectral mesh
        n = nel*nGLL**2       # Total number of interpolation points
        return nel, nGLL, n
    

    def _get(self):
        # Load the field and read the number of spectral elements, order, etc
        field = nek.readnek(self.vec_path)
        nel, nGLL, n = self._size(field)

        u = np.array([field.elem[i].vel[0, 0, j, k]
                   for i in range(nel) for j in range(nGLL) for k in range(nGLL)])
        v = np.array([field.elem[i].vel[1, 0, j, k]
                   for i in range(nel) for j in range(nGLL) for k in range(nGLL)])
        p = np.array([field.elem[i].pres[0, 0, j, k]
                   for i in range(nel) for j in range(nGLL) for k in range(nGLL)])
        
        vec = np.concatenate((u, v, p))
        return NekVector(vec)
    

    def _put(self, vec):
        """
        Modify exadata file to have the prescribed velocity 'vec'
                and write to Nek format
        """
        # Load the field and read the number of spectral elements, order, etc
        field = nek.readnek(self.save_template)
        nel, nGLL, n = self._size(field)
            
        u = vec.data_array[:n]
        v = vec.data_array[n:2*n]
        p = vec.data_array[2*n:]
        vec_idx = 0
        for i in range(nel):
            for j in range(nGLL):
                for k in range(nGLL):
                    field.elem[i].vel[0, 0, j, k] = u[vec_idx]
                    field.elem[i].vel[1, 0, j, k] = v[vec_idx]
                    field.elem[i].pres[0, 0, j, k] = p[vec_idx]
                    vec_idx += 1
        nek.writenek(self.vec_path, field)
       
    def get_size(self):
        field = nek.readnek(self.save_template)
        nel, nGLL, n = self._size(field)
        return n
    

class ComplexNekHandle(mr.VecHandle):
    def __init__(self, real_path, imag_path, base_handle=None, scale=None, save_template=None):
        """
        NekHandle for cases with complex data (e.g. DMD modes)
        Can load/save real/imag parts separately
        
        vec_path - where to load the vector from
        template - optional path to load metadata from somewhere besides the vec_path on saving
        """
        mr.VecHandle.__init__(self, base_handle, scale)
        if save_template is None:  save_template = real_path  # Default to using own path as template
        self.real_path = real_path
        self.imag_path = imag_path
        self.save_template = save_template
        
    def _size(self, field):
        nel = len(field.elem) # Number of spectral elements
        nGLL = field.elem[0].pos.shape[-1] # Order of the spectral mesh
        n = nel*nGLL**2       # Total number of interpolation points
        return nel, nGLL, n

    def _get_vector(self, vec_path):
        # Load the field and read the number of spectral elements, order, etc
        field = nek.readnek(vec_path)
        nel, nGLL, n = self._size(field)

        u = np.array([field.elem[i].vel[0, 0, j, k]
                   for i in range(nel) for j in range(nGLL) for k in range(nGLL)])
        v = np.array([field.elem[i].vel[1, 0, j, k]
                   for i in range(nel) for j in range(nGLL) for k in range(nGLL)])
        p = np.array([field.elem[i].pres[0, 0, j, k]
                   for i in range(nel) for j in range(nGLL) for k in range(nGLL)])
        
        vec = np.concatenate((u, v, p))
        return NekVector(vec)
        
    def _get(self):
        real_part = self._get_vector(self.real_path)
        imag_part = self._get_vector(self.imag_path)
        return real_part + 1j*imag_part
        
    def _put_vector(self, vec_path, vec):
        """
        Modify exadata file to have the prescribed velocity 'vec'
                and write to Nek format
        """
        # Load the field and read the number of spectral elements, order, etc
        field = nek.readnek(vec_path)
        nel, nGLL, n = self._size(field)
            
        u = vec.data_array[:n]
        v = vec.data_array[n:2*n]
        p = vec.data_array[2*n:]
        vec_idx = 0
        for i in range(nel):
            for j in range(nGLL):
                for k in range(nGLL):
                    field.elem[i].vel[0, 0, j, k] = u[vec_idx]
                    field.elem[i].vel[1, 0, j, k] = v[vec_idx]
                    field.elem[i].pres[0, 0, j, k] = p[vec_idx]
                    vec_idx += 1
        nek.writenek(vec_path, field)
        
    def _put(self, vec):
        real_part = NekVector(np.real(vec.data_array))
        imag_part = NekVector(np.imag(vec.data_array))
        self._put_vector(self.real_path, real_part)
        self._put_vector(self.imag_path, imag_part)
       
    def get_size(self):
        field = nek.readnek(self.save_template)
        nel, nGLL, n = self._size(field)
        return n
    
        
def mean(handles):
    from functools import reduce
    do_sum = lambda x1, x2: x1 + x2
    
    data = [h.get() for h in handles]
    return reduce(do_sum, data) * (1/len(data))


def project(mode_handles, data_handles, inner_product):
    """
    Computes projection coefficients for fields not necessarily the same as those used for POD
    """
    vec_space = mr.VectorSpaceHandles(inner_product)
    coeffs = vec_space.compute_inner_product_array(mode_handles, data_handles)
    return coeffs
