from LisbonTPMStool.mesh_functions import mesh_from_array
from LisbonTPMStool.mesh_functions import STL_from_mesh
from LisbonTPMStool import TPMS
from LisbonTPMStool.Utilities import gradient
import numpy as np
import os

def to_stl(SG, name='TPMS'):
    vertices, faces, v_normals = mesh_from_array(im=SG.im, dimensions=SG.dimensions)

    #Plot it
    mesh = [vertices, faces, 'oldlace', 1.0] #inpute to pass to PyVista_TriMeshes_plot

    path = os.path.dirname(os.path.realpath(__file__))
    STL_from_mesh(vertices, faces, name=name, path=path)

def main():
    """Only for debugging"""
    SG = TPMS('gyroid', dimensions=10.0, voxel_size=0.02)
    SG.cell_size_config(2)

    #Setup cell sizes to use
    initial_c = 0.1
    final_c = 0.9

    c_grad = gradient(domain=SG.domain, initial_value=initial_c, final_value=final_c, f=lambda x, y, z: x)

    #Apply it
    SG.level_set(c=c_grad)

    #Plot it
    mask = lambda x, y, z, im: np.where((x**2 + y**2 + z**2 <= np.power(np.amax(x), 2)), im, False)
    
    SG.im = mask(SG.domain[0], SG.domain[1], SG.domain[2], SG.im)

    to_stl(SG)


if __name__=="__main__":
    main()