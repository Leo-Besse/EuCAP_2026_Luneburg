from LisbonTPMStool import TPMS
from LisbonTPMStool.Utilities import gradient
from export_stl import to_stl
import numpy as np

def true_gradient(domain, gradient_function):
    """A small workaround for the behaviour of the gradient function in LisbonTPMS (see note at bottom of this script)"""
    initial_val = np.min(gradient_function(domain[0], domain[1], domain[2]))
    final_val = np.max(gradient_function(domain[0], domain[1], domain[2]))
    return gradient(domain, initial_val, final_val, gradient_function)

def inverse_MG(target, PLA, air):
    """Returns volume fraction of air given a target permittivity, according to the Maxwell Garnett approximation for air inclusions in a PLA matrix"""
    A = (target-PLA) / (target+2*PLA)
    B = (air-PLA) / (air+2*PLA)
    return A/B

def t_SHG(rho):
    """Quadratic fit relating the half-thickness of the gyroid t with the volume fraction rho (determined empirically)"""
    t = (-0.492 + np.sqrt(0.492**2 - 4*0.109*(0.0281 - rho))) / (2*0.109)
    t[t<0.2] = 0.2 # necessary to speed up synthesis; voxels outside the final domain with t < 0.2 add a lot of computational effort
    return t

def luneburg(x,y,z):
    """Returns target permittivity of a luneburg lens in cartesian coordinates"""
    r = (x**2 + y**2 + z**2)
    R = np.max(x**2)
    permittivity = 2 - (r/R)
    return permittivity

def MG_luneburg():
    """Returns TPMS object for the truncated Luneburg Lens"""
    unit_cell_size = 6e-3
    limit_of_manufacture = 0.905 # Value determined using the graph in plot_profile, assuming max ~14% printable volume fraction
    d = (6e-2) / limit_of_manufacture # s.t. truncated LL has OD of 6 cm

    draft_resolution = unit_cell_size/30
    fine_resolution = unit_cell_size/50

    SG = TPMS('gyroid', dimensions=d, voxel_size=fine_resolution)
    SG.cell_size_config(unit_cell_size)

    eps_air = 1
    eps_PLA = 2.55

    eps = true_gradient(SG.domain, luneburg)
    delta = 1-inverse_MG(eps, eps_PLA, eps_air)
    c_grad = t_SHG(delta)

    #Apply it
    SG.level_set(c=c_grad)

    #Plot it
    mask = lambda x, y, z, im: np.where(((x**2 + y**2 + z**2 <= (limit_of_manufacture*np.max(x))**2)), im, False)
    SG.im = mask(SG.domain[0], SG.domain[1], SG.domain[2], SG.im)

    mask = lambda x, y, z, im: np.where((y<0), im, False)
    SG.im = mask(SG.domain[0], SG.domain[1], SG.domain[2], SG.im)

    return SG


def plot_profile(epsilon, d):
    # plot a given permittivity gradient in one axis
    # The tpms domain is a (n x n x n) array, where n is the voxel dimension
    # In this case, we plot through the middle slice;
    r = d/2
    midpoint = np.size(epsilon[1,1,:]) // 2
    print(midpoint)
    eps_points = epsilon[midpoint, midpoint, midpoint:]
    x = np.linspace(0, 1, np.size(eps_points))
    
    import matplotlib.pyplot as plt
    plt.plot(x, eps_points)
    plt.show()


def get_N_Fourier(X, Y, Z, modes, coeffs, R0):
    X = np.asarray(X)
    Y = np.asarray(Y)
    Z = np.asarray(Z)

    # Per-point radius
    RHO = np.sqrt(X*X + Y*Y + Z*Z)

    # Accumulator (complex) with same shape as RHO
    acc = np.zeros(RHO.shape, dtype=np.complex128)

    for i, mode in enumerate(modes):
        k_rho = np.pi * mode / (2 * R0)
        acc += coeffs[i] * np.exp(1j * k_rho * RHO)

    return acc.real


def main():
    R0 = 3e-2 # Radius = 3cm

    # Optimization output with shape [n_max (magnitude of each mode) (phase of each mode)]
    x = [1.410811757628144, 0.839157502295643, 0.877016158265710, 0.309729688889705, 0.004419096763672, 0.013365365479853, 0.011162268046916, 0.099235951240807, 5.950399583340132, 2.602953709777274, 2.463560807199144, 0.465627103439559, 4.231313755656105]

    n_max = x.pop(0)
    n_min = 1.1
    modes = np.arange(len(x)//2) + 1
    c_mag = np.array(x[:len(x)//2])
    c_phase = np.array(x[len(x)//2:])

    coeffs = c_mag * np.exp(1j*c_phase)


    # Determine the scaling factors (see equations 2 and 3)
    xs = np.linspace(-R0, R0, 80)
    ys = np.linspace(-R0, R0, 80)
    zs = np.linspace(-R0, R0, 80)

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")

    Rarr = np.stack((X, Y, Z), axis=0) # shape (3, 80, 80, 80)
    Mask = np.linalg.norm(Rarr, axis=0) <= R0

    Nraw = get_N_Fourier(*(comp[Mask] for comp in Rarr), modes, coeffs, R0)
    smin = np.min(Nraw) 
    smax = np.max(Nraw)
    sf = (n_max-n_min)/(smax-smin) 
    coeffs = sf * coeffs
    c0 = -sf*smin+n_min

    Fourier_grad_function = lambda x, y, z: (
    get_N_Fourier(*np.meshgrid(x*2*R0/(2*np.pi),
                               y*2*R0/(2*np.pi),
                               z*2*R0/(2*np.pi),
                               indexing="ij"),
                  modes, coeffs, R0) + c0)**2

    unit_cell_size = 6e-3 # 6mm unit cell size -> ~40GHz max frequency (see https://doi.org/10.1109/USNC-URSI52151.2023.10237634)

    # voxel size
    draft_resolution = unit_cell_size/30
    fine_resolution = unit_cell_size/50

    SG = TPMS('gyroid', dimensions=2*R0, voxel_size=fine_resolution)
    SG.cell_size_config(unit_cell_size)

    eps_air = 1
    eps_PLA = 2.55

    eps = true_gradient(SG.domain, Fourier_grad_function)
    delta = 1-inverse_MG(eps, eps_PLA, eps_air)
    c_grad = t_SHG(delta)
    SG.level_set(c=c_grad)

    #Plot it
    mask = lambda x, y, z, im: np.where(((x**2 + y**2 + z**2 <= (np.max(x))**2)), im, False)
    SG.im = mask(SG.domain[0], SG.domain[1], SG.domain[2], SG.im)

    # In FFF, it is difficult to print a complete sphere in one print.
    # The recommended procedure it to print two hemispheres & attach them together, which is what we do here
    
    # Uncomment one of the two lines below to synthesize either hemisphere:
    # mask = lambda x, y, z, im: np.where((y>0), im, False) # 'plus' hemisphere
    mask = lambda x, y, z, im: np.where((y<0), im, False) # 'minus' hemisphere

    SG.im = mask(SG.domain[0], SG.domain[1], SG.domain[2], SG.im)

    return SG


if __name__ == "__main__":
    
    SG = main() # Synthesize LL hemisphere with optimized profile
    # SG = MG_luneburg() # Synthesize LL with truncated profile

    SG.im_visualize(save_fig=False) # Visualize the TPMS structure (requires openGL)
    # to_stl(SG, name="EUCAP_Optimized_LL_minus") # Write to .stl file in current directory






"""
Note: On the true_gradient function

By default, the gradient function in LisbonTPMS takes 4 inputs;

    gradient(domain, initial_value, final_value, gradient_function)

Here, initial_value and final_value are the upper and lower bounds for the half-thickness of the gyroid. 
The gradient_function is then scaled to these limits. 

However, this poses a problem, as the gradient function is applied to the whole cubic domain, 
including points which are outside of the final spherical lens.

To deal with this issue, we define a new true_gradient function, which only takes two inputs;

    true_gradient(domain, gradient_function)

This function directly applies the output of the gradient function as the half-thickness. 
The user should be careful that the output of this function is between ~0.2 to 1.5 everywhere in the intended domain
(as those are the limits for the gyroid half-thickness. Well, at least, 1.5 is; the lower bound depends how thin you can make the walls.)

"""