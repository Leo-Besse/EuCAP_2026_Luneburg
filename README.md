# EuCAP_2026_Luneburg

Paper: "Fourier Series-Based Optimization of Gradient Index Lens Permittivity Profiles Within Manufacturing Limits", presented at the 2026 20th European Conference on Antennas and Propagation (EuCAP). 

This repository contains the scripts used to synthesize the TPMS-based Luneburg lens reflector presented in the above paper.

## Requirements

This code is built on [LisbonTPMS](https://github.com/JorgeESantos/LisbonTPMS-tool), with a minor modification (see note below).

This library works for Python 3.10 up to 3.12, and is dependent Numpy, Scipy, scikit-image, PoreSpy and PyVista. See the [LisbonTPMS](https://github.com/JorgeESantos/LisbonTPMS-tool) repository for further details.

## Usage

Create a virtual environment and install LisbonTPMS by running the following in a terminal:

    python -m venv .venv
    source .venv/bin/activate   # on Windows: .venv\Scripts\activate
    pip install -r requirements.txt

After this, the TPMS Luneburg lens can be synthesized by running the **EuCAP_Luneburg.py** script. Two hemispheresas shown below

![Voxellized TPMS Luneburg Lens with sub-unit cell grading](LL_image.png)

Alternatively, the **LL\_optimized.3mf** file can be opened in PrusaSlicer to directly print the Luneburg lens.

## Note: on the true_gradient function

By default, the _gradient_ function in the LisbonTPMS library takes 4 inputs;

    gradient(domain, initial_value, final_value, gradient_function)

Here, _initial\_value_ and _final\_value_ are the upper and lower bounds for the half-thickness of the gyroid. The _gradient\_function_ is then scaled to these limits. However, this poses a problem, as the gradient function is applied to the whole cubic domain, including points which are outside of the final spherical lens.

To deal with this issue, we create a new _true\_gradient_ function, which only takes two inputs;

    true_gradient(domain, gradient_function)

This function directly applies the output of the gradient function as the half-thickness. The user should be careful that the output of this function is between 0.2-1.5 everywhere in the intended domain.