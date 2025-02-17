from cosmosis.datablock import names, option_section
from cosmosis.datablock.cosmosis_py import lib
from astropy.io import fits
import numpy as np
import astropy.units

DEFAULT_INPUT_SECTION = "shear_cl"
DEFAULT_OUTPUT_NAME = "bin_{}_{}"

def _build_rho_matrix(rho_stats):
    """
    Build the rho matrix from the rho_statistic data. This matrix is multiplied by the psf_leakage parameters to get the systematic contribution to the shear correlation functions.

    Parameters
    ----------
    rho_stats : astropy.table.Table
        Table containing the rho statistics
    
    Returns
    -------
    np.array
        A matrix containing the rho-statistics that allows to get the tau-statistics when multiplied by the psf_leakage parameters.
    """
    n_thetas = len(rho_stats['theta'])
    rho_matrix = np.empty((2 * n_thetas, 2))  # Preallocate the array for efficiency

    rho_matrix[:n_thetas, 0] = rho_stats["rho_0_p"]
    rho_matrix[:n_thetas, 1] = rho_stats["rho_2_p"]
    rho_matrix[n_thetas:, 0] = rho_stats["rho_2_p"]
    rho_matrix[n_thetas:, 1] = rho_stats["rho_1_p"]
    
    return rho_matrix

def _get_tau_theory(rho_matrix, alpha, beta):
    """
    Get the tau statistics given the rho matrix and the psf leakage parameters.

    Parameters
    ----------
    rho_matrix : np.array
        A matrix containing the rho-statistics that allows to get the tau-statistics when multiplied by the psf_leakage parameters.
    alpha : float
        The alpha psf leakage parameter.
    beta : float
        The beta psf leakage parameter.
    
    Returns
    -------
    np.array
        The tau statistics.
    """
    tau = np.dot(rho_matrix, np.array([alpha, beta]))

    return tau

#This file should be added to your cosmosis_standard_library following the path shear/xi_sys/xi_sys_psf.py
def setup(options):
    filename = options.get_string(option_section, 'data_file')
    rho_stats = fits.open(filename)['RHO_STATS'].data

    input_section = options.get_string(option_section, 'input_section', DEFAULT_INPUT_SECTION)

    theta = rho_stats['theta']
    rho_matrix = _build_rho_matrix(rho_stats)

    return theta, rho_matrix, input_section

def execute(block, config):

    theta, rho_matrix, input_section = config

    alpha, beta = block["psf_leakage_parameters", 'alpha'], block["psf_leakage_parameters", 'beta']
    #!!!Does not allow for tomography yet!!!
    tau_theory = _get_tau_theory(rho_matrix, alpha, beta)

    tau_0_p = tau_theory[:len(tau_theory) // 2]
    tau_2_p = tau_theory[len(tau_theory) // 2:]

    #Choose the bin value to go up to for the galaxy sample.
    if block.has_value(input_section, "nbin_a"):
        nbin_a = block[input_section, "nbin_a"]
    else:
        nbin_a = block[input_section, "nbin"]
    
    #The galaxies are always correlated with the same star catalog.
    #The star catalog does not have a redshift distribution.
    nbin_b = 1

    #Get metadata
    save_name = block.get_string(input_section, "save_name", default="")
    sample_a = block.get_string(input_section, "sample_a", default="")
    sample_b = block.get_string(input_section, "sample_b", default="")
    is_auto = block.get_bool(input_section, "is_auto", default=False)

    #Loop through bin pairs to write the output name.
    for i in range(nbin_a):
        for j in range(nbin_b):
            b1 = i+1
            b2 = j+1

            output_name = DEFAULT_OUTPUT_NAME.format(b1, b2) #Can be changed to allow other output names.

            #Save results back to cosmosis
            block["tau_0_plus", output_name] = tau_0_p
            block["tau_2_plus", output_name] = tau_2_p

    #Save metadata
    block["tau_0_plus", "nbin_a"] = nbin_a
    block["tau_0_plus", "nbin_b"] = nbin_b
    block["tau_0_plus", "sample_a"] = sample_a
    block["tau_0_plus", "sample_b"] = sample_b
    #CosmoSIS wants theta in radians
    arcmin_unit = astropy.units.arcmin
    theta = theta * arcmin_unit
    theta = theta.to(astropy.units.rad).value
    block["tau_0_plus", "theta"] = theta
    block["tau_0_plus", "sep_name"] = "theta"
    block["tau_0_plus", "save_name"] = save_name
    block["tau_0_plus", "is_auto"] = is_auto

    block["tau_2_plus", "nbin_a"] = nbin_a
    block["tau_2_plus", "nbin_b"] = nbin_b
    block["tau_2_plus", "sample_a"] = sample_a
    block["tau_2_plus", "sample_b"] = sample_b
    block["tau_2_plus", "theta"] = theta
    block["tau_2_plus", "sep_name"] = "theta"
    block["tau_2_plus", "save_name"] = save_name
    block["tau_2_plus", "is_auto"] = is_auto

    return 0