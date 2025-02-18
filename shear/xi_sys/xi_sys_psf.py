from cosmosis.datablock import names, option_section
from cosmosis.datablock.cosmosis_py import lib
from astropy.io import fits
import numpy as np

#This file should be added to your cosmosis_standard_library following the path shear/xi_sys/xi_sys_psf.py
def setup(options):
    filename = options.get_string(option_section, 'data_file')
    data = fits.open(filename)
    rho_stats_name = options.get_string(option_section, 'rho_stats_name')
    rho_stats = data[rho_stats_name].data

    return rho_stats

def execute(block, config):
    #### !!! Works only in the non-tomographic case !!! ####
    #### !!! Will be updated to the tomographic case soon !!! ####
    rho_stats = config

    alpha, beta = block["psf_leakage_parameters", 'alpha'], block["psf_leakage_parameters", 'beta']

    xi_sys_p = (
        alpha**2*rho_stats["rho_0_p"]
        + beta**2*rho_stats["rho_1_p"]
        + 2*alpha*beta*rho_stats["rho_2_p"]
    )

    xi_sys_m = (
        alpha**2*rho_stats["rho_0_m"]
        + beta**2*rho_stats["rho_1_m"]
        + 2*alpha*beta*rho_stats["rho_2_m"]
    )

    theta = rho_stats["theta"]

    #CosmoSIS wants theta in radians
    theta = np.radians(theta / 60.)

    block['xi_sys', 'shear_xi_plus'] = xi_sys_p
    block['xi_sys', 'shear_xi_minus'] = xi_sys_m
    block['xi_sys', 'theta'] = theta

    return 0