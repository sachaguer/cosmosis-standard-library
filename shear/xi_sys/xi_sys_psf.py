from cosmosis.datablock import names, option_section
from astropy.io import fits

def setup(options):
    filename = options.get_string('data_file')
    data = fits.open(filename)
    return data, filename

def execute(block, config):

    data, filename = config
    cosmo = names.cosmological_parameters
    xi_sys = names.add_xi_sys

    alpha, beta, eta = block[cosmo, 'alpha'], block[cosmo, 'beta'], block[cosmo, 'eta']

    xi_p, xi_m = data[block[xi_sys, 'data_sets'][0]], data[block[xi_sys, 'data_sets'][1]]

    rho_stats = data[block[xi_sys, 'rho_stats_name']].data

    xi_sys_p = (
        alpha**2*rho_stats["rho_0_p"]
        + beta**2*rho_stats["rho_1_p"]
        + eta**2*rho_stats["rho_3_p"]
        + 2*alpha*beta*rho_stats["rho_2_p"]
        + 2*beta*eta*rho_stats["rho_5_p"]
        + 2*alpha*eta*rho_stats["rho_4_p"]
    )

    xi_sys_m = (
        alpha**2*rho_stats["rho_0_m"]
        + beta**2*rho_stats["rho_1_m"]
        + eta**2*rho_stats["rho_3_m"]
        + 2*alpha*beta*rho_stats["rho_2_m"]
        + 2*beta*eta*rho_stats["rho_5_m"]
        + 2*alpha*eta*rho_stats["rho_4_m"]
    )

    xi_p.data["VALUE"] = xi_p.data["RAW_VALUE"]+xi_sys_p #SG: will add xi_sys for each sample. We need to keep in memory the value of the xi_sys
    xi_m.data["VALUE"] = xi_m.data["RAW_VALUE"]+xi_sys_m

    data.writeto(filename, overwrite=True)

    return 0