from dipole_lib import importDipoleData, monoExp
from fortran_dip import calc_corr
import numpy as np
from math import ceil
from scipy.constants import k, epsilon_0
import matplotlib.pyplot as plt
import scipy.optimize
import math
import os


def processAlltogether():
    # Define Constants
    conversion_factor = 3.33564 * 1e-30  # convert Debye to Cm
    T = 300  # thermodynamic temperature in Kelvin
    psPerFrame = 0.5  # ps between frames in dip file
    max_lag_ps = 500  # maximum evaluated point for correlation function in ps
    num_of_traj = 40  # number of trajectories
    amk_name = "cys"
    conc = "150"
    dip_folder = "/home/jirka/WORK/UFE/A1704/MD/" + amk_name + "/" + conc + "/DIP/"
    spec_folder = "/home/jirka/WORK/UFE/A1704/MD/" + amk_name + "/" + conc + "/spec/"
    plot_folder = "/home/jirka/WORK/UFE/A1704/PLOTS/"
    f_name_list = [amk_name + '_p_', amk_name + '_w_']
    #################
    # Create folders for output
    if not os.path.exists(spec_folder):
        os.mkdir(spec_folder)
    if not os.path.exists(plot_folder):
        os.mkdir(plot_folder)
    #################
    max_lag_frames = ceil(max_lag_ps / psPerFrame)
    corr_time = np.linspace(0, max_lag_ps, max_lag_frames + 1)
    nof_components = len(f_name_list)
    meanCorr = np.zeros((nof_components, nof_components, len(corr_time)))
    for i in range(num_of_traj):
        i += 1
        dip_part = []
        # import peptide and water part of dipole data
        for file_name in f_name_list:
            dipFileName = dip_folder + file_name + "%i.dip" % i
            # import dipole file
            dip_data, V = importDipoleData(dipFileName)
            dip_array = np.array(dip_data[['dip_x', 'dip_y', 'dip_z']])
            dip_array = np.transpose(dip_array)
            dip_part.append(dip_array)
        # calculate auto-correlation and cross-correlation functions
        for m in range(nof_components):
            for n in range(nof_components):
                corr_fun = calc_corr(dip_part[m], dip_part[n], max_lag_frames)
                # Convert to proper units
                constant = conversion_factor ** 2 / (V * 3 * k * T * epsilon_0)
                corr_fun = corr_fun * constant
                meanCorr[m, n] += corr_fun

    meanCorr /= num_of_traj
    freq = np.logspace(8, 11.5, 1048)
    omega = 2 * np.pi * freq
    total_susci = 0
    corr_fig, axs_corr = plt.subplots(nof_components**2)
    suscii_fig, axs_suscii = plt.subplots(2)
    axs_corr[0].set_title('Dipole moment correlation functions')
    axs_corr[nof_components**2 - 1].set_xlabel('Time [ps]')
    leg_text = []
    for m in range(nof_components):
        for n in range(nof_components):
            # Normalize correlation function
            A = meanCorr[m, n, 0]
            norm_corr = meanCorr[m, n] / A
            # Perform fit
            p0 = (10)  # start with values near those we expect
            params, cv = scipy.optimize.curve_fit(monoExp, corr_time, norm_corr, p0)
            tau_ps = params
            print('Static susceptibility of sample equals to %.2f' % A)
            print('Relaxation time of Debye process is for process %i, %i equal to %.2f ps (%.2f GHz)' % (m, n, tau_ps, 1e3/(tau_ps*2*math.pi)))
            # plot individual correlation functions
            idx = (m * nof_components + n)
            axs_corr[idx].plot(corr_time, meanCorr[m, n], 'r', linewidth=2.0)
            axs_corr[idx].plot(corr_time, A*monoExp(corr_time, tau_ps))
            axs_corr[idx].legend(['corr.fun component %i %i' % (m, n),'monoexp fit'])
            # Calculate the FT of monoexp decay function by analytical formula
            tau = tau_ps*1e-12
            susc = np.zeros(len(freq), dtype='complex_')
            for i in range(len(freq)):
                susc[i] = A - omega[i] * (A * tau**2 * omega[i] / (1 + tau**2 * omega[i]**2))
                susc[i] = susc[i] + 1j * (omega[i] * A * tau / (1 + tau**2 * omega[i]**2))
            total_susci += susc
            # Plot resulting susceptibility
            color_code = np.random.rand(3,)
            np.savetxt(spec_folder + "component_%i_%i.csv" % (m,n), np.vstack((freq, susc)).T,
                       header='Freq  Intensity', delimiter="\t")
            axs_suscii[0].plot(freq, susc.real, c=color_code)
            axs_suscii[1].plot(freq, susc.imag, c=color_code)
            axs_suscii[1].set_xlabel('Frequency [Hz]')
            axs_suscii[0].set_ylabel('Real susceptibility')
            axs_suscii[1].set_ylabel('Imag. susceptibility')
            axs_suscii[0].set_xscale('log')
            axs_suscii[1].set_xscale('log')
            leg_text.append("%s vs. %s (%.2f ps)" % (f_name_list[m], f_name_list[n], tau_ps))

    axs_suscii[0].plot(freq, total_susci.real)
    axs_suscii[1].plot(freq, total_susci.imag)
    leg_text.append("Total susceptibility")
    axs_suscii[1].legend(leg_text)
    plt.savefig(plot_folder + amk_name + conc + '.png', dpi=199)
    np.savetxt(spec_folder + "total.csv", np.vstack((freq, total_susci)).T,
               header='Freq  Intensity', delimiter="\t")
    plt.show()


processAlltogether()

