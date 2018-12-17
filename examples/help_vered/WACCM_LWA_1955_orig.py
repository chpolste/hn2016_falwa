from hn2016_falwa.wrapper import qgpv_eqlat_lwa # Module for plotting local wave activity (LWA) plots and
                        # the corresponding equivalent-latitude profile
from math import pi
from hn2016_falwa.oopinterface import curl_2D
from hn2016_falwa.utilities import static_stability,compute_qgpv_givenvort
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mtpltlb
import datetime as dt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from scipy import interpolate
from my_utils.vered_utilities import compute_vertical_derivative_3d_rho
from hn2016_falwa.wrapper import qgpv_eqlat_lwa_ncforce
#%matplotlib inline

# --- Parameters --- #
Earth_radius = 6.378e+6 # Earth's radius
omega = 7.292462e-5 #2*pi/86400   # Earth angular rotation velocity [1/s]
H = 7000. # scale-height
cp = 1004.
rr = 287.
P0 = 1000.
rkappa = rr/cp
deltaH = 1000.
kmax = 49
kmin = 0
kmax_zhalf = kmax*2-1
latlim_sh = 48
latlim_nh = 48
DAY = 86400.0

datadir_input = '/data_storage/veredsil_data/TEM_code/NAT3DO3/fill_miss/'
datadir_output = '/data_storage/veredsil_data/LWA_data/WACCM3DO3/'
#datadir_input = '/carbon_data/veredsil/WACCM3DO3_fill/fill_miss/'
#datadir_output = '/data_storage/veredsil_data/LWA_data/WACCM3DO3/'
             # Corresponding pressure level of the height array

# --- Datestamp ---
Start_date = dt.datetime(1990, 1, 1, 0, 0)
delta_t = dt.timedelta(hours=24)
Datestamp = [Start_date + delta_t*tt for tt in range(4)]

# --- Read in longitude and latitude arrays --- #
file_in_tmp = Dataset(datadir_input+'NATZMO3_U_1955_fill.nc', mode='r')
xlon = file_in_tmp.variables['lon'][:]
ylat = file_in_tmp.variables['lat'][:]
ylat = ylat[::-1]
plev = file_in_tmp.variables['lev'][:]
# up to 28 is 0.1mb
plev = plev[0:28]
clat = np.abs(np.cos(ylat * np.pi / 180.))  # cosine latitude
clat[0] = 0.0
clat[-1] = 0.0
nlon = xlon.size
nlat = ylat.size
nlev = plev.size
slat = (1.-clat**2)**0.5
flat = 2.*omega*slat
# --- Parameters needed to use the module HN2015_LWA --- #
dphi = (ylat[2] - ylat[1]) * pi / 180.  # Equal spacing between latitude grid points, in radian
dlambda = (xlon[2] - xlon[1]) * pi / 180.
area = 2. * pi * Earth_radius ** 2 * (np.cos(ylat[:, np.newaxis] * pi / 180.) * dphi) / float(nlon) * np.ones((nlat, nlon))
area = np.abs(area)  # To make sure area element is always positive (given floating point errors).
area[0,:] = 0
area[-1,:]= 0
zlog = -H * np.log(plev / P0)
unih = np.array([i * deltaH for i in range(kmax)])  # Uniform height
unihhalf = np.array([i * deltaH / 2. for i in range(2 * kmax - 1)])  # Uniform height at every 500m
unih_p = P0 * np.exp(-unih / H)  # Corresponding pressure level of the height array
unihhalf_p = P0 * np.exp(-unihhalf / H)
dmu = Earth_radius*clat*dphi
f1 = 2 * omega * np.sin(ylat * pi / 180)

# calculate potential temperature
#theta = potential_temperature(plev, T, P0=1000, height_dim = 1)
year_start = 1955
year_end = 1955

for iyear in range(year_start,year_end+1):
    # --- Load U,V,T--- #
    file_in_T = Dataset(datadir_input+'NATZMO3_T_'+str(iyear)+'_fill.nc', mode='r')
    file_in_U = Dataset(datadir_input+'NATZMO3_U_'+str(iyear)+'_fill.nc', mode='r')
    file_in_V = Dataset(datadir_input+'NATZMO3_V_'+str(iyear)+'_fill.nc', mode='r')
    file_in_SWR = Dataset(datadir_input + 'NATZMO3_QRSTOT_' + str(iyear) + '_fill.nc', mode='r')
    file_in_LWR = Dataset(datadir_input + 'NATZMO3_QRLTOT_' + str(iyear) + '_fill.nc', mode='r')
    file_in_TOT = Dataset(datadir_input + 'NATZMO3_PTTEND_' + str(iyear) + '_fill.nc', mode='r')

    
    # --- Read in the absolute vorticNATZMO3_T_ity field from the netCDF file --- #

    # --- Read in time array --- #
    #time_units = file_in_U.variables['time'].units
    #time_calendar = file_in_U.variables['time'].calendar
    ttime = file_in_U.variables['time'][:]
    ntime = ttime.size
    # --- Datestamp ---
    Start_date = dt.datetime(iyear, 1, 1, 0, 0)
    delta_t = dt.timedelta(hours=24)
    Datestamp = [Start_date + delta_t * tt for tt in range(ntime)]

    stat_stab_nh = np.zeros([ntime, kmax_zhalf])
    stat_stab_sh = np.zeros([ntime, kmax_zhalf])
    #QGPV = np.zeros([ntime, kmax, nlat, nlon])
    #LWA = np.zeros([ntime, kmax, nlat, nlon])
    #Qref = np.zeros([ntime, kmax, nlat])
    U_unih = np.zeros([ntime, kmax, nlat, nlon])
    THall_unih = np.zeros([ntime, kmax, nlat, nlon])

    # create file for full fields
    dataset = Dataset(datadir_output+'WACCM3DO3_LWA_'+str(iyear)+'_fix.nc', "w")#, format="NETCDF4_CLASSIC")
    level = dataset.createDimension('level', kmax)
    lat = dataset.createDimension('lat', nlat)
    lon = dataset.createDimension('lon', nlon)
    time = dataset.createDimension('time', None)
    times = dataset.createVariable('time', np.float32, ('time',))
    levels = dataset.createVariable('level', np.float32, ('level',))
    latitudes = dataset.createVariable('latitude', np.float32, ('lat',))
    longitudes = dataset.createVariable('longitude', np.float32, ('lon',))
    plevels = dataset.createVariable('p_level', np.float32, ('level',))
    latitudes[:] = ylat
    longitudes[:] = xlon
    levels[:] = unih
    plevels[:] = unih_p
    #times.units = time_units
    #times.calendar = Datestamp
    times[:] = ttime
    #file_th = dataset.createVariable('TH', np.float32, ('time','level', 'lat', 'lon'))
    #file_th.units = 'K'  # set the units attribute.
    file_qgpv = dataset.createVariable('QGPV', np.float32, ('time','level', 'lat', 'lon'))
    file_qgpv.units = '1/s'  # set the units attribute.
    file_LWA = dataset.createVariable('LWA', np.float32, ('time', 'level', 'lat', 'lon'))
    file_LWA.units = 'm/s'  # set the units attribute.
    file_Qref = dataset.createVariable('Qref', np.float32, ('time', 'level', 'lat'))
    file_Qref.units = '1/s'  # set the units attribute.
    #file_dzdiv = dataset.createVariable('dzdivr', np.float32, ('level', 'lat', 'lon'))
    #file_zdiv = dataset.createVariable('zdiv', np.float32, ('level', 'lat', 'lon'))
    file_U = dataset.createVariable('U', np.float32, ('time', 'level', 'lat', 'lon'))
    file_U.units = 'm/s'  # set the units attribute.

    # Create file for ZM variables
    dataset_zm = Dataset(datadir_output + 'WACCM3DO3_LWA_' + str(iyear) + '_zm_fix.nc', "w")  # , format="NETCDF4_CLASSIC")
    level = dataset_zm.createDimension('level', kmax)
    lat = dataset_zm.createDimension('lat', nlat)
    time = dataset_zm.createDimension('time', None)
    times = dataset_zm.createVariable('time', np.float32, ('time',))
    levels = dataset_zm.createVariable('level', np.float32, ('level',))
    plevels = dataset_zm.createVariable('p_level', np.float32, ('level',))
    latitudes = dataset_zm.createVariable('latitude', np.float32, ('lat',))
    latitudes[:] = ylat
    levels[:] = unih
    plevels[:] = unih_p
    level_half = dataset_zm.createDimension('level_half', 2 * kmax - 1)
    plevel_half = dataset_zm.createDimension('plevel_half', 2 * kmax - 1)
    levels_half = dataset_zm.createVariable('level_half', np.float32, ('level_half',))
    plevels_half = dataset_zm.createVariable('plevel_half', np.float32, ('plevel_half',))
    levels_half[:] = unihhalf
    plevels_half[:] = unihhalf_p
    # times.units = time_units
    # times.calendar = Datestamp
    times[:] = ttime
    file_th_zm = dataset_zm.createVariable('TH', np.float32, ('time', 'level', 'lat'))
    file_th_zm.units = 'K'  # set the units attribute.
    file_qgpv_zm = dataset_zm.createVariable('QGPV', np.float32, ('time', 'level', 'lat'))
    file_qgpv_zm.units = '1/s'  # set the units attribute.
    file_LWA_zm = dataset_zm.createVariable('LWA', np.float32, ('time', 'level', 'lat'))
    file_LWA_zm.units = 'm/s'  # set the units attribute.
    file_Qref_zm = dataset_zm.createVariable('Qref', np.float32, ('time', 'level', 'lat'))
    file_Qref_zm.units = '1/s'  # set the units attribute.
    file_U_zm = dataset_zm.createVariable('U', np.float32, ('time', 'level', 'lat'))
    file_U_zm.units = 'm/s'  # set the units attribute.
    #file_eps_zm = dataset_zm.createVariable('eps', np.float32, ('time', 'level', 'lat'))
    #file_eps_zm.units = ''  # set the units attribute.
    stat_stab_nh = dataset_zm.createVariable('stat_Nhalf',np.float32, ('time', 'level_half'))
    stat_stab_nh.units = 'K/m'  # set the units attribute.
    stat_stab_sh = dataset_zm.createVariable('stat_Shalf',np.float32, ('time', 'level_half'))
    stat_stab_sh.units = 'K/m'  # set the units attribute.
    file_SIGMA_SWR = dataset_zm.createVariable('SIGMA_SWR', np.float32, ('time', 'level', 'lat'))
    file_SIGMA_SWR.units = 'm/s'  # set the units attribute.
    file_SIGMA_LWR = dataset_zm.createVariable('SIGMA_LWR', np.float32, ('time', 'level', 'lat'))
    file_SIGMA_LWR.units = 'm/s'  # set the units attribute.
    file_SIGMA_TOT = dataset_zm.createVariable('SIGMA_TOT', np.float32, ('time', 'level', 'lat'))
    file_SIGMA_TOT.units = 'm/s'  # set the units attribute.
    #file_LWAswrtmp = dataset_zm.createVariable('LWA_swr_tmp', np.float32, ('time', 'level', 'lat'))
    #file_LWAswrtmp.units = 'm/s'  # set the units attribute.
    #file_LWAlwrtmp = dataset_zm.createVariable('LWA_lwr_tmp', np.float32, ('time', 'level', 'lat'))
    #file_LWAlwrtmp.units = 'm/s'  # set the units attribute.

    # Create file for ZM variables
    dataset_pv = Dataset(datadir_output + 'WACCM3DO3_QGPV_' + str(iyear) + '_model.nc', "w")  # , format="NETCDF4_CLASSIC")
    level = dataset_pv.createDimension('level', nlev)
    lat = dataset_pv.createDimension('lat', nlat)
    lon = dataset_pv.createDimension('lon', nlon)
    time = dataset_pv.createDimension('time', None)
    times = dataset_pv.createVariable('time', np.float32, ('time',))
    levels = dataset_pv.createVariable('level', np.float32, ('level',))
    zlevels = dataset_pv.createVariable('z_level', np.float32, ('level',))
    latitudes = dataset_pv.createVariable('latitude', np.float32, ('lat',))
    longitudes = dataset_pv.createVariable('longitude', np.float32, ('lon',))
    longitudes[:] = xlon
    latitudes[:] = ylat
    levels[:] = plev
    zlevels[:] = zlog
    # times.units = time_units
    # times.calendar = Datestamp
    times[:] = ttime
    file_qgpv_file = dataset_pv.createVariable('QGPV', np.float32, ('time', 'level', 'lat','lon'))
    file_qgpv_file.units = '1/s'  # set the units attribute.
    file_dzdiv = dataset_pv.createVariable('dzdiv', np.float32, ('time', 'level', 'lat','lon'))
    file_dzdiv.units = '1/s'  # set the units attribute.


    for it in range(0, ntime): #range(0, 2):
        U = file_in_U.variables['var131'][it,0:28,:,:]
        V = file_in_V.variables['var132'][it,0:28,:,:]
        T = file_in_T.variables['var130'][it,0:28,:,:]
        Q_SWR = file_in_SWR.variables['qrstot'][it,0:28,:,:]  # * DAY
        Q_LWR = file_in_LWR.variables['qrltot'][it,0:28,:,:]  # * DAY
        Q_TOT = file_in_TOT.variables['qrltot'][it,0:28,:,:]  # * DAY
        VOR = np.zeros_like((U))

        # calc absolute vorticity from u and v
        for iz in range(0, nlev - 1):
            # calculate absolute vorticity
            abs_vor_tmp = curl_2D(U[iz, :, :], V[iz, :, :], clat, dlambda, dphi, planet_radius=6.378e+6)
            #abs_vor_tmp = abs_vor_tmp.reshape(1, 1, nlat, nlon)
            VOR[iz, :, :] = abs_vor_tmp + 2*omega*np.sin(ylat[:,np.newaxis]*pi/180.)

        #for it in range(0, ntime-1):
        #f_TH = interpolate.interp1d(zlog, theta[it,:,:,:], axis=0, fill_value='extrapolate')
        f_T = interpolate.interp1d(zlog, T[:, :, :], axis=0, fill_value='extrapolate')
        f_avor = interpolate.interp1d(zlog, VOR[:,:,:], axis=0, fill_value='extrapolate')
        f_U = interpolate.interp1d(zlog, U[:,:,:], axis=0, fill_value='extrapolate')

        f_SWR = interpolate.interp1d(zlog, Q_SWR[ :, :, :], axis=0, fill_value='extrapolate')
        f_LWR = interpolate.interp1d(zlog, Q_LWR[:, :, :], axis=0, fill_value='extrapolate')
        f_TOT = interpolate.interp1d(zlog, Q_TOT[ :, :, :], axis=0, fill_value='extrapolate')

        #TH_unih = f_TH(unih)
        unihp = 1000. * np.exp(-unih / 7000.)  # Corresponding pressure level of the height array
        TH_unih = f_T(unih) * (1000. / unihp[:, np.newaxis, np.newaxis]) ** rkappa
        TH_unihhalf = f_T(unihhalf) * (1000. / unihhalf_p[:, np.newaxis, np.newaxis]) ** rkappa   # THETA at half steps (500m)

        file_th_zm[it, :, :] = np.mean(TH_unih, axis=-1)
        AVOR_unih = f_avor(unih)
        #U_unih[it, :, :, :] = f_U(unih)
        file_U[it, :, :, :]  = f_U(unih)

        # calculate static stability
        t0_N, t0_S, stat_stab_nh[it, :], stat_stab_sh[it, :] = static_stability(unihhalf, area, TH_unihhalf, s_et=nlat / 2, n_et=nlat / 2)
        #stat_stab_nh[it, :] = stat_Nhalf_tmp
        #stat_stab_sh[it, :] = stat_Shalf_tmp

        stat_cN = stat_stab_nh[it, ::2]  # every odd element of stat_Nhalf
        stat_cS = stat_stab_sh[it, ::2]  # every odd element of stat_Shalf

        t0_cN = t0_N[::2]  # every odd element of stat_Nhalf
        t0_cS = t0_S[::2]  # every odd element of stat_Shalf

        stat_cNpmhalf = stat_stab_nh[it, 1::2]  # every even element of stat_Nhalf
        stat_cSpmhalf = stat_stab_sh[it, 1::2]  # every even element of stat_Shalf

        #stat_half = np.zeros((kmax - 1, nlat))
        #stat_half[:, :nlat/2] = stat_cSpmhalf[:, np.newaxis] * np.ones((kmax-1, nlat/2))
        #stat_half[:, (-nlat/2+1):] = stat_cNpmhalf[:, np.newaxis] * np.ones((kmax-1, nlat/2))

        stat_stab = np.zeros((kmax, nlat))
        stat_stab[:, :nlat / 2] = stat_cS[:, np.newaxis] * np.ones((kmax, nlat / 2))
        stat_stab[:, (-nlat / 2 ):] = stat_cN[:, np.newaxis] * np.ones((kmax, nlat / 2))

        kmax_half = kmax*2-1
        stat_stab_half = np.zeros((kmax_half, nlat))
        ststabsh_tmp = stat_stab_sh[it,:]
        stat_stab_half[:, :nlat / 2] = ststabsh_tmp[:, np.newaxis] * np.ones((kmax_half, nlat / 2))
        ststabsh_tmp = stat_stab_nh[it, :]
        stat_stab_half[:, (-nlat / 2 ):] = ststabsh_tmp[:, np.newaxis] * np.ones((kmax_half, nlat / 2))

        #file_stat_Nhalf[it, :] = stat_cNpmhalf
        #stat_stab_sh[it,:] = stat_cSpmhalf

        # Compute QGPV
        file_qgpv[it, :, :, :], dzdiv = compute_qgpv_givenvort(omega, nlat, nlon, kmax, unih, ylat, AVOR_unih,
                                                       TH_unih, t0_cN, t0_cS, stat_cN, stat_cS)



        # interpolate qgpv and stretching term back to original grid and save to a different file
        f_qgpv_zlog = interpolate.interp1d(unih, file_qgpv[it, :, :, :], axis=0, fill_value='extrapolate')
        f_dzdiv_zlog = interpolate.interp1d(unih, dzdiv, axis=0, fill_value='extrapolate')
        file_qgpv_file[it,:,:,:] = f_qgpv_zlog(zlog)
        file_dzdiv[it,:,:,:] = f_dzdiv_zlog(zlog)


        # calculate etha from SWR
        #QrTERM_tmp = f_SWR(unih)
        #fdata = QrTERM_tmp / stat_stab[:, :, np.newaxis]
        #ETHA_SWR = compute_vertical_derivative_3d_rho(kmax, unih, ylat, fdata, scale_height=7000.)
        QrTERM_tmp = f_SWR(unihhalf)
        fdata = QrTERM_tmp / stat_stab_half[:, :, np.newaxis]
        ETHA_SWR = compute_vertical_derivative_3d_rho(kmax*2-1, unihhalf, ylat, fdata, scale_height=7000.)
        ETHA_SWR = ETHA_SWR * f1[np.newaxis, :, np.newaxis]
        ETHA_SWR = ETHA_SWR[::2,:,:]

        # calculate etha from LWR
        #QrTERM_tmp = f_LWR(unih)
        #fdata = QrTERM_tmp / stat_stab[:, :, np.newaxis]
        #ETHA_LWR = compute_vertical_derivative_3d_rho(kmax, unih, ylat, fdata, scale_height=7000.)
        QrTERM_tmp = f_LWR(unihhalf)
        fdata = QrTERM_tmp / stat_stab_half[:, :, np.newaxis]
        ETHA_LWR = compute_vertical_derivative_3d_rho(kmax*2-1, unihhalf, ylat, fdata, scale_height=7000.)
        ETHA_LWR = ETHA_LWR * f1[np.newaxis, :, np.newaxis]
        ETHA_LWR = ETHA_LWR[::2,:,:]

        # calculate etha from LWR
        QrTERM_tmp = f_TOT(unihhalf)
        fdata = QrTERM_tmp / stat_stab_half[:, :, np.newaxis]
        ETHA_TOT = compute_vertical_derivative_3d_rho(kmax*2-1, unihhalf, ylat, fdata, scale_height=7000.)
        ETHA_TOT = ETHA_TOT * f1[np.newaxis, :, np.newaxis]
        ETHA_TOT = ETHA_TOT[::2, :, :]

        #file_qgpv[it, :, :, :] = QGPV[it, :, :, :]

        #QGPV[it, :, :, :] = QGPV_tmp
        for ik in range(kmax):
            #file_Qref[it,ik,:], file_LWA[it,ik,:,:] = qgpv_eqlat_lwa(ylat, file_qgpv[it, ik, :, :], area,
             #                                            Earth_radius*clat*dphi)

            file_Qref[it, ik, :], file_LWA[it,ik,:,:], bigsigma_result = qgpv_eqlat_lwa_ncforce(ylat, file_qgpv[it, ik, :, :], ETHA_SWR[ik, :, :],
                                                                           area, dmu, nlat_s=None,
                                                                           n_points=None, planet_radius=6.378e+6)
            file_SIGMA_SWR[it, ik, :] = np.mean(bigsigma_result, axis=-1)

            qref_tmp, lwa_result, bigsigma_result = qgpv_eqlat_lwa_ncforce(ylat, file_qgpv[it, ik, :, :], ETHA_LWR[ik, :, :],
                                                                           area, dmu, nlat_s=None,
                                                                           n_points=None, planet_radius=6.378e+6)
            file_SIGMA_LWR[it, ik, :] = np.mean(bigsigma_result, axis=-1)
            qref_tmp, lwa_result, bigsigma_result = qgpv_eqlat_lwa_ncforce(ylat, file_qgpv[it, ik, :, :],
                                                                           ETHA_TOT[ik, :, :],
                                                                           area, dmu, nlat_s=None,
                                                                           n_points=None, planet_radius=6.378e+6)
            file_SIGMA_TOT[it, ik, :] = np.mean(bigsigma_result, axis=-1)

        file_Qref_zm[it, :, :] = file_Qref[it,:,:]
        file_LWA_zm[it, :, :] = np.mean(file_LWA[it,:,:,:], axis=-1)
        file_U_zm[it, :, :] = np.mean(file_U[it, :, :, :], axis=-1)
        file_qgpv_zm[it, :, :] = np.mean(file_qgpv[it, :, :, :], axis=-1)
        file_th_zm[it, :, :] = np.mean(THall_unih, axis=-1)

       
    dataset_zm.close()
    dataset.close()
    dataset_pv.close()

