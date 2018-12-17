import numpy as np
from numpy import dtype
from math import pi
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import datetime as dt
#matplotlib inline
from hn2016_falwa.oopinterface import QGField
import hn2016_falwa.utilities as utilities
import datetime as dt


# --- Load the zonal wind and QGPV at 240hPa --- #
#datadir_input = '/data_storage/veredsil_data/JRA55_fill/fill_miss/'
#datadir_output = '/data_storage/veredsil_data/LWA_data/JRA55/'
datadir_input = '/data_storage/veredsil_data/TEM_code/NAT3DO3/fill_miss/'
datadir_output = '/data_storage/veredsil_data/LWA_data/WACCM3DO3/'
u_file = Dataset(datadir_input+'NATZMO3_U_1955_fill.nc', mode='r')
v_file = Dataset(datadir_input+'NATZMO3_V_1955_fill.nc', mode='r')
t_file = Dataset(datadir_input+'NATZMO3_T_1955_fill.nc', mode='r')

time_array = u_file.variables['time'][:] #u_file.variables['time'][:]
#time_units = u_file.variables['time'].units
#time_calendar = u_file.variables['time'].calendar
ntimes = time_array.shape[0]

print('Dimension of time: {}'.format(time_array.size))

xlon = u_file.variables['lon'][:]

# latitude has to be in ascending order
ylat = u_file.variables['lat'][:]
if np.diff(ylat)[0]<0:
    print('Flip ylat.')
    ylat = ylat[::-1]

# pressure level has to be in descending order (ascending height)
plev = u_file.variables['lev'][0:28]
if np.diff(plev)[0]>0:
    print('Flip plev.')
    plev = plev[::-1]

nlon = xlon.size
nlat = ylat.size
nlev = plev.size

clat = np.cos(np.deg2rad(ylat))     # cosine latitude
p0 = 1000.                          # surface pressure [hPa]
height = np.arange(0,48001,1000)    # pseudoheight [m]
dz = 1000.                          # differential height element
dphi = np.diff(ylat)[0]*pi/180.     # differential latitudinal element
dlambda = np.diff(xlon)[0]*pi/180.  # differential latitudinal element
hh = 7000.                          # scale height
cp = 1004.                          # heat capacity of dry air
rr = 287.                           # gas constant
omega = 7.29e-5                     # rotation rate of the earth
aa = 6.378e+6                       # earth radius
kmax = 49                           # number of grid points for vertical extrapolation (dimension of height)
prefactor = 6500.                   # integrated sum of density from ground to aloft
npart = nlat                        # number of partitions to construct the equivalent latitude grids
maxits = 100000                     # maximum number of iteration in the SOR solver to solve for reference state
tol = 1.e-5                         # tolerance that define convergence of solution
rjac = 0.95                         # spectral radius of the Jacobi iteration in the SOR solver.
jd = nlat//2+1                      # (one plus) index of latitude grid point with value 0 deg


year_start = 1955
year_end = 1955

for iyear in range(year_start,year_end+1):
    tstamp = [dt.datetime(iyear,1,1,0,0) + dt.timedelta(seconds=6*3600) * tt for tt in range(ntimes)]
    plev_selected = 10 # selected pressure level to display
    tstep_selected = 0 # selected time step to display

    u_file = Dataset(datadir_input+'NATZMO3_U_'+str(iyear)+'_fill.nc', mode='r')
    v_file = Dataset(datadir_input+'NATZMO3_V_'+str(iyear)+'_fill.nc', mode='r')
    t_file = Dataset(datadir_input+'NATZMO3_T_'+str(iyear)+'_fill.nc', mode='r')
    time_array = u_file.variables['time'][:] #u_file.variables['time'][:]

    #ttime = file_in_U.variables['time'][:]
    ntime = time_array.size
    Start_date = dt.datetime(iyear, 1, 1, 0, 0)
    delta_t = dt.timedelta(hours=24)
    Datestamp = [Start_date + delta_t * tt for tt in range(ntime)]


#qref, uref, ptref = test_object.compute_reference_states()
                                    # This is to be input to fortran code. The index convention is different.
    # === Outputing files ===
    output_file = Dataset(datadir_output+'WACCM3DO3_LWAREF_'+str(iyear)+'.nc',"w") #'2005-01-23_to_2005-01-30_output.nc'
    #output_file = Dataset(output_fname, 'w')
    output_file.createDimension('levelist',kmax)
    output_file.createDimension('latitude',nlat)
    output_file.createDimension('longitude',nlon)
    output_file.createDimension('time',ntime)
    plevs = output_file.createVariable('levelist',dtype('float32').char,('levelist',)) # Define the coordinate variables
    lats = output_file.createVariable('latitude',dtype('float32').char,('latitude',)) # Define the coordinate variables
    lons = output_file.createVariable('longitude',dtype('float32').char,('longitude',))
    times = output_file.createVariable('time',dtype('int').char,('time',))
    plevs.units = 'hPa'
    lats.units = 'degrees_north'
    lons.units = 'degrees_east'
    #times.units = time_units
    #times.calendar = time_calendar
    plevs[:] = p0 * np.exp(-height/hh)
    lats[:]  = ylat
    lons[:]  = xlon
    times[:] = time_array
    qgpv = output_file.createVariable('qgpv',dtype('float32').char,('time','levelist','latitude','longitude'))
    qgpv.units = '1/s'


    interpolated_u = output_file.createVariable('interpolated_u',dtype('float32').char,('time','levelist','latitude','longitude'))
    interpolated_u.units = 'm/s'
    interpolated_v = output_file.createVariable('interpolated_v',dtype('float32').char,('time','levelist','latitude','longitude'))
    interpolated_v.units = 'm/s'
    interpolated_theta = output_file.createVariable('interpolated_theta',dtype('float32').char,('time','levelist','latitude','longitude'))
    interpolated_theta.units = 'K'
    qref = output_file.createVariable('qref',dtype('float32').char,('time','levelist','latitude'))
    qref.units = '1/s'
    uref = output_file.createVariable('uref',dtype('float32').char,('time','levelist','latitude'))
    uref.units = 'm/s'
    ptref = output_file.createVariable('ptref',dtype('float32').char,('time','levelist','latitude'))
    ptref.units = 'K'
    lwa = output_file.createVariable('lwa',dtype('float32').char,('time','levelist','latitude','longitude'))
    lwa.units = 'm/s'
    adv_flux_f1 = output_file.createVariable('Zonal advective flux F1',dtype('float32').char,('time','latitude','longitude'))
    adv_flux_f1.units = 'm**2/s**2'
    adv_flux_f2 = output_file.createVariable('Zonal advective flux F2',dtype('float32').char,('time','latitude','longitude'))
    adv_flux_f2.units = 'm**2/s**2'
    adv_flux_f3 = output_file.createVariable('Zonal advective flux F3',dtype('float32').char,('time','latitude','longitude'))
    adv_flux_f3.units = 'm**2/s**2'
    adv_flux_conv = output_file.createVariable('Zonal advective flux Convergence -Div(F1+F2+F3)',dtype('float32').char,('time','latitude','longitude'))
    adv_flux_conv.units = 'm/s**2'

    divergence_eddy_momentum_flux = output_file.createVariable('Eddy Momentum Flux Divergence',dtype('float32').char,('time','latitude','longitude'))
    divergence_eddy_momentum_flux.units = 'm/s**2'
    meridional_heat_flux = output_file.createVariable('Low-level Meridional Heat Flux',dtype('float32').char,('time','latitude','longitude'))
    meridional_heat_flux.units = 'm/s**2'
    lwa_baro = output_file.createVariable('lwa_baro',dtype('float32').char,('time','latitude','longitude'))
    lwa_baro.units = 'm/s'
    u_baro = output_file.createVariable('u_baro',dtype('float32').char,('time','latitude','longitude'))
    u_baro.units = 'm/s'


    for tstep in range(ntimes):
        uu = u_file.variables['var131'][tstep, 0:28, :, :]
        vv = v_file.variables['var132'][tstep, 0:28, :, :]
        tt = t_file.variables['var130'][tstep, 0:28, :, :]

        #U = file_in_U.variables['var131'][it, 0:28, :, :]
        #V = file_in_V.variables['var132'][it, 0:28, :, :]
        #T = file_in_T.variables['var130'][it, 0:28, :, :]

        qgfield_object = QGField(xlon, ylat, plev, uu, vv, tt)

        qgpv[tstep, :, :, :], interpolated_u[tstep, :, :, :], interpolated_v[tstep, :, :, :], \
        interpolated_theta[tstep, :, :, :], static_stability = qgfield_object.interpolate_fields()

        qref[tstep, :, (nlat // 2):], uref[tstep, :, (nlat // 2):], ptref[tstep, :, (nlat // 2):] = \
            qgfield_object.compute_reference_states(northern_hemisphere_results_only=True)

        f1, f2, f3, f4, div_tmp, mhf_tmp, lwa_baro_tmp, utmp, lwa_tmp = qgfield_object.compute_lwa_and_barotropic_fluxes()


        adv_flux_f1[tstep, (nlat // 2):, :] = f1
        adv_flux_f2[tstep, (nlat // 2):, : ] = f2
        adv_flux_f3[tstep, (nlat // 2):, :] = f3
        adv_flux_conv[tstep, (nlat // 2):, :] = f4
        divergence_eddy_momentum_flux[tstep, (nlat // 2):, :] = div_tmp
        meridional_heat_flux[tstep, (nlat // 2):, :] = mhf_tmp
        lwa_baro[tstep, (nlat // 2):, :] = lwa_baro_tmp
        u_baro[tstep, (nlat // 2):, :] = utmp
        lwa[tstep, :, (nlat // 2):, :] = lwa_tmp
         #   = qgfield_object.compute_lwa_and_barotropic_fluxes()

output_file.close()