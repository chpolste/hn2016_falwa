import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def plot_on_map(var, plot_title, xlon, ylat, eqidx=60, nhem=True, cmap='jet', figsize=(12, 8)):
    """
    Do a simple plot on a map in northern/southern hemisphere with the PlateCarree projection

    Args:
        var(np.ndarray): the 2D array that contains the variable to plot with size [nlat x nlon]
        plot_title(str): title of the plot
        xlon(np.array): a 1D array of longitude with size nlon
        ylat(np.array): a 1D array of latitude with size nlat of ascending latitude values
        eqidx(int): index of equator in the latitude array ylat
        nhem(bool): whether you are plotting the quantity over the northern hemisphere. Set it to False
            if you are plotting over the southern hemisphere.
        cmap(str): the name of the color map used. Default: 'jet'
        figsize(tuple of size 2): dimension of the figure size in inches
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=270))

    if nhem:
        half_var = var[eqidx:, :]
        cs = plt.contourf(xlon, ylat[eqidx:], half_var, 40, cmap=cmap,
                     transform=ccrs.PlateCarree())
    else:
        half_var = var[:eqidx, :]
        cs = plt.contourf(xlon, ylat[:eqidx], half_var, 40, cmap=cmap,
                     transform=ccrs.PlateCarree())

    plt.clim([0, half_var.max()])
    plt.colorbar(cs, orientation='horizontal', pad=0.05, aspect=50)
    plt.title(plot_title)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.xlocator = mticker.FixedLocator([-270, -180, -90, 0, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'k'}
    ax.coastlines()
    plt.show()