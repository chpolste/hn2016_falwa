import json
from datetime import date, datetime
from calendar import monthrange
import cdsapi


all_pressure_level = [
            '1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '125', '150', '175', '200',
            '225', '250', '300', '350', '400', '450', '500', '550', '600', '650', '700', '750', '775',
            '800', '825', '850', '875', '900', '925', '950', '975', '1000']
time_in_a_day = ["00:00", "06:00", "12:00", "18:00"]
grid_resolution = "1.0/1.0"


if __name__ == "__main__":
    cdsapi_client = cdsapi.Client()
    cdsapi_client.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'pressure_level': all_pressure_level,
            'variable': ["temperature", "u_component_of_wind", "v_component_of_wind"],
            'time': time_in_a_day,
            'grid': grid_resolution,
            'product_type': 'reanalysis',
            'date': "2017-12-01/2017-12-31",
            'format': 'netcdf'
        },
        'test_20171202to20171231.nc'
    )
