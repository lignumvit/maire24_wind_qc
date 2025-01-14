import sys
sys.path.append('..')
import qc_utils
import wind_utils
import numpy as np

nc_dict = qc_utils.open_nc('./')
flights = list(nc_dict.keys())
flight = flights[0]

# grab the coefficients used for AKRD and SSLIP
orig_akrd_coefs = nc_dict[flight].variables['AKRD'].getncattr('CalibrationCoefficients')
orig_sslip_coefs = nc_dict[flight].variables['SSLIP'].getncattr('CalibrationCoefficients')


def test_akrd_calc():
    # test_akrd_calc tests both the calc_mach_dry and calc_akrd functions by computing
    # AKRD and comparing it against that in the NetCDF flight files. The AKRD calculation
    # uses ADIFR, QCF, and PSF measurements, uncorrected for angle of attack and sideslip
    # (i.e. not QCFC and PSFC)

    # read data into dataframe
    data_df = qc_utils.read_nc(nc_dict[flight])

    # get reference akrd calculated by nimbus
    akrd = data_df['AKRD'].to_numpy()

    # compute akrd with code here
    ratio = data_df['ADIFR']/data_df['QCF']
    mach_dry = wind_utils.calc_mach_dry(data_df['QCF'], data_df['PSF'])

    x = np.array([ratio, mach_dry])
    akrd_test = wind_utils.calc_akrd(x, *orig_akrd_coefs, data_df['QCF'])

    #diff = akrd - akrd_test
    #print(np.nanmax(np.abs(diff)))
    # max diff was 2.2294e-6 deg

    # iterate through original and test akrd and ensure all are within 1e-5 deg
    for i in range(len(akrd)):
        diff = akrd[i] - akrd_test[i]
        if not np.isnan(akrd[i]):
            assert np.isclose(akrd[i], akrd_test[i], rtol=0, atol=1e-5)

def test_wind_calc():
    # read data into dataframe
    data_df = qc_utils.read_nc(nc_dict[flight])

    wind_calc.calc_winds(data_df, orig_akrd_coefs, orig_sslip_coefs)

if __name__ == "__main__":
    test_akrd_calc()
