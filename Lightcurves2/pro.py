from ciao_contrib.runtool import search_csc
from pyvo import DALFormatError, dal
from tqdm import tqdm  
from astropy import units
from astropy.coordinates import SkyCoord
from astropy.coordinates.name_resolve import NameResolveError
from ciao_contrib.runtool import dmcopy, dmextract, dmkeypar, dmlist, dmstat, new_pfiles_environment
import gzip 
import glob
from pathlib import Path
import os
from pandas import DataFrame
from astropy import io, table
import astropy.units as u
import numpy as np
from astropy.io import ascii, votable
import csv
import time 
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import matplotlib
from postage_stamp_plotter import CropBounds, plot_postagestamps
from astropy.time import Time
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from astropy.stats import bayesian_blocks
from astropy.timeseries import LombScargle
from scipy import stats

np.random.seed(0)

ENERGY_LEVELS = {
    "broad": "energy=150:7000",
    "ultrasoft": "energy=150:300",
    "soft": "energy=300:1200",
    "medium": "energy=1200:2000",
    "hard": "energy=2000:7000",
}

ITERATION_SIM = 10000

search_radius = 10 * units.arcmin #120
object_name = "01:32:36.83 +30:32:29.83" # in M33
object_is_pos = False
download_directory = './results/data_test' 
graphs_directory = './results/output_test'
dipflare_directory = './results/dipflare_test'
## NO LOGS CREATED HERE
bin_size = 500
significance_threshold = 3
min_avg_count_rate = 0.008
dip_threshold = 0.7
flare_threshold = 0.7
significance_of_dipflare_threshold = 0.2
start_with_downloading = True
already_considered_RA = './results/alreadyconsideredRA_test'
already_considered_DE = './results/alreadyconsideredDE_test'
create_just_big = False # only create csv files and graphs for things above the min_avg_count_rate
# graphs are created as .png here (but can also change to make them .svg)

# search_radius = 0.2 * units.arcmin #120
# object_name = "01:00:43.0 -72:11:33"
# object_is_pos = True
# download_directory = './results/data_goodgraph' 
# graphs_directory = './results/output_goodgraph'
# dipflare_directory = './results/dipflare_goodgraph'
# ## NO LOGS CREATED HERE
# bin_size = 500
# significance_threshold = 10
# min_avg_count_rate = 0.01
# dip_threshold = 0.7
# flare_threshold = 0.7
# significance_of_dipflare_threshold = 0.2
# start_with_downloading = True
# already_considered_RA = './results/alreadyconsideredRA_goodgraph' 
# already_considered_DE = './results/alreadyconsideredDE_goodgraph'
# create_just_big = False # only create csv files and graphs for things above the min_avg_count_rate
# # graphs are created as .png here (but can also change to make them .svg)


if start_with_downloading: 
    if not os.path.exists(download_directory):
        os.makedirs(download_directory)
    if not os.path.exists(already_considered_RA):
        os.makedirs(already_considered_RA)
    if not os.path.exists(already_considered_DE):
        os.makedirs(already_considered_DE)

    if not object_is_pos: 
        sky_coord = SkyCoord.from_name(object_name)
    else: 
        sky_coord = SkyCoord(object_name, unit=(u.hourangle, u.deg))

    csc_conesearch = "http://cda.cfa.harvard.edu/csc2scs/coneSearch"
    search_results = dal.conesearch(csc_conesearch, sky_coord, search_radius, verbosity=2)

    # https://cxc.cfa.harvard.edu/csc/columns/stack_alpha.html
    # significance = ratio of the flux measurement to its average error
    # average error = (high_bound_flux - low_bound_flux)/2

    print(len(search_results))

    significant_results = [
        result for result in search_results if result["significance"] >= significance_threshold
        ]
    significant_results_count = len(significant_results)
    print(f"Found {significant_results_count} sources meeting the significance threshold.")
    if significant_results_count == 0:
        exit()

    starttime=time.time()

    for source_count, source in enumerate(significant_results): #[::-1]
        right_ascension = source["ra"]
        declination = source["dec"]

        RApath = os.path.join(already_considered_RA, str(right_ascension) + '.txt')
        DEpath = os.path.join(already_considered_DE, str(declination) + '.txt')

        if os.path.exists(RApath) and os.path.exists(DEpath): 
            print('already exists')
            continue

        search_csc(
            pos=f"{right_ascension}, {declination}",
            radius="1.0",
            outfile="search-csc-outfile.tsv",
            radunit="arcsec",
            columns="",
            sensitivity="no",
            download="all",
            root=download_directory,
            bands="broad, wide",
            filetypes="regevt, reg",
            catalog="csc2",
            verbose="0",
            clobber="1",
        )
        if source_count%5==0:
            print(source_count)
            print(time.time()-starttime)
            starttime=time.time()

        file2write=open(RApath,'w')
        file2write.write('')
        file2write.close()
        file2write=open(DEpath,'w')
        file2write.write('')
        file2write.close()
        



## EXISTING FUNCTIONS IN LIGHTCURVES, MINOR MODIFICATION


def unzip_fits_files(observation_path: Path):
    """Files are downloaded in a GNU Zip format. This unzips all files in an observation
    directory into the FITS files we need."""
    gzip_files = observation_path.glob("*.gz")
    try:
        for gzip_file in gzip_files:
            with gzip.open(gzip_file, "rb") as gzipped_file:
                unzipped_data = gzipped_file.read()
                with open(gzip_file.with_suffix(""), "wb") as unzipped_file:
                    unzipped_file.write(unzipped_data)
                gzip_file.unlink()
    except gzip.BadGzipFile as error:
        raise RuntimeError(f"Could not unzip file {gzip_file.name}.") from error

def isolate_source_region(event_list: Path, source_region):
    """Restrict the event list to just the source region of the source we care about."""
    dmcopy(
        infile=f"{event_list}[sky=region({source_region})]",
        outfile=(outfile := f"{event_list.with_suffix('.src.fits')}"),
        clobber="yes",
    )
    return outfile


def adjust_bin_size(event_list, bin_size):
    """For ACIS, time resolution can be in the seconds in timed exposure mode, as compared to in
    the microseconds for HRC. Thus we must round the bin_size to the time resolution."""
    time_resolution = float(dmkeypar(infile=str(event_list), keyword="TIMEDEL", echo=True))
    return bin_size // time_resolution * time_resolution


def extract_lightcurves(event_list, bin_size):
    outfiles = []
    bin_size = adjust_bin_size(event_list, bin_size)
    for light_level, energy_range in ENERGY_LEVELS.items():
        dmextract(
            infile=f"{event_list}[{energy_range}][bin time=::" f"{bin_size}]",
            outfile=(outfile := f"{event_list}.{light_level}.lc"),
            opt="ltc1",
            clobber="yes",
        )
        outfiles.append(Path(outfile))
    return outfiles


def filter_lightcurve_columns(lightcurves):
    outfiles = []
    for lightcurve in lightcurves:
        dmlist(
            infile=f"{lightcurve}"
            f"[cols time,count_rate,count_rate_err,counts,exposure,area]",
            opt="data,clean",
            outfile=(outfile := f"{lightcurve}.ascii"),
        )
        outfiles.append(Path(outfile))
    return outfiles
    

def get_lightcurve_data(lightcurves: list[Path]):
    lightcurve_data: dict[str, DataFrame] = {
        energy_level: table.Table.read(lightcurve, format="ascii").to_pandas() 
        for energy_level, lightcurve in zip(ENERGY_LEVELS.keys(), lightcurves)
    }
    # Trim zero exposure points
    
    for energy_level, lightcurve_dataframe in lightcurve_data.items():
        lightcurve_data[energy_level] = lightcurve_dataframe[
            lightcurve_dataframe["EXPOSURE"] != 0
        ]
    return lightcurve_data



def get_images(event_list, region_event_list, sky_coords_image=None, detector_coords_image=None):
    """Gets the images from the event list, one in sky coordinates, the other in detector
    coordinates. This is useful to be able to track if a lightcurve dip corresponds with a
    source going over the edge of the detector. The image is cropped otherwise we would be
    dealing with hundreds of thousands of blank pixels. The cropping bounds are written to the
    FITS file for later use in plotting."""

    dmstat(infile=f"{region_event_list}[cols x,y]")
    sky_bounds = CropBounds.from_strings(*dmstat.out_min.split(","), *dmstat.out_max.split(","))
    sky_bounds.double()
    dmcopy(
        infile=f"{event_list}"
        f"[bin x={sky_bounds.x_min}:{sky_bounds.x_max}:0.5,"
        f"y={sky_bounds.y_min}:{sky_bounds.y_max}:0.5]",
        outfile=(sky_coords_image := f"{region_event_list}.skyimg.fits"),
    )
    with io.fits.open(sky_coords_image, mode="append") as hdu_list:
        hdu_list.append(sky_bounds.to_hdu())

    dmstat(infile=f"{region_event_list}[cols detx,dety]")
    detector_bounds = CropBounds.from_strings(
        *dmstat.out_min.split(","), *dmstat.out_max.split(",")
    )
    detector_bounds.add_padding(x_padding=5, y_padding=5)
    dmcopy(
        infile=f"{event_list}"
        f"[bin detx={detector_bounds.x_min}:{detector_bounds.x_max}:0.5,"
        f"dety={detector_bounds.y_min}:{detector_bounds.y_max}:0.5]",
        outfile=(detector_coords_image := f"{region_event_list}.detimg.fits"),
    )
    with io.fits.open(detector_coords_image, mode="append") as hdu_list:
        hdu_list.append(detector_bounds.to_hdu())
    
    return sky_coords_image, detector_coords_image


def get_lightcurve_counts(lightcurve_data):
    return int(lightcurve_data["broad"]["COUNTS"].sum())


def create_csv(lightcurve_data, file_relativeS, graphs_directory):
    """Create CSV with columns for each energy level."""
    combined_data = DataFrame({
        "time": lightcurve_data["broad"]["TIME"],
        "count_rate": lightcurve_data["broad"]["COUNT_RATE"],
        "counts": lightcurve_data["broad"]["COUNTS"],
        "count_error": lightcurve_data["broad"]["COUNT_RATE_ERR"],
        "ultrasoft_count_rate": lightcurve_data["ultrasoft"]["COUNT_RATE"],
        "soft_count_rate": lightcurve_data["soft"]["COUNT_RATE"],
        "medium_count_rate": lightcurve_data["medium"]["COUNT_RATE"],
        "hard_count_rate": lightcurve_data["hard"]["COUNT_RATE"],
        "ultrasoft_counts": lightcurve_data["ultrasoft"]["COUNTS"],
        "soft_counts": lightcurve_data["soft"]["COUNTS"],
        "medium_counts": lightcurve_data["medium"]["COUNTS"],
        "hard_counts": lightcurve_data["hard"]["COUNTS"],
        "exposure": lightcurve_data["broad"]["EXPOSURE"],
        "area": lightcurve_data["broad"]["AREA"],
    })
    file_relative0, file_relative1 = file_relativeS
    output_csv = os.path.join(graphs_directory, file_relative0, str(file_relative1)+'.csv')
    combined_data.to_csv(output_csv, index=False)
    return output_csv




def create_plot(lightcurve_data: dict[str, DataFrame], binsize, file_relativeS, graphs_directory):
    """Generate a plt plot to model the lightcurves."""
    matplotlib.use("svg")

    # chandra time initialization
    chandra_mjd_ref = 50814.0
    time_seconds = lightcurve_data["broad"]["TIME"]
    initial_time = time_seconds.min()
    final_time = time_seconds.max()
    initial_time_days = initial_time / 86400.0
    observation_mjd = chandra_mjd_ref + initial_time_days
    observation_date = Time(observation_mjd, format='mjd').to_datetime()
    readable_date = observation_date.strftime('%Y-%m-%d %H:%M:%S')
    # integer counts
    integer_counts = lightcurve_data["broad"]["COUNTS"].round().astype(int) 
    # time_kiloseconds = time_seconds / 1000
    zero_shifted_time_kiloseconds = (time_seconds - initial_time) / 1000
    observation_duration = zero_shifted_time_kiloseconds.max()
    # glvary_times, glvary_probs = self.get_glvary_data(glvary_file)
    width = 12 * (500 / binsize if binsize < 500 else 1)
    nrows = 7
    # figure 
    figure, (broad_plot, separation_plot, counts_plot, hr_plot, cumulative_counts_plot, bayesian_blocks_plot, lomb_scargle_plot) = plt.subplots(
        nrows=nrows, ncols=1, figsize=(width, nrows*3), constrained_layout=True
    )

    # figure, (broad_plot, separation_plot, counts_plot, hr_plot, bayesian_blocks_plot, lomb_scargle_plot, glvary_plot) = plt.subplots(
    #     nrows=7, ncols=1, figsize=(12, 21), constrained_layout=True
    # )

    # count rate plot 
    broad_plot.errorbar(
        x=zero_shifted_time_kiloseconds,
        y=lightcurve_data["broad"]["COUNT_RATE"],
        yerr=lightcurve_data["broad"]["COUNT_RATE_ERR"],
        color="red",
        marker="s",
        markerfacecolor="black",
        markersize=4,
        ecolor="black",
        markeredgecolor="black",
        capsize=3,
    )        
    broad_plot.set_xlim([0, observation_duration])
    broad_plot.set_title("Broadband Count Rate", fontsize=14)
    broad_plot.set_ylabel("Count Rate (counts/s)", fontsize=12)
    broad_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
    broad_plot.xaxis.set_major_locator(MultipleLocator(5))  
    broad_plot.xaxis.set_minor_locator(MultipleLocator(1))
    broad_plot.tick_params(axis='both', which='major', labelsize=10)
    broad_plot.text(0.995, 1.13, f"Start: {readable_date}",
        transform=broad_plot.transAxes,
        fontsize=10, ha='right', va='top', bbox=dict(facecolor='white', alpha=0.7))

    # separated light plot
    separated_light_level_colors = {"ultrasoft": "green", "soft": "red", "medium": "gold", "hard": "blue"}

    # loop through light levels
    for light_level, color in separated_light_level_colors.items():
        separation_plot.plot(
            zero_shifted_time_kiloseconds,
            lightcurve_data[light_level]["COUNT_RATE"],
            color=color,
            label=light_level.capitalize(),
        )
    separation_plot.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.19),
        ncol=4,
        frameon=False,
        fontsize=12,
    )
    separation_plot.set_xlim([0, observation_duration])
    separation_plot.set_title("Separated Energy Band Count Rates", fontsize=14, y = 1.15)
    separation_plot.set_ylabel("Count Rate (counts/s)", fontsize=12)
    separation_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
    separation_plot.xaxis.set_major_locator(MultipleLocator(5))  
    separation_plot.xaxis.set_minor_locator(MultipleLocator(1))
    separation_plot.tick_params(axis='both', which='major', labelsize=10)
    
    # counts plot
    counts_plot.plot(
        zero_shifted_time_kiloseconds,
        lightcurve_data["broad"]["COUNTS"],
        color="purple",
        label="Broadband Counts",
        marker='o',
        markersize=4,
    )
    counts_plot.set_xlim([0, observation_duration])
    counts_plot.set_title("Counts in Broadband", fontsize=14, y = 1.05)
    counts_plot.set_xlabel("Time (kiloseconds)", fontsize=12)
    counts_plot.set_ylabel("Counts", fontsize=12)
    counts_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
    counts_plot.xaxis.set_major_locator(MultipleLocator(5))  
    counts_plot.xaxis.set_minor_locator(MultipleLocator(1))
    counts_plot.tick_params(axis='both', which='major', labelsize=10)
    total_counts=round(float(lightcurve_data["broad"]["COUNTS"].sum()), 3)
    average_count_rate=round(float(lightcurve_data["broad"]["COUNT_RATE"].mean()), 3)
    min_count_rate=round(float(lightcurve_data["broad"]["COUNT_RATE"].min()), 3)
    max_count_rate=round(float(lightcurve_data["broad"]["COUNT_RATE"].max()), 3)

    text_str = (
        f"Total Counts: {total_counts:.3f}          Average Count Rate: {average_count_rate:.3f}\n"
        f"Min Count Rate: {min_count_rate:.3f}          Max Count Rate: {max_count_rate:.3f}"
    )
    counts_plot.text(0.995, 1.2, text_str,
                    transform=counts_plot.transAxes,
                    fontsize=10, ha='right', va='top', bbox=dict(facecolor='white', alpha=0.7))


    # calculate hardness ratios
    hard_counts = lightcurve_data["hard"]["COUNTS"]
    soft_counts = lightcurve_data["soft"]["COUNTS"]
    medium_counts = lightcurve_data["medium"]["COUNTS"]
    ultrasoft_counts = lightcurve_data["ultrasoft"]["COUNTS"]

    total_counts_hs = soft_counts + hard_counts
    hr_hs = (hard_counts - soft_counts) / total_counts_hs
    total_counts_ms = soft_counts + medium_counts
    hr_ms = (medium_counts - soft_counts) / total_counts_ms
    total_counts_smh = soft_counts + medium_counts + hard_counts
    hr_smh = (soft_counts - (medium_counts + hard_counts)) / total_counts_smh
    total_counts_mhsu = ultrasoft_counts + soft_counts + medium_counts + hard_counts
    hr_mhsu = ((medium_counts + hard_counts) - (soft_counts + ultrasoft_counts) ) / total_counts_mhsu

    # hardness ratio plot
    hr_plot.plot(
        zero_shifted_time_kiloseconds,
        hr_hs,
        color="blue",
        label="HR_HS (H-S)/(H+S)",
        marker='^',
        markersize=4,
    )

    hr_plot.plot(
        zero_shifted_time_kiloseconds,
        hr_ms,
        color="orange",
        label="HR_MS (M-S)/(M+S)",
        marker='v',
        markersize=4,
    )

    hr_plot.plot(
        zero_shifted_time_kiloseconds,
        hr_smh,
        color="green",
        label="HR_SMH (S-(M+H))/(S+M+H)",
        marker='o',
        markersize=4,
    )

    hr_plot.plot(
        zero_shifted_time_kiloseconds,
        hr_mhsu,
        color="purple",
        label="HR_HMSU ((M+H)-(S+U)))/(M+H+S+U)",
        marker='1',
        markersize=4,
    )

    hr_plot.set_xlim([0, observation_duration])
    hr_plot.set_title("Hardness Ratios", fontsize=14, y = 1.29)
    hr_plot.set_xlabel("Time (kiloseconds)", fontsize=12)
    hr_plot.set_ylabel("Hardness Ratio", fontsize=12)
    hr_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
    hr_plot.xaxis.set_major_locator(MultipleLocator(5))
    hr_plot.xaxis.set_minor_locator(MultipleLocator(1))
    hr_plot.tick_params(axis='both', which='major', labelsize=10)
    hr_plot.legend(loc="upper center", bbox_to_anchor=(0.5, 1.32), ncol=2, frameon=False, fontsize=12)

    # jimmy's cumulative counts plot
    cumulative_counts = np.cumsum(integer_counts)
    cumulative_counts_plot.plot(zero_shifted_time_kiloseconds, cumulative_counts, color='magenta')
    cumulative_counts_plot.set_xlim([0, observation_duration])
    cumulative_counts_plot.set_title("Cumulative Counts", fontsize=14)
    cumulative_counts_plot.set_xlabel("Time (kiloseconds)", fontsize=12)
    cumulative_counts_plot.set_ylabel("Cumulative Counts", fontsize=12)
    cumulative_counts_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
    cumulative_counts_plot.xaxis.set_major_locator(MultipleLocator(5))
    cumulative_counts_plot.xaxis.set_minor_locator(MultipleLocator(1))
    cumulative_counts_plot.tick_params(axis='both', which='major', labelsize=10)

    # bayesian blocks plot
    bayesian_count_rate = lightcurve_data["broad"]["COUNT_RATE"]
    bins = bayesian_blocks(zero_shifted_time_kiloseconds, bayesian_count_rate, fitness='measures')
    median_bin_size = np.median(np.diff(bins))

    # Binned count rate light curve
    # binned_time = np.arange(zero_shifted_time_kiloseconds.min(), zero_shifted_time_kiloseconds.max(), median_bin_size)
    # binned_count_rate, _ = np.histogram(zero_shifted_time_kiloseconds, bins=binned_time, weights=bayesian_count_rate)
    # binned_count_rate_err, _ = np.histogram(zero_shifted_time_kiloseconds, bins=binned_time, weights=lightcurve_data["broad"]["COUNT_RATE_ERR"])

    # bayesian_blocks_plot.set_xlim([0, observation_duration])
    # bayesian_blocks_plot.hist(zero_shifted_time_kiloseconds, bins=bins, color='lightblue', edgecolor='black', alpha=0.7)
    # bayesian_blocks_plot.set_title("Bayesian Blocks Segmentation", fontsize=14)
    # bayesian_blocks_plot.set_xlabel("Time (seconds)", fontsize=12)
    # bayesian_blocks_plot.set_ylabel("Counts Rate", fontsize=12)
    # bayesian_blocks_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
    # bayesian_blocks_plot.xaxis.set_major_locator(MultipleLocator(5))
    # bayesian_blocks_plot.xaxis.set_minor_locator(MultipleLocator(1))
    # bayesian_blocks_plot.tick_params(axis='both', which='major', labelsize=10)

    # # Overplot binned count rate light curve
    # bayesian_blocks_plot.errorbar(
    #     x=binned_time[:-1] + median_bin_size/2,  # center of bins
    #     y=binned_count_rate,
    #     yerr=binned_count_rate_err,
    #     color="red",
    #     marker="s",
    #     markerfacecolor="black",
    #     markersize=4,
    #     ecolor="black",
    #     markeredgecolor="black",
    #     capsize=3,
    # )
    # bayesian blocks plot
    bins = bayesian_blocks(zero_shifted_time_kiloseconds, integer_counts, fitness='measures')
    bayesian_blocks_plot.set_xlim([0, observation_duration])
    bayesian_blocks_plot.hist(zero_shifted_time_kiloseconds, bins=bins, color='lightblue', edgecolor='black', alpha=0.7)
    bayesian_blocks_plot.set_title("Bayesian Blocks Segmentation", fontsize=14)
    bayesian_blocks_plot.set_xlabel("Time (kiloseconds)", fontsize=12) ## OK EDIT TO BE KILOSECONDS
    bayesian_blocks_plot.set_ylabel("Counts Rate", fontsize=12)
    bayesian_blocks_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
    bayesian_blocks_plot.xaxis.set_major_locator(MultipleLocator(5))
    bayesian_blocks_plot.xaxis.set_minor_locator(MultipleLocator(1))
    bayesian_blocks_plot.tick_params(axis='both', which='major', labelsize=10)

    # # lomb scargle periodicity plot
    # frequency, power = LombScargle(zero_shifted_time_kiloseconds, lightcurve_data["broad"]["COUNTS"]).autopower()
    # period = 1 / frequency
    # lomb_scargle_plot.set_xlim([0, 1 / frequency.max()])
    # lomb_scargle_plot.plot(period, power, color='darkblue')
    # lomb_scargle_plot.set_title("Lomb-Scargle Periodogram", fontsize=14)
    # lomb_scargle_plot.set_xlabel("Period (sec)", fontsize=12)
    # lomb_scargle_plot.set_ylabel("Power", fontsize=12)
    # lomb_scargle_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
    # lomb_scargle_plot.xaxis.set_major_locator(MultipleLocator(0.2))
    # lomb_scargle_plot.xaxis.set_minor_locator(MultipleLocator(0.1))
    # lomb_scargle_plot.tick_params(axis='both', which='major', labelsize=10)

    # lomb scargle periodicity plot
    frequency, power = LombScargle(zero_shifted_time_kiloseconds, integer_counts).autopower()
    lomb_scargle_plot.set_xlim([0, frequency.max()])
    lomb_scargle_plot.plot(frequency, power, color='darkblue')
    lomb_scargle_plot.set_title("Lomb-Scargle Periodogram", fontsize=14)
    lomb_scargle_plot.set_xlabel("Frequency (1/sec)", fontsize=12)
    lomb_scargle_plot.set_ylabel("Power", fontsize=12)
    lomb_scargle_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
    lomb_scargle_plot.xaxis.set_major_locator(MultipleLocator(0.2))
    lomb_scargle_plot.xaxis.set_minor_locator(MultipleLocator(0.1))
    lomb_scargle_plot.tick_params(axis='both', which='major', labelsize=10)

    # glvary_plot.plot(
    #     glvary_times, 
    #     glvary_probs, 
    #     color='darkgreen', 
    #     label='GLVAR Probabilities'
    # )
    # glvary_plot.set_xlim([0, observation_duration])
    # glvary_plot.set_title("Gregory-Laredo Variability Analysis", fontsize=14)
    # glvary_plot.set_xlabel("Time (kiloseconds)", fontsize=12)
    # glvary_plot.set_ylabel("Probability", fontsize=12)
    # glvary_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
    # glvary_plot.xaxis.set_major_locator(MultipleLocator(5))
    # glvary_plot.xaxis.set_minor_locator(MultipleLocator(1))
    # glvary_plot.tick_params(axis='both', which='major', labelsize=10)

    figure.suptitle(
        f"Lightcurve in Broadband and Separated Energy Bands (Binsize of {binsize}s)"
    )
    # plt.savefig(svg_data := StringIO(), bbox_inches="tight")
    file_relative0, file_relative1 = file_relativeS
    # plt.savefig(os.path.join(graphs_directory, file_relative0, str(file_relative1)+'.svg'),  bbox_inches="tight")
    plt.savefig(os.path.join(graphs_directory, file_relative0, str(file_relative1)+'.png'),  bbox_inches="tight")
    plt.close(figure)
    # return svg_data







## NEW FUNCTIONS


def filter_low_count(lightcurve_data, cutoff):
    """
    Returns average count rate if it's above cutoff, otherwise returns "small".
    """
    multipl = np.multiply(np.array(lightcurve_data["broad"]["COUNT_RATE"]), np.array(lightcurve_data["broad"]["EXPOSURE"]))
    mean_cr = np.sum(multipl)/(lightcurve_data["broad"]["EXPOSURE"].sum())
    if mean_cr<cutoff:
        return 'small'
    return mean_cr

def get_intervals(dips, max_len):
    """
    Returns positions and sizes of dips, excluding those at edges of the lightcurves.
    """
    dip_positions = []

    end_i = 0
    start_i = 0
    for i in dips: 
        if i == end_i + 1: #we made a step of size 1, continue building interval
            end_i = i
        else: # interval is done
            if start_i != 0: # don't consider any interval that starts from start (could also be 0,0 being added if first dip starts later)
                dip_positions += [[start_i, end_i]]
            start_i = i
            end_i = i
    
    if end_i != max_len and start_i != 0: # don't consider any interval that reaches exact end
        dip_positions += [[start_i, end_i]] # add last interval
    

    if dip_positions == []:
        return [], None
    dip_lengths = [d[1]-d[0]+1 for d in dip_positions]
    return dip_positions, dip_lengths

    
def get_intervals_all(dips):
    """
    Returns positions and sizes of dips, including those at edges of the lightcurves.
    """
    dip_positions = []

    end_i = 0
    start_i = 0
    for i in dips: 
        if i == end_i + 1: #we made a step of size 1, continue building interval
            end_i = i
        else: # interval is done
            if end_i != 0: # don't consider 0,0 being added, do allow dips to start from 0
                dip_positions += [[start_i, end_i]]
            start_i = i
            end_i = i

    dip_positions += [[start_i, end_i]] # add last interval
    
    if dip_positions == []:
        return [], None 
    dip_lengths = [d[1]-d[0]+1 for d in dip_positions]
    return dip_positions, dip_lengths


def automatic_dip_flare(lightcurve_data, cutoff_dip, cutoff_flare, mean_cr):
        """
        Determines stretches of points on the lightcurve consistently below the dip cutoff 
        or above the flare cutoff (as measured in percent deviation from the average count rate).
        
        For each dip/flare, determines 4 values: 

            1. The probability that a deviation of at least X length occured randomly assuming
            that the percent deviation from the mean was normally distributed with the same std dev as in 
            the real data. This is a conservative estimate of the probability that the dip/flare occured due to 
            randomness (i.e, type 1 error, for the hypothesis: 
            H0: There is no dip in this lightcurve
            vs.  
            HA: There is a dip somewhere in this lightcurve.
            )
            This analysis is done with a simulation. Analytical calculation is not feasible.

            2. Same as above (#1), but simulates the number of photons/bin (from there dividing by bin size to calculate 
            photon count rate and then calculating deviation from mean) with Poisson where 
            lambda = (total photons)/(total time) * (bin_size)

            3 & 4. Determines, "Assuming this is a dip/flare, how lucky were we to find it?": excludes points below the 
            dip/flare cutoffs and recalculates the distribution average and std for Normal- and Poisson-based simulations.
        
        """
        percent_diff = np.array((lightcurve_data["broad"]['COUNT_RATE']-mean_cr)/mean_cr)
        std_diff = np.std(percent_diff) #we consider each bin to contribute equally, regardless of it's size to the distribution


        dips = []
        flares = []
        
        for j in range(len(percent_diff)):
            if percent_diff[j]<-cutoff_dip:
                dips += [j]
            elif percent_diff[j]>cutoff_flare:
                flares += [j]
        
        if len(dips) == 0 and len(flares) == 0: 
            return None, None, None, None
        
        if len(dips) != 0:
            dip_positions,  dip_lengths = get_intervals(dips, len(percent_diff))
        else: 
            dip_positions = []
            dip_lengths = None

        if len(flares) != 0:
            flare_positions, flare_lengths = get_intervals(flares, len(percent_diff))
        else: 
            flare_positions = []
            flare_lengths = None


        no_dipflare_data = np.delete(lightcurve_data["broad"]["COUNT_RATE"], dip_positions + flare_positions, axis=0)
        no_dipflare_exposures = np.delete(np.array(lightcurve_data["broad"]["EXPOSURE"]), dip_positions + flare_positions, axis=0)
        mean_data_no_dipflare = np.sum(np.multiply(no_dipflare_data, no_dipflare_exposures))/(np.sum(no_dipflare_exposures))
        
        no_dipflare_diff = np.array((no_dipflare_data-mean_data_no_dipflare)/mean_data_no_dipflare)
        no_dipflare_std_diff = np.std(no_dipflare_diff)
        
        s_all_N = np.random.normal(0, std_diff, ITERATION_SIM*len(percent_diff)) # simulation data normal, all data
        s_part_N = np.random.normal(0, no_dipflare_std_diff, ITERATION_SIM*len(percent_diff)) # simulation data normal, no dip/flare data

        mode_bin_size =  stats.mode(lightcurve_data["broad"]["EXPOSURE"])[0]


        poisson_all_cr = np.random.poisson(lam=mean_cr*mode_bin_size, size=ITERATION_SIM*len(percent_diff))/mode_bin_size
        poisson_part_cr= np.random.poisson(lam=mean_data_no_dipflare*mode_bin_size, size=ITERATION_SIM*len(percent_diff))/mode_bin_size

        s_all_P = (poisson_all_cr-np.mean(poisson_all_cr))/np.mean(poisson_all_cr)
        s_part_P = (poisson_part_cr-np.mean(poisson_part_cr))/np.mean(poisson_part_cr)


        dip_lengths_all_N = []
        dip_lengths_part_N = []
        flares_lengths_all_N = []
        flares_lengths_part_N = []

        dip_lengths_all_P = []
        dip_lengths_part_P = []
        flares_lengths_all_P = []
        flares_lengths_part_P = []

        for xx in range(ITERATION_SIM):

            dips_all_N = []
            flares_all_N  = []

            dips_part_N = []
            flares_part_N  = []

            dips_all_P = []
            flares_all_P  = []

            dips_part_P = []
            flares_part_P  = []

            for j in range(len(percent_diff)):
                if s_all_N[xx*len(percent_diff)+j]<-cutoff_dip:
                    dips_all_N  += [j]
                elif s_all_N[xx*len(percent_diff)+j]>cutoff_flare:
                    flares_all_N += [j]

                if s_part_N[xx*len(percent_diff)+j]<-cutoff_dip:
                    dips_part_N  += [j]
                elif s_part_N[xx*len(percent_diff)+j]>cutoff_flare:
                    flares_part_N += [j]

                if s_all_P[xx*len(percent_diff)+j]<-cutoff_dip:
                    dips_all_P  += [j]
                elif s_all_P[xx*len(percent_diff)+j]>cutoff_flare:
                    flares_all_P += [j]

                if s_part_P[xx*len(percent_diff)+j]<-cutoff_dip:
                    dips_part_P  += [j]
                elif s_part_P[xx*len(percent_diff)+j]>cutoff_flare:
                    flares_part_P += [j]

            
            if len(dips_all_N) != 0:
                dip_lengths_all_N += [np.max(get_intervals_all(dips_all_N)[1])] # record len longest dip occured
            if len(dips_part_N) != 0:
                dip_lengths_part_N += [np.max(get_intervals_all(dips_part_N)[1])]
            if len(flares_all_N) != 0:
                flares_lengths_all_N += [np.max(get_intervals_all(flares_all_N)[1])]
            if len(dips_part_N) != 0:
                flares_lengths_part_N += [np.max(get_intervals_all(flares_part_N)[1])]

            if len(dips_all_P) != 0:
                dip_lengths_all_P += [np.max(get_intervals_all(dips_all_P)[1])] # record len longest dip occured
            if len(dips_part_P) != 0:
                dip_lengths_part_P += [np.max(get_intervals_all(dips_part_P)[1])]
            if len(flares_all_P) != 0:
                flares_lengths_all_P += [np.max(get_intervals_all(flares_all_P)[1])]
            if len(dips_part_P) != 0:
                flares_lengths_part_P += [np.max(get_intervals_all(flares_part_P)[1])]




        if len(dips) != 0:
            if len(dip_lengths_all_N) == 0: 
                prob_dip_lengths_all_N = np.zeros(len(dip_lengths))
            if len(dip_lengths_part_N) == 0: 
                prob_dip_lengths_part_N = np.zeros(len(dip_lengths))
        if len(flares) != 0:
            if len(flares_lengths_all_N) == 0: 
                prob_flares_lengths_all_N = np.zeros(len(flare_lengths))
            if len(flares_lengths_part_N) == 0: 
                prob_flares_lengths_part_N = np.zeros(len(flare_lengths))
            
        if len(dips) != 0:
            if len(dip_lengths_all_P) == 0: 
                prob_dip_lengths_all_P = np.zeros(len(dip_lengths))
            if len(dip_lengths_part_P) == 0: 
                prob_dip_lengths_part_P = np.zeros(len(dip_lengths))
        if len(flares) != 0:
            if len(flares_lengths_all_P) == 0: 
                prob_flares_lengths_all_P = np.zeros(len(flare_lengths))
            if len(flares_lengths_part_P) == 0: 
                prob_flares_lengths_part_P = np.zeros(len(flare_lengths))


        dip_lengths_all_N = np.array(dip_lengths_all_N)
        dip_lengths_part_N = np.array(dip_lengths_part_N)
        flares_lengths_all_N = np.array(flares_lengths_all_N)
        flares_lengths_part_N = np.array(flares_lengths_part_N)

        dip_lengths_all_P = np.array(dip_lengths_all_P)
        dip_lengths_part_P = np.array(dip_lengths_part_P)
        flares_lengths_all_P = np.array(flares_lengths_all_P)
        flares_lengths_part_P = np.array(flares_lengths_part_P)


        if dip_lengths: 
            prob_dip_lengths_all_N = []
            prob_dip_lengths_part_N = []
            prob_dip_lengths_all_P = []
            prob_dip_lengths_part_P = []
            for l in dip_lengths: 
                prob_dip_lengths_all_N += [np.sum(np.where(dip_lengths_all_N>=l, 1, 0))/ITERATION_SIM]
                prob_dip_lengths_part_N += [np.sum(np.where(dip_lengths_part_N>=l, 1, 0))/ITERATION_SIM]
                prob_dip_lengths_all_P += [np.sum(np.where(dip_lengths_all_P>=l, 1, 0))/ITERATION_SIM]
                prob_dip_lengths_part_P += [np.sum(np.where(dip_lengths_part_P>=l, 1, 0))/ITERATION_SIM]
        else: 
            prob_dip_lengths_all_N = None
            prob_dip_lengths_part_N = None
            prob_dip_lengths_all_P = None
            prob_dip_lengths_part_P = None

        if flare_lengths: 
            prob_flares_lengths_all_N = []
            prob_flares_lengths_part_N = []
            prob_flares_lengths_all_P = []
            prob_flares_lengths_part_P = []
            for l in flare_lengths: 
                prob_flares_lengths_all_N += [np.sum(np.where(flares_lengths_all_N>=l, 1, 0))/ITERATION_SIM]
                prob_flares_lengths_part_N += [np.sum(np.where(flares_lengths_part_N>=l, 1, 0))/ITERATION_SIM]
                prob_flares_lengths_all_P += [np.sum(np.where(flares_lengths_all_P>=l, 1, 0))/ITERATION_SIM]
                prob_flares_lengths_part_P += [np.sum(np.where(flares_lengths_part_P>=l, 1, 0))/ITERATION_SIM]
        else: 
            prob_flares_lengths_all_N = None
            prob_flares_lengths_part_N = None
            prob_flares_lengths_all_P = None
            prob_flares_lengths_part_P = None
        
        print(dip_positions, dip_lengths) 
        print(prob_dip_lengths_all_N, prob_dip_lengths_part_N, prob_dip_lengths_all_P, prob_dip_lengths_part_P)

        print(flare_positions, flare_lengths)
        print(prob_flares_lengths_all_N, prob_flares_lengths_part_N, prob_flares_lengths_all_P, prob_flares_lengths_part_P)

        return (dip_positions, dip_lengths), (prob_dip_lengths_all_N, prob_dip_lengths_part_N, prob_dip_lengths_all_P, prob_dip_lengths_part_P), (flare_positions, flare_lengths), (prob_flares_lengths_all_N, prob_flares_lengths_part_N, prob_flares_lengths_all_P, prob_flares_lengths_part_P)






Path(dipflare_directory).mkdir(parents=True, exist_ok=True)
bigenough = os.path.join(dipflare_directory, "bigenough")
dirdip = os.path.join(dipflare_directory, "dip")
dirflare = os.path.join(dipflare_directory, "flare")

Path(bigenough).mkdir(exist_ok=True)
Path(dirdip).mkdir(exist_ok=True)
Path(dirflare).mkdir(exist_ok=True)

if not os.path.exists(graphs_directory):
    os.makedirs(graphs_directory)


for subdir, dirs, files in os.walk(download_directory):
    # for file in files:
    directory = os.path.join(subdir)
    dirPath = Path(directory)

    file_relative = dirPath.relative_to(download_directory)

    file_relative = os.path.normpath(file_relative)
    file_relative_split = file_relative.split(os.sep)
    
    Path(os.path.join(graphs_directory, file_relative_split[0])).mkdir(exist_ok=True)

    filename = str(file_relative).replace('/', '_')
    print(filename)
    
    unzip_fits_files(dirPath)
    event_list_file, source_region_file = None, None
    for data_product_file in dirPath.glob("*.fits"):
        stem = data_product_file.stem.lower()
        if stem.endswith("regevt3"):
            event_list_file = data_product_file
        if stem.endswith("reg3"):
            source_region_file = data_product_file
    if not event_list_file or not source_region_file:
        continue # the directory was not deep enough
    
    print('cleaning: ', directory)
    source_region = Path(source_region_file)
    event_list = Path(event_list_file)

    try: 
        region_event_list = isolate_source_region(event_list, source_region)
        lightcurves = extract_lightcurves(region_event_list, bin_size)
        filtered_lightcurves = filter_lightcurve_columns(lightcurves)
        lightcurve_data = get_lightcurve_data(filtered_lightcurves)
    except:
        print('data does not have photon energy, skipping')
        continue
        # PREVIOUS ERROR: dmextract (CIAO 4.16.0): 
        # #dsDMBLOCKOPENERR -- ERROR: Failed to open the DM data block in dataset 'datadir_small_cloud_000500-723000_rad30_incompl/2CXOJ005407.2-722522/07279_000/hrcf07279_000N020_r0043_regevt3.src.fits[energy=300:7000]'.  
        # DM Error: 'Error: Could not find identifier energy.'.
    
 
    filter_result = filter_low_count(lightcurve_data, min_avg_count_rate)
    
    if filter_result != 'small':
        print('analyzing')
        file2write=open(os.path.join(bigenough, str(filename)+'.txt'),'w')
        file2write.write(str(filename))
        file2write.close()
        
        dips, dipsprob, flares, flaresprob = automatic_dip_flare(lightcurve_data, dip_threshold, flare_threshold, filter_result)


        if dips: 
            if dips[1]:
                ## significance threshold considers poisson with data removed (generally the lowest one)
                ## considers longest interval (min prob)
                if min(dipsprob[3]) < significance_of_dipflare_threshold:
                    print('dip')

                    file2write=open(os.path.join(dirdip, str(min(dipsprob[3]))+'__'+str(filename)+'.txt'),'w')
                    file2write.write(str(filename)+ '\n')
                    file2write.write('DIP POSITIONS | DIP LENGTHS' + '\n')
                    file2write.write(str(dips[0]) + '  |  ' + str(dips[1]) + '\n')
                    file2write.write('DIP PROBABILITIES: 1) Normal all data, 2) Normal no dip/flare, 3) Poisson all data, 4) Poisson no dip/flare' + '\n')
                    for f in range(4): 
                        file2write.write(str(dipsprob[f])+ '\n')
                    
                    file2write.close()

        if flares: 
            if flares[1]:
                if min(flaresprob[3]) < significance_of_dipflare_threshold:
                    print('flare')

                    file2write=open(os.path.join(dirflare, str(min(flaresprob[3]))+'__'+str(filename)+'.txt'),'w')
                    file2write.write(str(filename)+ '\n')
                    file2write.write('FLARE POSITIONS | FLARE LENGTHS' + '\n')
                    file2write.write(str(flares[0]) + '  |  ' + str(flares[1]) + '\n')
                    file2write.write('FLARE PROBABILITIES: 1) Normal all data, 2) Normal no dip/flare, 3) Poisson all data, 4) Poisson no dip/flare' + '\n')
                    for f in range(4): 
                        file2write.write(str(flaresprob[f])+ '\n')
                    
                    file2write.close()


    else: 
        print('too small')

    if create_just_big: 
        if filter_result != 'small':
            create_csv(lightcurve_data, file_relative_split, graphs_directory)
            create_plot(lightcurve_data, bin_size, file_relative_split, graphs_directory)    
    else:
        create_csv(lightcurve_data, file_relative_split, graphs_directory)
        create_plot(lightcurve_data, bin_size, file_relative_split, graphs_directory)


