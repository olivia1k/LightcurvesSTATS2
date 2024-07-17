"""Mihir Patankar [mpatankar06@gmail.com]"""
import uuid
from abc import ABC, abstractmethod
from io import StringIO
from pathlib import Path
import threading

import numpy as np

import matplotlib
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from astropy import io, table
from astropy.timeseries import LombScargle
from astropy.stats import bayesian_blocks
from astropy.time import Time

# pylint: disable-next=no-name-in-module
from ciao_contrib.runtool import dmcopy, dmextract, dmkeypar, dmlist, dmstat, new_pfiles_environment, glvary
from matplotlib import pyplot as plt
from pandas import DataFrame

from data_structures import LightcurveParseResults, Message, ObservationData, ObservationHeaderInfo
from postage_stamp_plotter import CropBounds, plot_postagestamps

## OK ADD
import numpy as np
import os
ITERATION_SIM = 10000
## OK ADD END

class ObservationProcessor(ABC):
    """Base class for observation processor implementations for different Chandra instruments."""

    def __init__(self, data_products, binsize, source_name, dirlabels, cutoff_dip, cutoff_flare, significance_threshold_count_rate, cutoff, message_collection_queue=None, counts_checker=None):
        self.event_list = Path(data_products.event_list_file)
        self.source_region = Path(data_products.source_region_file)
        self.detector_coords_image, self.sky_coords_image = None, None
        self.message_collection_queue = message_collection_queue
        self.counts_checker = counts_checker
        self.binsize = binsize

        ## OK ADD PARAMETERS + CREATE DIRS
        self.source_name = source_name
        self.cutoff = cutoff
        self.cutoff_dip = cutoff_dip
        self.cutoff_flare = cutoff_flare
        self.significance_threshold_count_rate = significance_threshold_count_rate
        self.dirlabels = dirlabels
        Path(self.dirlabels).mkdir(parents=True, exist_ok=True)

        self.bigenough = os.path.join(self.dirlabels, "bigenough")
        Path(self.bigenough).mkdir(exist_ok=True)
        self.dirdip = os.path.join(self.dirlabels,"dip")
        Path(self.dirdip).mkdir(exist_ok=True)
        self.dirflare = os.path.join(self.dirlabels,"flare")
        Path(self.dirflare).mkdir(exist_ok=True)
        

    def process(self):
        """Sequence in which all steps of the processing routine are called."""

        message_uuid = uuid.uuid4()
        with new_pfiles_environment():
            observation_id = dmkeypar(infile=f"{self.event_list}", keyword="OBS_ID", echo=True)
            prefix = f"Observation {observation_id}: "

            filename = str(self.source_name)+'_'+str(observation_id)

            def status(status):
                self.message_collection_queue.put(Message(f"{prefix}{status}", message_uuid))

            status("Isolating source region...")
            region_event_list = self.isolate_source_region(self.event_list, self.source_region)
            status("Extracting lightcurves...")
            lightcurves = self.extract_lightcurves(region_event_list, self.binsize)
            status("Copying columns...")
            filtered_lightcurves = self.filter_lightcurve_columns(lightcurves)

            status("Checking counts...")
            lightcurve_data = self.get_lightcurve_data(filtered_lightcurves)
            status('Got Lightcurve data')
            self.counts_checker.queue.put(self.get_lightcurve_counts(lightcurve_data))
            self.counts_checker.queue.join()
            if self.counts_checker.cancel_event.is_set():
                return None

## OK ADD
            status('Preparing to filter')
            filter_result = self.filter_low_count(lightcurve_data, self.cutoff)
            status("filter done "+str(filter_result))
            if filter_result != 'small':
                status("not small")
                file2write=open(self.bigenough + "/" +str(filename)+'.txt','w')
                file2write.write(str(filename))
                file2write.close()
                status("Dip/flare detection...")
                dipflare = self.automatic_dip_flare(lightcurve_data, self.cutoff_dip, self.cutoff_flare, filter_result)
                status("Dip/flare done")
                prob_dip_random = float(dipflare[2])
                prob_flare_random = float(dipflare[3])   

                if dipflare[0] > 0 and prob_dip_random < self.significance_threshold_count_rate:
                    file2write=open(self.dirdip + "/" +str(prob_dip_random)+'__'+str(filename)+'.txt','w')
                    file2write.write(str(prob_dip_random) + '\n')
                    file2write.write(str(dipflare[0])+ '\n')
                    file2write.write(str(filename))
                    file2write.close()

                if dipflare[1] >0 and prob_flare_random < self.significance_threshold_count_rate:
                    file2write=open(self.dirflare + "/" +str(prob_flare_random)+'__'+str(filename)+'.txt','w')
                    file2write.write(str(prob_flare_random)+ '\n')
                    file2write.write(str(dipflare[1])+ '\n')
                    file2write.write(str(filename))
                    file2write.close()

## OK ADD END   
            
            # status("Checking counts...")
            # lightcurve_data = self.get_lightcurve_data(filtered_lightcurves)
            # # Function to check counts
            # def check_counts():
            #     self.counts_checker.queue.put(self.get_lightcurve_counts(lightcurve_data))
            #     self.counts_checker.queue.join()

            # # Create and start the counts checking thread
            # counts_thread = threading.Thread(target=check_counts)
            # counts_thread.start()

            # # Wait for up to 10 seconds for the counts checking to complete
            # counts_thread.join(timeout=30)

            # # Check if the thread is still running after the timeout
            # if counts_thread.is_alive():
            #     status("Counts checking exceeded 30 seconds. Cancelling operation.")
            #     self.counts_checker.cancel_event.set()
            #     return None

            # # Check if the cancel event was set by the counts checker
            # if self.counts_checker.cancel_event.is_set():
            #     return None

            status("Retrieving images...")
            self.get_images(region_event_list)
            status("Plotting lightcurves...")
            results = self.plot(lightcurve_data)
            status("Plotting lightcurves... Done")

        return results

    @abstractmethod
    def extract_lightcurves(self, event_list, binsize):
        """Extract lightcurve(s) from an event list, one should pass in one with a specific source
        region extracted."""

    @staticmethod
    @abstractmethod
    def filter_lightcurve_columns(lightcurves: list[Path]):
        """Filter lightcurve(s) to get the columns we care about and format them."""

    @staticmethod
    @abstractmethod
    def get_lightcurve_data(lightcurves: list[Path]):
        """Return the data from the lightcurve files in a python object."""

    @staticmethod
    def get_lightcurve_counts(lightcurve_data):
        """Return the total counts in a lightcurve so they can be verified to meet the threshold."""

    @abstractmethod
    def plot(self, lightcurve_data) -> LightcurveParseResults:
        """Plots lightcurve(s). Returns what will then be returned by the thread pool future."""

    @staticmethod
    def isolate_source_region(event_list: Path, source_region):
        """Restrict the event list to just the source region of the source we care about."""
        dmcopy(
            infile=f"{event_list}[sky=region({source_region})]",
            outfile=(outfile := f"{event_list.with_suffix('.src.fits')}"),
            clobber="yes",
        )
        return outfile

    def get_images(self, region_event_list):
        """Gets the images from the event list, one in sky coordinates, the other in detector
        coordinates. This is useful to be able to track if a lightcurve dip corresponds with a
        source going over the edge of the detector. The image is cropped otherwise we would be
        dealing with hundreds of thousands of blank pixels. The cropping bounds are written to the
        FITS file for later use in plotting."""

        dmstat(infile=f"{region_event_list}[cols x,y]")
        sky_bounds = CropBounds.from_strings(*dmstat.out_min.split(","), *dmstat.out_max.split(","))
        sky_bounds.double()
        dmcopy(
            infile=f"{self.event_list}"
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
            infile=f"{self.event_list}"
            f"[bin detx={detector_bounds.x_min}:{detector_bounds.x_max}:0.5,"
            f"dety={detector_bounds.y_min}:{detector_bounds.y_max}:0.5]",
            outfile=(detector_coords_image := f"{region_event_list}.detimg.fits"),
        )
        with io.fits.open(detector_coords_image, mode="append") as hdu_list:
            hdu_list.append(detector_bounds.to_hdu())
        self.sky_coords_image, self.detector_coords_image = sky_coords_image, detector_coords_image

    def get_observation_details(self):
        """Gets keys from the header block detailing the observation information."""
        return ObservationHeaderInfo(
            instrument=dmkeypar(infile=f"{self.event_list}", keyword="INSTRUME", echo=True),
            observation_id=dmkeypar(infile=f"{self.event_list}", keyword="OBS_ID", echo=True),
            region_id=dmkeypar(infile=f"{self.event_list}", keyword="REGIONID", echo=True),
            start_time=dmkeypar(infile=f"{self.event_list}", keyword="DATE-OBS", echo=True),
            end_time=dmkeypar(infile=f"{self.event_list}", keyword="DATE-END", echo=True),
        )


## OK ADD
    @staticmethod
    def filter_low_count(lightcurve_data, cutoff):
        """
        Returns average count rate if it's above cutoff, otherwise returns "small".
        """
        multipl = np.multiply(np.array(lightcurve_data["broad"]["COUNT_RATE"]), np.array(lightcurve_data["broad"]["EXPOSURE"]))
        mean_cr = np.sum(multipl)/(lightcurve_data["broad"]["EXPOSURE"].sum())
        if mean_cr<cutoff:
            return 'small'
        return mean_cr

    @staticmethod
    def automatic_dip_flare(lightcurve_data, cutoff_dip, cutoff_flare, mean_cr):
        """
        Determines the longest number of points on the lightcurve consistently below the dip cutoff 
        or above the flare cutoff (as measured in percent deviation from the average count rate).
        Determines the probability that a deviation of at least that length occured randomly in the event 
        that the percent deviation from the mean was normally distributed with the same std dev as in 
        the real data. This is a conservative estimate of the probability that the dip/flare occured due to 
        randomness (i.e, type 1 error, for the hypothesis: 
        H0: There is no dip in this lightcurve
        vs.  
        HA: There is a dip somewhere in this lightcurve.)
        """
        percent_diff = np.array((lightcurve_data["broad"]['COUNT_RATE']-mean_cr)/mean_cr)
        std_diff = np.std(percent_diff) #we consider each bin to contribute equally, regardless of it's size to the distribution
               
        k = 0
        k_max = 0 # longest dip
        f = 0
        f_max = 0 # longest flare
        for j in range(len(percent_diff)):
            if percent_diff[j]<-cutoff_dip:
                k+=1
                if k>k_max: 
                    k_max = k
                else: 
                    k=0
            if percent_diff[j]>cutoff_flare:
                f+=1
                if f>f_max: 
                    f_max = f
                else: 
                    f=0

        if k_max == 0 and f_max ==0: # no flare or dip
            return 0, 0, 0, 0

        s = np.random.normal(0, std_diff, ITERATION_SIM*len(percent_diff)) # simulation data

        accidental_dip = 0
        accidental_flare = 0

        for xx in range(ITERATION_SIM):
            Tk = 0
            Tk_max = 0
            Tf = 0
            Tf_max = 0
            stop_dip = False
            stop_flare = False
            for j in range(len(percent_diff)): 

                if not stop_dip and s[xx*len(percent_diff)+j]<-cutoff_dip:
                    Tk+=1
                    if Tk>Tk_max: 
                        Tk_max = Tk
                    else: 
                        Tk=0
                if not stop_flare and s[xx*len(percent_diff)+j]>cutoff_flare:
                    Tf+=1
                    if Tf>Tf_max: 
                        Tf_max = Tf
                    else: 
                        Tf=0

                if not stop_dip and Tk_max>=k_max:
                    accidental_dip+=1
                    stop_dip = True
                if not stop_flare and Tf_max>=f_max:
                    accidental_flare+=1
                    stop_flare = True

                if stop_dip and stop_flare: 
                    continue

        return k_max, f_max, accidental_dip/ITERATION_SIM, accidental_flare/ITERATION_SIM


## OK ADD END





class AcisProcessor(ObservationProcessor):
    """Processes observations produced by the ACIS (Advanced CCD Imaging Spectrometer) instrument
    aboard Chandra."""

    ENERGY_LEVELS = {
        "broad": "energy=150:7000",
        "ultrasoft": "energy=150:300",
        "soft": "energy=300:1200",
        "medium": "energy=1200:2000",
        "hard": "energy=2000:7000",
    }

    @staticmethod
    def adjust_binsize(event_list, binsize):
        """For ACIS, time resolution can be in the seconds in timed exposure mode, as compared to in
        the microseconds for HRC. Thus we must round the binsize to the time resolution."""
        time_resolution = float(dmkeypar(infile=str(event_list), keyword="TIMEDEL", echo=True))
        return binsize // time_resolution * time_resolution

    def extract_lightcurves(self, event_list, binsize):
        outfiles = []
        self.binsize = self.adjust_binsize(event_list, binsize)
        for light_level, energy_range in AcisProcessor.ENERGY_LEVELS.items():
            dmextract(
                infile=f"{event_list}[{energy_range}][bin time=::" f"{self.binsize}]",
                outfile=(outfile := f"{event_list}.{light_level}.lc"),
                opt="ltc1",
                clobber="yes",
            )
            outfiles.append(Path(outfile))
        return outfiles
    
    # def extract_glvary(self, event_list, binsize, effile):
    #     outfiles = []
    #     self.binsize = self.adjust_binsize(event_list, binsize)
        
    #     # Define the infile and outfile names
    #     infile = f"{event_list}"
    #     outfile = f"{event_list}.gl_prob.fits"
    #     lcfile = f"{event_list}.lc_prob.fits"
        
    #     # Call the glvary tool
    #     glvary(
    #         infile=infile,
    #         outfile=outfile,
    #         lcfile=lcfile,
    #         effile=effile,  # path to the efficiency file
    #         clobber="yes"
    #     )
        
    #     # Append the generated lightcurve file to the outfiles list
    #     outfiles.append(Path(lcfile))
    #     outfiles.append(Path(outfile))
        
    #     return outfiles
    
    # @staticmethod
    # def get_glvary_data(glvary_file):
    #     with io.fits.open(glvary_file) as hdul:
    #         glvary_data = hdul[1].data
    #         times = glvary_data['TIME']
    #         probabilities = glvary_data['PROB']
    #     return times, probabilities

    @staticmethod
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

    @staticmethod
    def get_lightcurve_data(lightcurves: list[Path]):
        lightcurve_data: dict[str, DataFrame] = {
            energy_level: table.Table.read(lightcurve, format="ascii").to_pandas()
            for energy_level, lightcurve in zip(AcisProcessor.ENERGY_LEVELS.keys(), lightcurves)
        }
        # Trim zero exposure points
        for energy_level, lightcurve_dataframe in lightcurve_data.items():
            lightcurve_data[energy_level] = lightcurve_dataframe[
                lightcurve_dataframe["EXPOSURE"] != 0
            ]
        return lightcurve_data

    # def get_lightcurve_data(lightcurves: list[Path]):
    #     lightcurve_data: dict[str, DataFrame] = {}
    #     for energy_level, lightcurve in zip(AcisProcessor.ENERGY_LEVELS.keys(), lightcurves):
    #         try:
    #             # Attempt to read the file and convert to pandas DataFrame
    #             table_data = table.Table.read(lightcurve, format="ascii")
    #             lightcurve_dataframe = table_data.to_pandas()
    #             # Ensure the DataFrame has an "EXPOSURE" column
    #             if "EXPOSURE" in lightcurve_dataframe.columns:
    #                 # Filter out rows with zero exposure
    #                 lightcurve_data[energy_level] = lightcurve_dataframe[lightcurve_dataframe["EXPOSURE"] != 0]
    #             else:
    #                 print(f"Warning: 'EXPOSURE' column not found in {lightcurve}")
    #         except FileNotFoundError:
    #             print(f"Error: File {lightcurve} not found.")
    #             continue
    #         except Exception as e:
    #             print(f"Error processing {lightcurve}: {e}")
    #             continue
    #     return lightcurve_data

    @staticmethod
    def get_lightcurve_counts(lightcurve_data):
        return int(lightcurve_data["broad"]["COUNTS"].sum())

    @staticmethod
    def create_csv(lightcurve_data):
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
        output_csv = StringIO()
        combined_data.to_csv(output_csv, index=False)
        return output_csv

    def plot(self, lightcurve_data):
        # The type casts are important as the data is returned by CIAO as NumPy data types.
        observation_data = ObservationData(
            average_count_rate=float(round(lightcurve_data["broad"]["COUNT_RATE"].mean(), 3)),
            total_counts=self.get_lightcurve_counts(lightcurve_data),
            total_exposure_time=float(round(lightcurve_data["broad"]["EXPOSURE"].sum(), 3)),
            raw_start_time=int(lightcurve_data["broad"]["TIME"].min()),
        )
        # This data is just so people can view the exact numerical data that was plotted.
        # output_plot_data = StringIO(lightcurve_data["broad"].to_string())
        return LightcurveParseResults(
            observation_header_info=self.get_observation_details(),
            observation_data=observation_data,
            plot_csv_data= self.create_csv(lightcurve_data),
            plot_svg_data=self.create_plot(lightcurve_data, self.binsize),
            postagestamp_png_data=plot_postagestamps(
                self.sky_coords_image, self.detector_coords_image
            ),
        )

    def create_plot(self, lightcurve_data: dict[str, DataFrame], binsize):
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
        plt.savefig(svg_data := StringIO(), bbox_inches="tight")
        plt.close(figure)
        return svg_data

class HrcProcessor(ObservationProcessor):
    """Processes observations produced by the HRC (High Resolution Camera) instrument aboard Chandra."""

    # Define energy level as a constant for one band
    ENERGY_RANGE = "broad=100:10000"  # 0.1 - 10.0 keV

    def adjust_binsize(self, event_list):
        """
        Adjusts the binsize according to the HRC time resolution.

        Parameters:
        event_list (str): Path to the event list file.

        Returns:
        float: Adjusted binsize.
        """
        try:
            time_resolution = float(dmkeypar(infile=event_list, keyword="TIMEDEL", echo=True))
            return max(self.binsize // time_resolution * time_resolution, time_resolution)
        except Exception as e:
            raise ValueError(f"Failed to adjust binsize: {e}")

    def extract_lightcurves(self, event_list):
        """
        Extracts lightcurves from the event list file.

        Parameters:
        event_list (str): Path to the event list file.

        Returns:
        pathlib.Path: Path to the extracted lightcurve file.
        """
        adjusted_binsize = self.adjust_binsize(event_list)
        output_file = Path(f"{event_list}.broad.lc")
        dmextract(
            infile=f"{event_list}[{HrcProcessor.ENERGY_RANGE}][bin time=::{adjusted_binsize}]",
            outfile=str(output_file),
            opt="ltc1",
            clobber="yes",
        )
        return output_file

    def filter_lightcurve_columns(self, lightcurve):
        """
        Filters and extracts specific columns from lightcurve files.

        Parameters:
        lightcurve (pathlib.Path): Lightcurve file path.

        Returns:
        pathlib.Path: Path to the filtered lightcurve file.
        """
        output_file = lightcurve.with_suffix('.ascii')
        dmlist(
            infile=f"{lightcurve}[cols time,count_rate,count_rate_err,counts,exposure,area]",
            opt="data,clean",
            outfile=str(output_file),
        )
        return output_file

    def get_lightcurve_data(self, lightcurve):
        """
        Reads lightcurve data from file and returns as a DataFrame.

        Parameters:
        lightcurve (pathlib.Path): Lightcurve file path.

        Returns:
        pandas.DataFrame: Lightcurve data.
        """
        data = table.Table.read(lightcurve, format="ascii").to_pandas()
        # Filter out zero exposure points
        return data[data["broad"]["EXPOSURE"] != 0]

    def get_lightcurve_counts(self, lightcurve_data):
        """
        Sums the total counts from the lightcurve data.

        Parameters:
        lightcurve_data (pandas.DataFrame): Lightcurve data.

        Returns:
        int: Total counts.
        """
        return int(lightcurve_data["broad"]["COUNTS"].sum())

    @staticmethod
    def create_csv(lightcurve_data):
        """Create CSV with columns for each energy level.

        Parameters
        ----------
        lightcurve_data : dict
            Dictionary of DataFrames containing lightcurve data.

        Returns
        -------
        StringIO
            StringIO object with CSV data.
        """
        combined_data = DataFrame({
            "time": lightcurve_data["broad"]["TIME"],
            "count_rate": lightcurve_data["broad"]["COUNT_RATE"],
            "counts": lightcurve_data["broad"]["COUNTS"],
            "count_error": lightcurve_data["broad"]["COUNT_RATE_ERR"],
            "exposure": lightcurve_data["broad"]["EXPOSURE"],
            "area": lightcurve_data["broad"]["AREA"],
        })
        output_csv = StringIO()
        combined_data.to_csv(output_csv, index=False)
        return output_csv

    def plot(self, lightcurve_data):
        """
        Generates plots for the lightcurve data.

        Parameters:
        lightcurve_data (pandas.DataFrame): Lightcurve data.

        Returns:
        LightcurveParseResults: Results including observation data and plots.
        """
        observation_data = ObservationData(
            average_count_rate=float(lightcurve_data["broad"]["COUNT_RATE"].mean()),
            total_counts=self.get_lightcurve_counts(lightcurve_data),
            total_exposure_time=float(lightcurve_data["broad"]["EXPOSURE"].sum()),
            raw_start_time=int(lightcurve_data["broad"]["TIME"].min()),
        )

        # output_plot_data = StringIO(lightcurve_data.to_string())
        return LightcurveParseResults(
            observation_header_info=self.get_observation_details(),
            observation_data=observation_data,
            plot_csv_data=self.create_csv(lightcurve_data),
            plot_svg_data=self.create_plot(lightcurve_data, self.binsize),
            postagestamp_png_data=plot_postagestamps(
                self.sky_coords_image, self.detector_coords_image
            ),
        )

    @staticmethod
    def create_plot(lightcurve_data: dict[str, DataFrame], binsize):
        """Generate a plt plot to model the lightcurves.

        Parameters
        ----------
        lightcurve_data : dict
            Dictionary of DataFrames containing lightcurve data.
        binsize : float
            Binsize used for the lightcurve.

        Returns
        -------
        StringIO
            StringIO object with SVG data of the generated plot.
        """
        plt.switch_backend("svg")
        initial_time = lightcurve_data["broad"]["TIME"].min()
        zero_shifted_time_ks = (lightcurve_data["broad"]["TIME"] - initial_time) / 1000
        observation_duration = zero_shifted_time_ks.max()

        fig, (count_rate_plot, counts_plot) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), constrained_layout=True)

        # Plot for count rates
        count_rate_plot.errorbar(
            x=zero_shifted_time_ks,
            y=lightcurve_data["broad"]["COUNT_RATE"],
            yerr=lightcurve_data["broad"]["COUNT_RATE_ERR"],
            color="blue",
            marker="o",
            markerfacecolor="black",
            markersize=4,
            ecolor="black",
            markeredgecolor="black",
            capsize=3,
        )
        count_rate_plot.set_xlim([0, observation_duration])
        count_rate_plot.set_title("Broadband Count Rate", fontsize=14)
        count_rate_plot.set_ylabel("Count Rate (counts/s)", fontsize=12)
        count_rate_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
        count_rate_plot.xaxis.set_major_locator(MultipleLocator(5))
        count_rate_plot.xaxis.set_minor_locator(MultipleLocator(1))
        count_rate_plot.tick_params(axis='both', which='major', labelsize=10)

        # Plot for counts
        counts_plot.plot(
            zero_shifted_time_ks,
            lightcurve_data["broad"]["COUNTS"],
            color="purple",
            label="Broadband Counts",
            marker='o',
            markersize=4,
        )
        counts_plot.set_xlim([0, observation_duration])
        counts_plot.set_title("Counts in Broadband", fontsize=14)
        counts_plot.set_xlabel("Time (kiloseconds)", fontsize=12)
        counts_plot.set_ylabel("Counts", fontsize=12)
        counts_plot.grid(True, which='both', linestyle='--', linewidth=0.5)
        counts_plot.xaxis.set_major_locator(MultipleLocator(5))
        counts_plot.xaxis.set_minor_locator(MultipleLocator(1))
        counts_plot.tick_params(axis='both', which='major', labelsize=10)

        fig.suptitle(f"Lightcurve in Broadband (Binsize of {binsize}s)", fontsize=16)
        svg_data = StringIO()
        plt.savefig(svg_data, bbox_inches="tight")
        plt.close(fig)
        return svg_data
