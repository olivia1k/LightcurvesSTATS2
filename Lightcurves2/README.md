In difference to the Lightcurves original code, here an additional folder is created (the DipFlare Directory). 
This directory has 3 subdirectories of dips and flare >statistical significance (< p-value cutoff) and all lightcurves with average count rate > cutoff. The subdirectories contain .txt. files whose names include the source and obs-id of the relevant lightcurves. For the dips and flares the probabilities calculated by both statistical analyses for Normal and Poisson distribution assumptions (see below) and the time-intervals of the dips/flares are reported as well.

To run code without multiprocessing (somewhat slower but does not suddenly stop working, unlike Lightcurves original):
`python pro.py`
Define parameters at the top of the code (there is no argparse). `pro.py` builds txt files to keep track of already downloaded sources. In the downloading phase, the code can be stopped and restarted + started in multiple terminals with changes to source list ([::-1] to go backwards for example) to emulate multi-processing. Note that if ever the threads try to download the same thing at the same time, the program will break.

`pro.py` also has the option to make csv and images only for data that meets the average counts/sec threshold and records flares/dips only for sources-observations if >=1 meets a threshold probability (determined by the Poisson distribution, removing dips/flares).

There are 2 options for statistical significance analysis:

1) Test the hypothesis “H0: there is no dip, k points in a row below dip threshold happened by accident” versus the alternative “HA: there is a dip, the points below threshold were not an accident.” This is done by considering all points in the lightcurve in the calculation of the noise distribution and simulating data based on that. 2 simulations are done: assuming a Normal distribution and assuming Poisson. 

2) Ask "Assuming this is a dip, how lucky were we to find it?" This is done by exlcuding the points below threshold for calculating the noise distribution. Again, Normal and Poisson distributions are considered.

See the `automatic_dip_flare(...)` function in `pro.py` and `lightcurve_processing.py` for more details.

Random seed is provided.

Note that this code does not consider a changing baseline. Many observations of interest are short, and considering error with respect to a moving average removes our ability to consider points early/late in the lightcurves. This code exclude potential dips/flares at the edges of the observations (i.e., which start or end at the first or last positions). 

Also note that if there is a significant dip ("low outlier points"), what would be a flare relative to the total average may no longer a flare after the dip points been removed from consideration. (The reverse is also true.) Std decreases with removing low/high outliers (potential dips/flares), so in this case, misleadingly, the predicted probability of observing an "extreme" value randomly decreases. Thus, what are actually not flares may appear to be very significant by the analysis with extreme point removal.

<br /> 

<br />
Future work: using all observations for Lomb Scargle periodicity analyis. 

Consider calculating probabilites based on "real dip extent" instead of just the theshold value: however, this would take a bit longer. Could consider saving just for dips of interest. <br />

<br />

Instructions from Lightcurves original code: 

CIAO 4.15 must be installed as well as the other project dependencies.

In the correct conda env:
`python main.py`
or
`python main.py --no-gui`

The server for manipulation can be started with:
`python server.py`\
The user will be then prompted for an absolute path to the index file.\
Example: `/Users/mihirpatankar/Documents/Projects/Lightcurves/output/M3-2023-12-08_18:07:31/index.html`

You can also run this tool on a batch of objects.
First you must put your list of objects in `batch_run/objects_list.txt` (one object per line).
Then run `python batch_run/batch_run.py`
The batch progress will be continuously output to `batch_run/current_progress.txt`