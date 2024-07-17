In difference to the Lightcurves original code, here an additional folder is created (the DipFlare Directory). 
This directory has 3 subdirectories of dips and flare >statistical significance (< p-value cutoff) and all lightcurves with average count rate > cutoff. The subdirectories contain .txt. files whose names include the source and obs-id of the relevant lightcurves. For the dips and flares the exact p-value is reported as well. If there are >1 dips, the longest dip is considered. 

To run code without multiprocessing (somewhat slower but does not suddenly stop working, unlike Lightcurves original):
`python pro.py`
Define parameters at the top of the code (there is no argparse). `pro.py` builds txt files to keep track of already downloaded sources. In the downloading phase, the code can be stopped and restarted + started in multiple terminals with changes to source list ([::-1] to go backwards for example) to emulate multi-processing. Note that if ever the threads try to download the same thing at the same time, the program will break.

Disucssion of analysis of statistical significance:

Without loss of generality, this disucssion also applies to flares.

For a given lightcurve, we want to test the hypothesis “H0: there is no dip, the k points in a row below threshold happened by accident” versus the alternative “HA: there is a dip, the points below threshold were not an accident”. The standard way to do this is to assume that H0 is true and to find the probability (conditional on that H0 is true) that any other randomly collected sample will result in a lightcurve as extreme or even more extreme as our lightcurve (this is called “rejection region” for the test). This probability is called “p-value”. If such probability is small (in statistics, the “small” means below p=5% or p=10%), we can reject H0 and can claim that, most likely (with 95% or 90% confidence for p=5% and 10% respectively) HA is true.

How we define “rejection region”:

If our lightcurve has k points in a row below the dip cutoff, the “as extreme or even more extreme” result would be that a randomly collected sample (conditional on that H0 is true) will have at least k points in a row below the cutoff. In our case, the “randomly collected sample” is just a random simulation of the lightcurve with computed sample mean (0) and standard deviation. This is done with a numpy built in function for random sampling, we repeat 10000 experiments. No seed is provided. 

Why we use all points to compute standard deviation (not just all that are not part of the potential dip):

It is important that we conduct the test assuming H0 is true, this is why we need to use all points in computing sample standard deviation. If we exclude points that are below the dip cutoff while computing standard deviation, the test will not be valid since the random simulation of the lightcurve must be done under the assumption that H0 is true (i.e., that the observed dip is not really a dip but is due to the random variability of the points). 

If we calculated the probability with a distribution whose standard deviation would be based on just the not-potential-dip points, we would be asking the question "Assuming this is a dip, how lucky were we to find it?" instead of "What is the chance this is not a real dip, the k points happened by accident?"

Future work: all observations for Lomb Scargle. 

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