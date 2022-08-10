Author: Adele Hardie

Email: adele.hardie@ed.ac.uk

#### Requirements:
* pyemma
* A number of MD trajectories

## Markov State Models
Now that all of the required seeded MD simulations have run, we can use the trajectory data to construct a Markov State Model (MSM). [PyEMMA](http://emma-project.org/latest/) is a python library for this, and has their own [tutorials](http://emma-project.org/latest/tutorial.html). There is a lot to MSMs that could span a few workshops on its own, but an [example notebook](03_msm_full.ipynb) is available. For now, here is the result of one build using the data obtained through the methods covered previously:

<img src="figures/msm_final.png" width=500>

The model indicates that this particular way of modelling PTP1B with peptide substrate results in catalytically active conformations 20% of the time. Returning to the idea of allosteric inhibition, if a second model, built for a system including an allosteric binder of interest, showed a decrease in active conformation probability, it would suggest that this binder has inhibition potential.