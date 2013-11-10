Package to run the WH 2 lepton analysis.

Getting started
---------------

Set the area up using the commands in doc/setup_area.sh

Running the reference cutflow
-----------------------------
First you need to create the necessary list of input files.
If you are not running from the uci machines, then you need to change
the base input path in the script.

```
cd run
./python/makeList.py --also-placeholders
```

Run the code:

```
SusySel \
 -i filelist/Sherpa_CT10_lllnu_WZ.txt \
 -o /tmp/Sherpa_CT10_lllnu_WZ.root \
 -s Sherpa_CT10_lllnu_WZ \
 | tee SusySel.log
```

Making plots
------------

The histograms are defined and filled in SusyPlot.


```
export ${TAG}=Oct_10 # this is some tag to identify today's results
./python/submitJobs.py --susyplot -t ${TAG}

# merge the outputs when the jobs are done
./python/mergeOutput.py -v out/susyplot/

./python/plotHistos.py \
  --input-dir out/susyplot/ \
  --tag ${TAG} \
  --test
```

davide.gerbaudo@gmail.com
2013
