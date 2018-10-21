The cross sections in this directory were made with the small script, WIMP_script.m, in this directory.

If you run WIMP_script.m with matlab, it will produce a text file with the recoil energy in keV and a number proportional to the WIMP differential rate.  This information can be turned into a WIMP PDF by normalizing the curve.

Note that this script has several dependencies:
 - "limitcode.git" code, aailable from the git respository on nero.
 - "DetDesign.git", also available from the git repository on nero.
 - make MatLab aware of the code in these repositories, e.g. with the commands in startup.m found in this directory.

The main code called is limitcode.git/highstat/drdq_corr.m.

I used  highstat/rundata/CDMSliteR2/scripts/computeLimit_bestFit_si_CDMSliteR2.m as a guide, which does not set the "haloPar" argument.
