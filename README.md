# Estimating the "nuclear-recoil fano factor" based on yield widths

This is code to estimate the "nuclear-recoil fano factor", or spread in the number of electron-hole pairs produced by a nuclear recoil depositing a fixed amount of energy.  We use the yield bands reported in the 2004 Edelweiss paper.

See the directory `paper_notebooks` for a summary of the work and details about the calculations of the results.

# Python environment

The notebooks use quite a few libraries; you probably don't want to install them piecemeal!  If you're using the Anaconda python distribution you can set up an environment with the below conda commands.  If you're not using Anaconda python, consider switching?

```
conda env create -f nr_fano_env.yml
conda activate nr_fano
```

# Running tests

If you've made changes to the code, you can check that you haven't broken anything by:

```
cd paper_notebooks
py.test
```

Note that the test will take a while, about half an hour.

If you've installed the nr_fano environment then you should be able to run these commands without installing any additional packages.

