# pastisML-tess
Code to produce simulated TESS lightcurves based on PASTIS

## Dependency
The package uses the PASTIS package (not public yet; sorry!)

The PASTIS code is available on request, and the main paper describing it is:

```
@ARTICLE{pastis,
       author = {{D{\'\i}az}, R.~F. and {Almenara}, J.~M. and {Santerne}, A. and
         {Moutou}, C. and {Lethuillier}, A. and {Deleuil}, M.},
        title = "{PASTIS: Bayesian extrasolar planet validation - I. General framework, models, and performance}",
      journal = {\mnras},
     keywords = {methods: statistical, techniques: photometric, techniques: radial velocities, planetary systems, Astrophysics - Earth and Planetary Astrophysics},
         year = "2014",
        month = "Jun",
       volume = {441},
       number = {2},
        pages = {983-1004},
          doi = {10.1093/mnras/stu601},
archivePrefix = {arXiv},
       eprint = {1403.6725},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2014MNRAS.441..983D},
      }
```
    
## Usage
Clone the repo
```
git clone https://github.com/ckm3/pastis-dev.git
```
Enter the folder
```
cd pastis-dev/
```
Create an conda environment with the given environment.yml
```
conda env create -f environment.yml
```
Activate the environment
```
conda activate pastis-env
```
Install pastis
```
cd pastis-core
pip install -e .
```

Also download necessary lib files [https://storage.cuikaiming.com/share/pastis-lib.tar](https://) and put them under the lib folder.
The stellar sample can be found at the shared folder `/storage/tess/armstrong/spoc_ffi`
 