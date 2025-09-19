# pastis-dev

This code simulates TESS planetary and false positive scenarios.

It combines [PASTIS](https://github.com/exord/pastis) and [pastisML-tess](https://github.com/exord/pastisML-tess) into a standalone repository to simulate TESS light curves for planets and astrophysical false positives.

Compared to the original PASTIS and pastisML-tess, we simplified the code to speed up simulations and modified some distributions of physical parameters such as distance, eccentricity, mass-radius, etc. Simulated planetary and false positive signals are used in our series of papers. 

If you find this code useful, please cite the original PASTIS paper and include the link to this GitHub repository.

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
Create an conda environment with the given environment.yml, and if you are going to install other packages with conda, make sure you are using conda-forge channel.
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

Also download necessary lib files [https://storage.cuikaiming.com/share/pastis-lib.tar]([https://](https://storage.cuikaiming.com/share/pastis-lib.tar)) and put them under the lib folder.
The stellar sample can be found at the shared folder `/storage/tess/armstrong/spoc_ffi`
 