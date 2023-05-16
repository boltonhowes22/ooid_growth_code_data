# ooid_growth_code_data
data and code for ooid growth manuscript


To remake figures within the manuscript

Figure 1: Data comes from binned_xDD_ooids.csv
Figure 3a: run sphericity_plots.m
Figure 3b: figure from Geyman et al. 2022
Figure 3c: Run sadler_comp.m
            - uses some helping functions
            - dating.csv is just the age portion of the sisalv dataset. See that publication for the rest of the dataset. 
Figure 6: Run ooid_growth.m
            - uses some additional functions
            -  most important of which are fitOoids.m (contains all the spherical harmonics). Requires the slepian alpha toolbox which can be downloaded here: https://github.com/csdms-contrib/slepian_alpha
Figure 7: alk_LG_modern.m
Figure 8: ooid_potential.m
Figure 9: resampling.ipynb and the files within the size_resampling directory
