# Model-Community-Atmospheric-Model3

This repo document my script with NCAR CAM3, which is a spectral AGCM. Comparing to CAM5, it is much smaller and require much less computation resources.

A brief introduction can be found [here](https://github.com/QQFRaphael/Note-Numerical-Methods/blob/master/CESM%E4%BD%BF%E7%94%A8%E6%8C%87%E5%8D%97.md)

go: build and config model

jobcam3: submit the job

plot dir: some example script. 

I usually extract PS, Z3, hyam, hybm, hyai and hybi using CDO and then write code to calculate.

`cdo selname,PS,Z3 infile.nc outfile.nc` 

`cdo selname,hyam,hybm,hyai,hybi infile.nc outfile.nc`