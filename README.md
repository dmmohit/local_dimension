# local_dimension
Quantifies the geometry at different regions of a gridded cosmological field

"counter.c" counts the no. of cells with nonzero value and sum of values inside spheres of different radii within the range specified by rmin, rmax and a step size.  
"fit.py" fits the selected quantity with a power law A*R^D to calculate the local dimension at each grid cell.  
The results are analysed using "fit_results.ipynb".

"input.txt" specifies the path of the map file and the input parameters to "counter.c". The map must be a binary file in the following format:  
Grid dimension along x-axis (ni) - int64  
Grid dimension along y-axis (nj) - int64  
Grid dimension along z-axis (nk) - int64  
Grid spacing in Mpc - float32  
Values of the field in an array of size (ni * nj * nk) - float32

Run counter.c using  
$ gcc -fopenmp -o counter counter.c  
$ ./counter  

lmfit module must be installed before running "fit.py". Run it using  
$ python3 fit.py <r_min> <r_max> <step_size> <input_directory (same as output directory of counter)> <output_filename> <max. no. of threads to use>  
All distances must be in grid units.  
