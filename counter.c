#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/stat.h>
#include<time.h>
#include<omp.h>

#define sqr(x) (x*x)

void main()
{
	FILE *fin, *fmap, **fout;	// parameter file, map and set of output files
	char mapfname[256], outpath[256], *outdir, outfname[256];	// filename of map, output directory path and output filename (dummy)
	long ni, nj, nk, nijk;	// grid dimensions (x, y, z grid coordinates), grid volume (grid units ^ 3)
	int rmin, rmax, nr, dr;	// minimum and maximum radii in grid units, no. of steps in r, step size in radius
	float gs, rminmpc, rmaxmpc, rmpc, I0;	// grid spacing; minimum, maximum and current radii (in Mpc); threshold intensity
	float *I;	// intensity grid
	int ind, index, i, j, k, l, ii, jj, kk, m, n;
	int r, r0, count, *sizes, **xvals, **yvals, **zvals;	// r value (grid units); previous r value (grid units);
															// counter for no. of cells within the annulus bounded by r0 and r;
															// array to store the counts; x, y and z coordinates of cells within the annulus;
	int nonzerocount, *nonzero, *centres;	// no. of nonzero cells; indices of nonzero cells; indices of sampled centres
	int i0, j0, k0, dist2;	// x, y, z coordinates of the chosen origin; counter for different r; Cartesian distance squared
	int **Nvals, Ncount, nthreads, ncentres;	// max no. of threads, no. of centres at which local dimension is to be determined
	float **Lvals, Lcount;
	float start_time, end_time;
	time_t t;
	
	start_time = omp_get_wtime();

	if ((fin = fopen("input.txt", "r")) == NULL)
	{
		fprintf(stderr, "Cannot open input file.\n");
		exit(1);
	}
	
	if (fscanf(fin, "%s %d %d %d %f %s %d %d", mapfname, &rmin, &rmax, &nr, &I0, outpath, &nthreads, &ncentres) != 8)
	{
        fprintf(stderr, "Error reading data from input file.\n");
        exit(1);
    }
	fclose(fin);
	
	outdir = (char*) malloc((strlen(outpath)+1) * sizeof(char));	// removing extra characters as output filename has to be appended
	strcpy(outdir, outpath);
		
	if ((fmap = fopen(mapfname, "rb")) == NULL)    
  	{
		fprintf(stderr, "Cannot open %s.\n", mapfname);
		exit(1);
  	}
  	
  	fread(&ni, sizeof(long), 1, fmap);
  	fread(&nj, sizeof(long), 1, fmap);
  	fread(&nk, sizeof(long), 1, fmap);
	fread(&gs, sizeof(float), 1, fmap);
	
	rminmpc = rmin*gs;	// conversion from grid to Mpc units
	rmaxmpc = rmax*gs;
	dr = (int) (rmax - rmin)/nr;
	
	printf("map filename: %s\n", mapfname);
	printf("ni = %ld,\tnj = %ld,\tnk = %ld\n", ni, nj, nk);
	printf("grid spacing = %.2f Mpc\n", gs);
    printf("rmin = %d grid units (%.2f Mpc)\n", rmin, rminmpc);
    printf("rmax = %d grid units (%.2f Mpc)\n", rmax, rmaxmpc);
    printf("nr = %d\n", nr);
    printf("dr = %d\n", dr);
    printf("threshold intensity = %f (map units)\n", I0);
	printf("bright centres to sample = %d\n", ncentres);
    printf("output directory: %s\n", outdir);
    
	nijk = ni*nj*nk;
	I = (float*) malloc(nijk * sizeof(float));
  	
  	fread(I, sizeof(float), nijk, fmap);
	fclose(fmap);

	printf("map data read successfully.\n");
	
	// counting the exact number of steps to be taken in r which need not be the same as the input value
	// useful in case the given nr is not divisible by the radius range
	nr = -1;
	for(r = rmin; r <= rmax; r += dr)
		nr++;
	
	sizes = (int*) malloc((nr+1) * sizeof(int));
	xvals = (int**) malloc((nr+1) * sizeof(int*));
	yvals = (int**) malloc((nr+1) * sizeof(int*));
	zvals = (int**) malloc((nr+1) * sizeof(int*));
	i0 = 0; j0 = 0; k0 = 0;	// setting the origin at (0,0,0)
	r0 = 0;
	
	// storing the cell indices lying in the region bounded by r0 and r centred at the chosen origin
	for(r = rmin, n = 0; r <= rmax; r += dr, n++)
	{
		count = 0;
		// looping within the cube of side 2r containing the sphere of radius r
		for(i = i0-r; i <= i0+r; i++)
			for(j = j0-r; j <= j0+r; j++)
				for(k = k0-r; k <= k0+r; k++)
				{
					ii = i - i0;
					jj = j - j0;
					kk = k - k0;
					dist2 = sqr(ii) + sqr(jj) + sqr(kk);
					// In case of the distance being exactly equal to r, it will be counted for the next r value.
					count += ((dist2 >= sqr(r0)) && (dist2 < sqr(r))) ? 1 : 0;
					//printf("%d %d %d\n", dist2, sqr(r0), sqr(r));
				}
			
		// printf("%d\t%d\n", n, r);
        xvals[n] = (int*) malloc(count * sizeof(int));	// I did not want to allocate any extra memory
		yvals[n] = (int*) malloc(count * sizeof(int));	// so I ran the loop to count and allocate exactly the memory required
        zvals[n] = (int*) malloc(count * sizeof(int));	// and running the loop again to store the coordinates.
        // printf("%d\n", count);
		sizes[n] = count;
        	
        count = 0;
		for(i = i0-r; i <= i0+r; i++)
			for(j = j0-r; j <= j0+r; j++)
				for(k = k0-r; k <= k0+r; k++)
				{
					ii = i - i0;
					jj = j - j0;
					kk = k - k0;
					dist2 = sqr(ii) + sqr(jj) + sqr(kk);
					if((dist2 >= sqr(r0)) && (dist2 < sqr(r)))
					{
						xvals[n][count] = i;
						yvals[n][count] = j;
						zvals[n][count] = k;
						count ++;
					}
				}
        	
		r0 = r;
	}

	printf("cells within each sphere identified.\n");

	nonzerocount = 0;
	for (index = 0; index < nijk; index++)
		nonzerocount += (I[index] > I0) ? 1 : 0;
		
	nonzero = (int*) malloc(nonzerocount * sizeof(int));
	centres = (int*) malloc(ncentres * sizeof(int));
	
	m = 0;
	for (index = 0; index < nijk; index++)	{
		if (I[index] > I0)
			nonzero[m++] = index;
	}

	srand(time(&t));
	for (m = 0; m < ncentres; m++)
		centres[m] = nonzero[rand() % nonzerocount];

	printf("sampled bright centres.\n");
	
    struct stat st;	// checking if the output directory exists
    if (stat(outdir, &st) == -1) // Directory doesn't exist
    {
        if (mkdir(outdir, 0777) == -1)	// checking rwx permissions
        {
            perror("mkdir");
            exit(1);
        }
        printf("output directory created: %s\n", outdir);
    } 
    else
        printf("output directory already exists.\n");
        
	fout = (FILE**) malloc((nr+1) * sizeof(FILE*));
    
	for(r = rmin, n = 0; r <= rmax; r += dr, n++)
	{
		rmpc = r*gs;
		sprintf(outfname, "%sN_values_%.2fMpc", outdir, rmpc);
		fout[n] = fopen(outfname, "wb");
	}

	Nvals = (int**) malloc(ncentres * sizeof(int*));
	Lvals = (float**) malloc(ncentres * sizeof(float*));
	for(l = 0; l < ncentres; l++)	{
		Nvals[l] = (int*) malloc((nr+1) * sizeof(int));
		Lvals[l] = (float*) malloc((nr+1) * sizeof(float));
	}

	#pragma omp parallel for num_threads(nthreads) private(i,j,k,ind,Ncount,Lcount,r,n,m,ii,jj,kk,index)
	for (l = 0; l < ncentres; l++)
	{
		ind = centres[l];
		k = ind % nk;
		ind /= nk;
		j = ind % nj;
		i = ind / nj;
		Ncount = 0;
		Lcount = 0.0;
		for(r = rmin, n = 0; r <= rmax; r += dr, n++)
		{
			Nvals[l][n] = Ncount;
			Lvals[l][n] = Lcount;
			for(m = 0; m < sizes[n]; m++)
			{
				// coordinate transformation (adding grid dimensions to keep the values non-negative)
				ii = xvals[n][m] + i - i0 + ni;
				jj = yvals[n][m] + j - j0 + nj;
				kk = zvals[n][m] + k - k0 + nk;
				// applying periodic boundary conditions
				ii -= ni*(ii/ni);
				jj -= nj*(jj/nj);
				kk -= nk*(kk/nk);
				index = (ii*nj + jj)*nk + kk;
				if(I[index] > I0)
				{
					Nvals[l][n] += 1;
					Lvals[l][n] += I[index];
				}
			}
			Ncount = Nvals[l][n];
			Lcount = Lvals[l][n];
		}
	}

	for(l = 0; l < ncentres; l++)	
	{
		for(n = 0; n < nr+1; n++)
		{
			fwrite(&i, 1, sizeof(int), fout[n]);
			fwrite(&j, 1, sizeof(int), fout[n]);
			fwrite(&k, 1, sizeof(int), fout[n]);
			fwrite(&Nvals[l][n], 1, sizeof(int), fout[n]);
			fwrite(&Lvals[l][n], 1, sizeof(float), fout[n]);
		}
		free(Nvals[l]);
		free(Lvals[l]);
	}
	
	for(r = rmin, n = 0; r <= rmax; r += dr, n++)	{
		fclose(fout[n]);
		free(xvals[n]); free(yvals[n]); free(zvals[n]);
	}

	// for(l = 0; l < ncentres; l++)	{
	// 	free(Nvals[l]);
	// 	free(Lvals[l]);
	// }
	free(xvals); free(yvals); free(zvals);
	free(Nvals); free(Lvals);
	free(sizes); free(nonzero); free(centres);
		
	end_time = omp_get_wtime();
	printf("Executed in %fs\n", end_time - start_time);
}