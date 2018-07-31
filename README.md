*Last edit 7/30/18*

Analysis consists of using three programs: checkrescode.py, spectraplot.py, and residplot.py
--------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------




checkrescode.py information
----------------------------
Usage:   
./checkrescode.py PSFfilename FRAMEfilename

Output:

- OUTfilename.csv  
A csv file with a header of Row Number, Fiber ID, Legendre Polynomial Coefficient 1 Difference, Legendre Polynomial Coefficient 2 Difference, Legendre Polynomial Coefficient 2 Difference. Row Number is to be used as input for the spectraplot.py script. It represents the row of the fiber in the original WSIGMA array. The Legendre Polynomial Coefficient Differences are (Coefficient value - median coefficient across the 500 fibers). It is used as a measure of how 'bad' a coefficient is. The larger the value, the worse the coefficient. This file is used to determine which fibers to plot. Check the row number of the fiber you want to plot!

- Medians  
The WSIGMA array median values used are printed to stdout. Output is in an array of 3 values corresponding to [p0, p1, p2] coefficients.

- Standard deviation values  
The standard deviation values used are printed to stdout. Output is in an array of 3 values corresponding to [p0, p1, p2] coefficients.

- Number of 'bad' Coefficients  
The number of 'bad' coefficients for each parameter are printed to stdout. Output is in an array of 3 values corresponding to [p0, p1, p2] coefficients. A coefficient is 'bad' if it is >2RMS.

- Rows (Fibers) with bad parameters  
If a row (fiber) of the WSIGMA array contains a bad coefficient, it is sorted and the row number is printed to stdout. Three lists are printed, 'The rows with a bad param0', 'The rows with a bad param1', 'The rows with a bad param2'.

- Ideal Fiber Rows  
A list with all row numbers (fibers) of the WSIGMA array where all three coefficients are <1RMS from the median. List is printed to stdout

- Fiber IDs for 2 bad parameters  
A list with all Fibers that have 2 'bad' coefficients. Printed to stdout. These are the FIBER IDs, NOT the WSIGMA array row number.

- Fiber IDs for 3 bad parameters  
A list with all Fibers that have 3 'bad' coefficients. Printed to stdout. These are the FIBER IDs, NOT the WSIGMA array row number.


Example:

	./checkrescode.py psf-r0-00003891.fits frame-r0-00003891.fits


spectraplot.py information
----------------------------
Usage:  
./spectraplot.py FRAMEfilename idealfiberrownumber nearbyfiberrownumber worstfiberrownumber  

Output:  
A plot of the overlaid sprectrums of the fibers identified by their row in the WSIGMA array (output of ./checkrescode.py).

Input:  
- idealfiberrownumber  
The row number of the ideal fiber for which you want to plot the spectrum. Chosen close to row number of worst fiber.

- nearbyfiberrownumber  
The row number of a fiber close to the chosen ideal which has only ONE bad parameter.

- worstfiberrownumber  
The row number of the worst fiber as identified in the csv output file from checkrescode.py  

Example:

	./spectraplot.py frame-r0-00003891.fits 429 430 432


residplot.py information
-------------------------
Usage:  
./residplot.py framefile idealfiberrow worstfiberrow  

Output:  
-outfilename.eps  
A plot of worst fiber spectrum over ideal fiber spectrum with residual subplot.

Example:

	./residplot.py frame-r0-00003891.fits 429 432


