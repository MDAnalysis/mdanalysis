Running the PepTso example
--------------------------

Run the example from this directory::

   mkdir dcd dat
   convolve-mtss-rotamers.py --resid 47 330  --histogramBins 0 80 1  --clashDistance 2.2  \
          --output "dat/peptso-xrd"  --dcdfilename "dcd/peptso-xrd-47-330" \
          peptso.gro 

Compare the output to the reference::

   diff reference/peptso-xrd-47-330.dat dat/peptso-xrd-47-330.dat

You can also look at the PDF of the distance histogram and compare to
reference/peptso-xrd-47-330.pdf.

The above test can be run by executing the script ``test_mtssl.sh``.




  
