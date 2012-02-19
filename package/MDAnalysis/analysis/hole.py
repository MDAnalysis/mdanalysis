import numpy
import matplotlib
import pylab
import glob
import os

import tempfile

def hole_out(hole_output,run,length,RMSD_array):
    with open (hole_output,"r") as hole:
#count frames
         for i in range(length+1):
#"given end point is never part of the generated list";i will start with 0  http://docs.python.org/tutorial/controlflow.html
             print i
             reading_data=False
             previous_line =" "
# start looking for input for frame with reading_data="False" and previous_line with an empty string
             for line in hole:
                 line = line.strip()
                 if not reading_data:
#analysing output for one frame; looking for header section for HOLE profiles of this frame
                    if line.startswith("Starting calculation for position number"):
                       coord_list=[]
                       radius_list=[]
                       fields = line.split()
                       hole_profile_no = int(fields[5])
#looking for the start of the actual HOLE profile
# RMSD data starts numbering at 0, hole with 1.
                    if previous_line.startswith("cenxyz.cvec") and hole_profile_no==(i+1):
                       reading_data=True
                       print "reading in hole profile for frame "
                       print hole_profile_no
                    else:
                         reading_data=False
                         previous_line=line
                 if reading_data:
                    fields = line.split()
#encountering an empty line when reading in data means that the frame is completely read in.
                    if len(line)==0:
                       number_rows = len(radius_list)
                       coord_1D = numpy.array(coord_list).reshape(number_rows, 1)
                       radius_1D = numpy.array(radius_list).reshape(number_rows, 1)
                       frame_1D=numpy.zeros(number_rows)
                       frame_1D.fill(hole_profile_no)
                       frame_hole_output= numpy.column_stack((frame_1D,coord_1D,radius_1D))
                       frame_hole_output=frame_hole_output.astype(float)
# save a profile for each frame, in principle one doesn`t need to save these profiles; they are saved anyway once stacked together
# for now useful for testing
# convert str to float, can use savetxt
                       frame_hole_txt = run+"_"+str(hole_profile_no)+".txt"
                       numpy.savetxt(frame_hole_txt,frame_hole_output)
                      # a tmp folder for each trajectory
                       if not os.path.exists(run):
                          os.makedirs(run)
                       os.rename(frame_hole_txt,str(run)+"/"+frame_hole_txt)
                       print "finished with frame "+str(hole_profile_no)+" file in "+run+"/"+frame_hole_txt
# add RMSD; this step could omitted when transfering the script to MDAnalysis
                       frame_hole_RMSD_array=master_ar(run_array,hole_profile_no,frame_hole_output)
                       if (i+1)==1:
                          run_hole_array =numpy.array(frame_hole_RMSD_array)
                       else:
                            run_hole_array=numpy.vstack((run_hole_array,frame_hole_RMSD_array))
                            if (i+1)==length:
                               print "finished with run :"+run
                               return run_hole_array
                           #return finished array
                       break #should break out of for loop. I don`t want read the rest of the file
                    else:
                         coord,radius=fields[0],fields[1]
                         coord_list.append(coord)
                         radius_list.append(radius)




def master_ar(run_arrray,frame_number,frame_hole_output): # this "master array is simpler than in a previous implementation
    print frame_number-1
    print " (frame-1 for RMSD)"
    RMSD = run_array[frame_number-1,3]
    print RMSD
    RMSD_column = numpy.zeros((len(frame_hole_output)),dtype=float)
    RMSD_column.fill(RMSD)
    frame_hole_RMSD_array=numpy.hstack((frame_hole_output,RMSD_column[:,numpy.newaxis])) #check
    #http://stackoverflow.com/questions/4158388/numpy-concatenating-multidimensional-and-unidimensional-arrays
    return frame_hole_RMSD_array

path_RMSD="/sansom/hfs0/trin1736/DIMS_storage/GBSWmemb/DIMS_imp_80_p_open_h_f_4_f/RMSD/"
path_hole="/sansom/hfs0/trin1736/DIMS_storage/GBSWmemb/DIMS_imp_80_p_open_h_f_4_f/out/"

for infile in glob.glob(os.path.join(path_RMSD, '*_25-10-11.output.npy') ):
    run_array = numpy.load(infile)
    run_name=run_array[0,0]
    hole_output= path_hole+"/hole_output_30-10-11_"+run_name+".out"
    length= run_array[-1,2]+1 # RMSD data starts numbering at 0, hole with 1. maybe have a test function that just counts frames in hole ouput files
    hole_RMSD_array=hole_out(hole_output, run_name,length,run_array)
    #add RMSD data
    numpy.savetxt("RMSD_hole_"+run_name+".txt",hole_RMSD_array)

