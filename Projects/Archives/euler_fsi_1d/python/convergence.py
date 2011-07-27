#!/usr/bin/python
import sys
import os
from math import log

def Interpolate(x1,t1,x2,t2,t):
    alpha=(t-t1)/(t2-t1)
    x=(1-alpha)*x1+alpha*x2
    return x

def Floatize(list_input):
    return [float(entry) for entry in list_input]

def Intize(list_input):
    return [int(entry) for entry in list_input]

def Read_Position_At_Time(filename,time_input):
    filehandle=open(filename,'r')
    done=False
    (time_last,position_last,velocity_last)=(0,0,0)
    obtained_previous_sample=False
    while(not done):
        (time,position,velocity)=Floatize(filehandle.readline().strip().split(" "))
        if(float(time)<time_input):
            time_last=time
            position_last=position
            velocity_last=velocity
            obtained_previous_sample=True
        else:
            done=True
            if((time>time_input) & obtained_previous_sample):
                print "Interpolating time values between time %f and %f"%(time_last,time)
                position=Interpolate(position_last,time_last,position,time,time_input)
                velocity=Interpolate(velocity_last,time_last,velocity,time,time_input)
            
    print "position=",position,"velocity=",velocity
    return (position,velocity)
    filehandle.close()

def Write_To_Gnuplot_File(filename,dimension,data_list,separator):
    filehandle=open(filename,'w')

    data=["" for i in range(0,dimension)]
    number_of_entries=len(data_list[0])
    for i in range(0,number_of_entries):
        for d in range(dimension):
            data[d]=str(data_list[d][i])
        filehandle.write(separator.join(data)+"\n")

    filehandle.close()

def Convergence(resolutions,directories,time,gnuplot_output_file):
    position_list=["" for i in directories]
    velocity_list=["" for i in directories]
    for dir_index,directory in enumerate(directories):
        filename=directory+"/common/gnuplot_data.dat"
        (position_list[dir_index],velocity_list[dir_index])=Read_Position_At_Time(filename,time)
    
    position_error_magnitudes=[abs(position-position_list[-1]) for position in position_list[:-1]]
    velocity_error_magnitudes=[abs(velocity-velocity_list[-1]) for velocity in velocity_list[:-1]]
     
    log_error_positions=[log(position_error) for position_error in position_error_magnitudes]
    log_error_velocities=[log(velocity_error) for velocity_error in velocity_error_magnitudes]
    log_resolutions=[log(resolution) for resolution in resolutions[:-1]]

    Write_To_Gnuplot_File(gnuplot_output_file,3,[log_resolutions,log_error_positions,log_error_velocities]," ")

    print "positions=",position_list
    print "velocities=",velocity_list
    print "position_error_magnitudes=",position_error_magnitudes
    print "velocity_error_magnitudes=",velocity_error_magnitudes
    print "log_error_positions=",log_error_positions
    print "log_error_velocities=",log_error_velocities
    print "log_resolutions=",log_resolutions

def Dump_Plot_To_PNG(gnuplot_file,png_file):
    gnuplot_dump_file_command=""
    gnuplot_dump_file_command+='\n'+'unset key'
    gnuplot_dump_file_command+='\n'+r'f(x)=m*x+a'
    gnuplot_dump_file_command+='\n'+'fit f(x) \"%s\" using 1:2 via m,a'%gnuplot_file
    gnuplot_dump_file_command+="\n"+'plot \"%s\" using 1:2'%gnuplot_file
    gnuplot_dump_file_command+='\n'+r'replot f(x)'
    gnuplot_dump_file_command+='\n'+'set terminal png'
    gnuplot_dump_file_command+="\n"+'set output \"%s\"'%png_file
    gnuplot_dump_file_command+='\n'+'replot'
    print "gnuplot_dump_file_command=",gnuplot_dump_file_command
    gnuplot_command="gnuplot -persist <<< \'%s\'"%gnuplot_dump_file_command
    print "gnuplot_command=",gnuplot_command
    os.system(gnuplot_command)

def main():
    if(len(sys.argv)<6):
        print "Usage: %s <time> <gnuplot_output_file> <plot_png_file> <directory_pattern> <resolutionlist>"%sys.argv[0]
        sys.exit(1)
    time=float(sys.argv[1])
    gnuplot_output_file=sys.argv[2]
    plot_png_file=sys.argv[3]
    directory_pattern=sys.argv[4]
    resolutions=Intize(sys.argv[5:])

    directories=[directory_pattern%resolution for resolution in resolutions]

    print "******************************************"
    print "time=",time
    print "gnuplot_output_file=",gnuplot_output_file
    print "plot_png_file=",plot_png_file
    print "directory_pattern=",directory_pattern
    print "resolutions=",resolutions
    print "directories=",directories
    print "******************************************"

    Convergence(resolutions,directories,time,gnuplot_output_file)
    Dump_Plot_To_PNG(gnuplot_output_file,plot_png_file)

if __name__=="__main__":
    main()


#./convergence.py .9 Test_2_201_6401_time_pt9.dat Test_2_201_6401_time_pt9.png ../Standard_Tests/Test_2__Resolution_%d_semiimplicit_mass_1.000000/ 201 401 801 1601 3201 6401
#./convergence.py .008 Test_3_201_3201_time_pt008.dat Test_3_201_3201_time_pt008.png ../Standard_Tests/Test_3__Resolution_%d_semiimplicit_mass_6.000000 201 401 801 1601 3201
#./convergence.py 4 Test_5_201_6401_time_4.dat Test_5_201_6401_time_4.png ../Standard_Tests/Test_5__Resolution_%d_semiimplicit_mass_1.000000/ 201 401 801 1601 3201 6401
