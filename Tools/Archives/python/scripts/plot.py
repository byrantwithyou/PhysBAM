#!/usr/bin/python

import sys
import os
import readline
#from physbam import *
from subprocess import *
import Gnuplot, Gnuplot.funcutils
from optparse import OptionParser
import re

class PLOT:
    def Get_Directory_Base_Name(self,directory_path):
        directory_name=os.path.dirname(directory_path+"/") # To remove any trailing slashes. An extra '/' is passed for the case when 'directory' didn't have any '/' initially (else dirname gets confused)
        directory_base_name=os.path.basename(directory_name)
        return directory_base_name

    def Remove_Preceding_Dashes(self,string_with_dash):
        for i in range(0,len(string_with_dash)):
            if(string_with_dash[i]=='-'): pass
            else: break
        return string_with_dash[i:]

    def Set_Max_Min(self,value,min_value=None,max_value=None):
        return_list=[]
        if(min_value!=None):
            if(value<min_value): min_value=value
            return_list.append(min_value)
        if(max_value!=None):
            if(value>max_value): max_value=value
            return_list.append(max_value)
        return return_list

    def Interpolate(self,x1,t1,x2,t2,t):
        alpha=(t-t1)/(t2-t1)
        x=(1-alpha)*x1+alpha*x2
        return x

    #def PhysBAM_Array_1D_Into_List(self,array_physbam):
    #    """Takes a physbam array and returns domain indices and a python list containing array data"""

    #    domain_indices_box=array_physbam.Domain_Indices()

    #    domain_indices=range(domain_indices_box.min_corner.x,domain_indices_box.max_corner.x+1)

    #    list_array=[0 for i in domain_indices]
    #    array_index=Vi1()
    #    list_index=0
    #    for i in domain_indices:
    #        array_index.x=i
    #        list_array[list_index]=array_physbam.__getitem__(array_index)
    #        list_index+=1
    #    return (domain_indices,list_array)

    #def Get_Node_Locations_from_PhysBAM_Grid_1D(self,grid_physbam,node_indices):
    #    """Reads a physbam grid and returns a python list of the node positions"""

    #    node_locations=[0 for i in node_indices]
    #    list_index=0
    #    for i in node_indices:
    #        node_locations[list_index]=grid_physbam.X(i).x
    #        list_index+=1
    #    return node_locations

    #def Read_Grid_And_Array_Data(self,directory,variable_name,frame_number,stream_type):
    #    grid_file=directory+"/grid.%d.gz"%frame_number
    #    array_file=directory+"/"+variable_name+".%d.gz"%frame_number

    #    if(stream_type=='float'): grid_physbam=GRID_Vf1()
    #    else: grid_physbam=GRID_Vd1()
    #    Read_From_File(stream_type,grid_file,grid_physbam)

    #    if(stream_type=='float'): array_physbam=ARRAYS_1D_f()
    #    else: array_physbam=ARRAYS_1D_d()
    #    Read_From_File(stream_type,array_file,array_physbam)

    #    (domain_indices,list_array)=self.PhysBAM_Array_1D_Into_List(array_physbam)
    #    node_locations=self.Get_Node_Locations_from_PhysBAM_Grid_1D(grid_physbam,domain_indices)

    #    return (node_locations,list_array)


    def Read_Rigid_Body_Simplicial_Object_Particle_Positions(self,directory,rigid_body_id,frame_number,stream_type):
        tmp_output_file_name="/tmp/plot_tmp_files/rigid_body_simplicial_positions.txt"
        physbam_rigid_body_read_tool=os.environ["PHYSBAM"]+"/Tools/rigid_body_simplicial_particles_dump/rigid_body_simplicial_particles_dump_nocona"
        physbam_rigid_body_read_cmd="%s %s -o %s -rbid %d -frame %d -d %d"%(physbam_rigid_body_read_tool,directory,tmp_output_file_name,rigid_body_id,frame_number,self.dimension)
        if(stream_type=="double"):
            physbam_rigid_body_read_cmd+=" -double"
        os.system(physbam_rigid_body_read_cmd)

        simplicial_positions_file=open(tmp_output_file_name,'r')
        number_of_particles_line=simplicial_positions_file.readline()
        (temp_dump_string,number_of_particles)=number_of_particles_line.split("=")
        number_of_particles=int(number_of_particles)
        particle_position_list=[0 for i in range(0,number_of_particles)]
        for i in range(0,number_of_particles):
            particle_position_line=simplicial_positions_file.readline()
            (temp_dump_string,particle_position)=particle_position_line.split("=")
            particle_position=particle_position.strip()
            particle_position=particle_position.strip("[]")
            particle_position_list[i]=float(particle_position)

        simplicial_positions_file.close()
        os.system("rm -f %s"%tmp_output_file_name)
        return particle_position_list

    def Read_Deformable_Object_Particle_Positions(self,directory,frame_number,stream_type):
        tmp_output_file_name="/tmp/plot_tmp_files/deformable_body_simplicial_positions.txt"
        physbam_deformable_body_read_tool=os.environ["PHYSBAM"]+"/Tools/particle_dump/particle_dump_nocona"
        physbam_deformable_body_read_cmd="%s -d %d %s/%d/deformable_object_particles.gz"%(physbam_deformable_body_read_tool,self.dimension,directory,frame_number)
        if(stream_type=="double"):
            physbam_deformable_body_read_cmd+=" -double"
        physbam_deformable_body_read_cmd+=" |grep 'X ='|sed s,'X = \[',,g|sed s,],,g > %s"%tmp_output_file_name
        os.system(physbam_deformable_body_read_cmd)

        # Assuming just one deformable body
        simplicial_positions_file=open(tmp_output_file_name,'r')
        particle_position_first=float(simplicial_positions_file.readline().strip())
        particle_position_second=float(simplicial_positions_file.readline().strip())
        simplicial_positions_file.close()
        os.system("rm -f %s"%tmp_output_file_name)
        return (particle_position_first,particle_position_second)

    def Remove_Data_Inside_Simplicial_Object(self,particle_position_list,grid_array_data):
        (node_locations,list_array)=grid_array_data
        node_locations_new=[]
        list_array_new=[]
        for index,node_location in enumerate(node_locations):
            if(not ((node_location>=particle_position_list[0]) & (node_location<=particle_position_list[1]))):
                node_locations_new.append(node_location)
                list_array_new.append(list_array[index])

        return (node_locations_new,list_array_new)

    def Zero_Out_Data_Inside_Simplicial_Object(self,particle_position_list,grid_array_data):
        (node_locations,list_array)=grid_array_data
        for index,node_location in enumerate(node_locations):
            if((node_location>=particle_position_list[0]) & (node_location<=particle_position_list[1])):
                list_array[index]=0

        return (node_locations,list_array)

    def Linearly_Interpolate_Data_Inside_Simplicial_Object(self,particle_position_list,grid_array_data):
        (node_locations,list_array)=grid_array_data
        # finds fluid cell index just before and just after the rigid body
        for index,node_location in enumerate(node_locations):
            if(node_location<particle_position_list[0]):
                start_fluid_index=index
            if(node_location>particle_position_list[1]):
                end_fluid_index=index
                break

        # interpolate data
        start_location=node_locations[start_fluid_index]
        start_data_value=list_array[start_fluid_index]
        end_location=node_locations[end_fluid_index]
        end_data_value=list_array[end_fluid_index]
        for index,node_location in enumerate(node_locations):
            if((node_location>=particle_position_list[0]) & (node_location<=particle_position_list[1])):
                list_array[index]=self.Interpolate(start_data_value,start_location,end_data_value,end_location,node_location)

        return (node_locations,list_array)

    def Get_Rigid_Body_Gnuplot_Data_MinMax_Bounds(self,particle_position_list,rigid_thickness):
        tmp_rigid_rectable_file_name="/tmp/plot_tmp_files/tmp_rigid_rectable.dat"
        rigid_body_width=particle_position_list[1]-particle_position_list[0]
        self.gnuplot_plotter("set parametric")
        self.gnuplot_plotter("f(x)=(x<=1)?%f+x*%f:%f+(1-x)*%f"%(particle_position_list[0],rigid_body_width,particle_position_list[1],rigid_body_width))
        self.gnuplot_plotter("g(y)=(y<=1)?%f:%f"%(-rigid_thickness,rigid_thickness))
        self.gnuplot_plotter("set terminal table")
        self.gnuplot_plotter("set output \"%s\""%tmp_rigid_rectable_file_name)
        self.gnuplot_plotter("plot [0:2] f(t),g(t) with filledcurve")
        self.gnuplot_plotter("set terminal x11")
        gnuplot_data=Gnuplot.File(tmp_rigid_rectable_file_name,title="",with="filledcurve lt 3")
        return (gnuplot_data,-rigid_thickness,rigid_thickness)

    def Read_Grid_And_Array_Data_From_Gnuplot_File(self,gnuplot_file_name):
        # assumes that file writes data in x,data format
        gnuplot_file=open(gnuplot_file_name,'r')
        node_locations=[]
        list_array=[]
        for line in gnuplot_file:
            (node_location,value)=line.split()
            node_locations.append(float(node_location))
            list_array.append(float(value))

        gnuplot_file.close()

        return (node_locations,list_array)

    def Read_Grid_And_Array_Data(self,directory,variable_name,frame_number,stream_type):
        tmp_output_directory="/tmp/plot_tmp_files/%s"%self.Get_Directory_Base_Name(directory)
        physbam_to_gnuplot_converter=os.environ["PHYSBAM"]+"/Tools/physbam2gnuplot/physbam2gnuplot_nocona"
        physbam_to_gnuplot_cmd="%s %s -o %s -v %s -start_frame %d -last_frame %d -dimension %d"%(physbam_to_gnuplot_converter,directory,tmp_output_directory,variable_name,frame_number,frame_number,self.dimension)
        if(stream_type=="double"):
            physbam_to_gnuplot_cmd+=" -double"
        os.system(physbam_to_gnuplot_cmd)

        gnuplot_file_name=tmp_output_directory+"/%s.%d"%(variable_name,frame_number)
        grid_array_data=self.Read_Grid_And_Array_Data_From_Gnuplot_File(gnuplot_file_name)

        #os.system('rm -f %s'%gnuplot_file_name)
        return grid_array_data

    def Get_Plot_Title(self,directory,verbose_plot):
        directory_base_name=self.Get_Directory_Base_Name(directory)
        if(verbose_plot==1):
            if("explicit" in directory_base_name):
                plot_title="Explicit"
            elif("semiimplicit" in directory_base_name):
                plot_title="Semiimplicit"
            else:
                plot_title="unknown"
        elif(verbose_plot==2):
            plot_title=directory_base_name
        elif(verbose_plot==3):
            plot_title="%d"%(len(node_locations)-6) # Setting to the resolution of the grid. Subtracting 6 to remove the 6 ghost nodes.
        else:
            plot_title=""

        return plot_title

    def Get_Gnuplot_Data_And_MinMax_Bounds(self,directory,variable_name,frame_number,plot_style,plot_width,stream_type,verbose_plot):
        grid_array_data=self.Read_Grid_And_Array_Data(directory,variable_name,frame_number,stream_type)
        if(self.display_rbody):
            particle_position_list=self.Read_Rigid_Body_Simplicial_Object_Particle_Positions(directory,1,frame_number,stream_type)
            grid_array_data=self.Linearly_Interpolate_Data_Inside_Simplicial_Object(particle_position_list,grid_array_data)
        if(self.display_deformable_body):
            (particle_position_fist,particle_position_second)=self.Read_Deformable_Object_Particle_Positions(directory,frame_number,stream_type)
            grid_array_data=self.Remove_Data_Inside_Simplicial_Object([particle_position_fist,particle_position_second],grid_array_data)

        (node_locations,list_array)=grid_array_data

        min_array_value=min(list_array)
        max_array_value=max(list_array)

        plot_title=self.Get_Plot_Title(directory,verbose_plot)
        gnuplot_data=Gnuplot.Data(node_locations,list_array,title=plot_title,with="%s lw %d"%(plot_style,plot_width))
        if(not(self.set_input_minmax_xrange)): # If not given as input set the xrange here
            self.minxrange=node_locations[0]
            self.maxxrange=node_locations[-1]
        return (gnuplot_data,min_array_value,max_array_value)
    
    def Get_Gnuplot_Data_And_MinMax_Bounds_For_All_Directories(self,directories,variable_name,frame_number,plot_styles,plot_widths,stream_type,verbose_plot):
        """ This function is called to get the data which is plotted. Overwrite this function in an inherited class to change the plots."""
        gnuplot_data_list=[]
        for dir_index,directory in enumerate(directories):
            (gnuplot_data,min_array_value,max_array_value)=self.Get_Gnuplot_Data_And_MinMax_Bounds(directory,variable_name,frame_number,plot_styles[dir_index],plot_widths[dir_index],stream_type,verbose_plot)
            gnuplot_data_list.append(gnuplot_data)

            if(dir_index==0):
                min_plot_value=min_array_value
                max_plot_value=max_array_value
            else:
                min_plot_value,=self.Set_Max_Min(value=min_array_value,min_value=min_plot_value)
                max_plot_value,=self.Set_Max_Min(value=max_array_value,max_value=max_plot_value)     

        if(self.display_rbody or self.display_deformable_body):
            for directory in directories:
                if(self.display_rbody):
                    particle_position_list=self.Read_Rigid_Body_Simplicial_Object_Particle_Positions(directory,1,frame_number,stream_type)
                else:
                    (particle_position_fist,particle_position_second)=self.Read_Deformable_Object_Particle_Positions(directory,frame_number,stream_type)
                    particle_position_list=[particle_position_fist,particle_position_second]

                if(min_plot_value<0):
                    plot_height=max_plot_value-min_plot_value
                else:
                    plot_height=max_plot_value

                (gnuplot_data_rigid,min_array_value,max_array_value)=self.Get_Rigid_Body_Gnuplot_Data_MinMax_Bounds(particle_position_list,plot_height*self.rigid_relative_thickness)
                gnuplot_data_list.append(gnuplot_data_rigid)

                min_plot_value,=self.Set_Max_Min(value=min_array_value,min_value=min_plot_value)
                max_plot_value,=self.Set_Max_Min(value=max_array_value,max_value=max_plot_value)     


        return (gnuplot_data_list,min_plot_value,max_plot_value)

    def Hardcoded_Directories(self,base_directory_path,test_number):
        #directories_format="Test_%d__Resolution_400_eno_order-2_rk_order-3_implicit_rk_semiimplicit,Test_%d__Resolution_800_eno_order-2_rk_order-3_implicit_rk_semiimplicit,Test_%d__Resolution_1600_eno_order-2_rk_order-3_implicit_rk_semiimplicit,Test_%d__Resolution_400_eno_order-2_rk_order-3_implicit_rk_density_weighted_ENO_semiimplicit,Test_%d__Resolution_800_eno_order-2_rk_order-3_implicit_rk_density_weighted_ENO_semiimplicit,Test_%d__Resolution_1600_eno_order-2_rk_order-3_implicit_rk_density_weighted_ENO_semiimplicit,Test_%d__Resolution_400_eno_order-2_rk_order-3_implicit_rk_velocity_weighted_ENO_semiimplicit,Test_%d__Resolution_800_eno_order-2_rk_order-3_implicit_rk_velocity_weighted_ENO_semiimplicit,Test_%d__Resolution_1600_eno_order-2_rk_order-3_implicit_rk_velocity_weighted_ENO_semiimplicit,Test_%d__Resolution_400_eno_order-2_rk_order-3_semiimplicit,Test_%d__Resolution_800_eno_order-2_rk_order-3_semiimplicit,Test_%d__Resolution_1600_eno_order-2_rk_order-3_semiimplicit,Test_%d__Resolution_400_eno_order-2_rk_order-3_density_weighted_ENO_semiimplicit,Test_%d__Resolution_800_eno_order-2_rk_order-3_density_weighted_ENO_semiimplicit,Test_%d__Resolution_1600_eno_order-2_rk_order-3_density_weighted_ENO_semiimplicit,Test_%d__Resolution_400_eno_order-2_rk_order-3_velocity_weighted_ENO_semiimplicit,Test_%d__Resolution_800_eno_order-2_rk_order-3_velocity_weighted_ENO_semiimplicit,Test_%d__Resolution_1600_eno_order-2_rk_order-3_velocity_weighted_ENO_semiimplicit,Test_%d__Resolution_400_eno_order-2_rk_order-3_explicit,Test_%d__Resolution_800_eno_order-2_rk_order-3_explicit,Test_%d__Resolution_1600_eno_order-2_rk_order-3_explicit,Test_%d__Resolution_400_eno_order-2_rk_order-3_density_weighted_ENO_explicit,Test_%d__Resolution_800_eno_order-2_rk_order-3_density_weighted_ENO_explicit,Test_%d__Resolution_1600_eno_order-2_rk_order-3_density_weighted_ENO_explicit,Test_%d__Resolution_400_eno_order-2_rk_order-3_velocity_weighted_ENO_explicit,Test_%d__Resolution_800_eno_order-2_rk_order-3_velocity_weighted_ENO_explicit,Test_%d__Resolution_1600_eno_order-2_rk_order-3_velocity_weighted_ENO_explicit"
        directories_format="Test_%d__Resolution_401_eno_order-2_rk_order-3_implicit_rk_density_weighted_ENO_semiimplicit,Test_%d__Resolution_401_eno_order-2_rk_order-3_explicit"

        return directories_format.replace("Test_%d","%s/Test_%d"%(base_directory_path,test_number))

    def Setup_Parser(self,usage):
        self.parser=OptionParser(usage)
        self.parser.add_option("-v","--variable",action="store",type="string",dest="variable_name",default="density",help="name of the variable to be plotted. Default=%default")
        self.parser.add_option("--rigid",action="store_true",dest="display_rbody",help="display rigid body")
        self.parser.add_option("--deformable",action="store_true",dest="display_deformable_body",help="display deformable body")
        self.parser.add_option("-f","--frame",action="store",type="int",dest="frame_number",default=0,help="Frame number. Default=%default")
        self.parser.add_option("--rig_rel_thick",action="store",type="int",dest="rigid_relative_thickness",default=.01,help="Rigid relative thickness. Default=%default")
        self.parser.add_option("-s","--style",action="store",type="string",dest="plot_style",default="lines",help="style for plotting. Default=%default")
        self.parser.add_option("--fontsize",action="store",type="int",dest="fontsize",default=40,help="font size for writing into file")
        self.parser.add_option("--pw",action="store",type="int",dest="plot_width",default=1,help="style width. Default=%default")
        self.parser.add_option("--verbosity",action="store",type="int",dest="verbose_plot",default=0,help="verbosity level for plot. Default=%default")
        self.parser.add_option("--xlabel",action="store_true",dest="xlabel",help="display xlabel")
        self.parser.add_option("--ylabel",action="store_true",dest="ylabel",help="display ylabel")
        self.parser.add_option("--minxrange",action="store",type="float",dest="minxrange",help="set minimum x range. If set, maximum x range will default to 1")
        self.parser.add_option("--maxxrange",action="store",type="float",dest="maxxrange",help="set maximum x range. If set, minimum x range will default to 0")
        self.parser.add_option("-c","--cmd",action="store",type="string",dest="command",default="",help="Command")
        self.parser.add_option("--path",action="store",type="string",dest="base_directory_path",default=".",help="Base directory path for using hardcoded directories and writing output files. Default=%default")
        self.parser.add_option("-n",action="store",type="int",dest="test_number",default=1,help="Test number for using hardcoded directories. Default=%default")
        self.parser.add_option("--float",action="store_const",const="float",dest="stream_type",help="set stream type to float. Default value is float.")
        self.parser.add_option("--double",action="store_const",const="double",dest="stream_type",help="set stream type to double. Default value is float.")

    def Parse_Command_Line(self,argument_list):
        (options,args)=self.parser.parse_args(argument_list)

        self.variable_name=options.variable_name
        self.display_rbody=options.display_rbody
        self.display_deformable_body=options.display_deformable_body
        self.rigid_relative_thickness=options.rigid_relative_thickness
        self.plot_style=options.plot_style
        self.plot_width=options.plot_width
        self.verbose_plot=options.verbose_plot
        self.xlabel=options.xlabel
        self.ylabel=options.ylabel
        self.frame_number=options.frame_number
        self.fontsize=options.fontsize
        self.base_directory_path=options.base_directory_path
        self.test_number=options.test_number
        self.command=options.command
        self.minxrange=0
        self.maxxrange=1
        if((options.minxrange!=None) | (options.maxxrange!=None)): #Set default values
            self.set_input_minmax_xrange=True
        else:
            self.set_input_minmax_xrange=False
        if(options.minxrange!=None):
            self.minxrange=options.minxrange
        if(options.maxxrange!=None):
            self.maxxrange=options.maxxrange

        if(options.stream_type):
            self.stream_type=options.stream_type
        else:
            print "stream_type not specified. Using float"
            self.stream_type='float'

        if(len(args)):
            self.directories=args
        else:
            self.directories=self.Hardcoded_Directories(self.base_directory_path,self.test_number).split(",")

        return (options,args)

    def Print_Arguments(self):
        print "---------------------------------"
        print "variable_name =",self.variable_name
        print "frame_number =",self.frame_number
        print "fontsize=",self.fontsize
        print "base_directory_path =",self.base_directory_path
        print "test_number =",self.test_number
        print "directories =",self.directories
        print "stream_type =",self.stream_type
        print "---------------------------------"

    def __init__(self,argument_list):
        #get parameters from argument list
        usage="usage: %prog [options] <directory 1> <directory 2> ... <directory n>"
        self.Setup_Parser(usage)
        self.Parse_Command_Line(argument_list)
        self.Print_Arguments()
        
        # set up plotter
        self.gnuplot_plotter=Gnuplot.Gnuplot(debug=1)
        if(self.xlabel): self.gnuplot_plotter.xlabel('x')
        if(self.ylabel): self.gnuplot_plotter.ylabel(self.variable_name)

        # set up regexes for parsing input data
        self.shortoption_re=re.compile(r'(?P<sopt>\b[a-zA-Z])\b')
        self.longoption_re=re.compile(r'(?P<lopt>\b[a-zA-Z]{2,})\b')

        #dimension
        self.dimension=1

        # Misc params
        self.padding_factor=.05
        self.data_range_min=1e13
        self.data_range_max=-1e13

        # set up parser for interactive input
        self.parser_interactive=OptionParser("You can ignore the - and -- in the following options")
        self.parser_interactive.add_option("-d","--dir",action="store",type="int",dest="selected_directory_index",help="select a directory")
        self.parser_interactive.add_option("-v","--var",action="store",type="string",dest="variable_name",help="change variable to be plotted")
        self.parser_interactive.add_option("-g","--goto",action="store",type="int",dest="frame_number",help="goto frame number")
        self.parser_interactive.add_option("-r","--reset",action="store_const",const=0,dest="frame_number",help="goto frame 0")
        self.parser_interactive.add_option("-s","--step",action="store_const",const=1,dest="increment",help="goto next frame")
        self.parser_interactive.add_option("-b","--backstep",action="store_const",const=-1,dest="increment",help="goto previous frame")
        self.parser_interactive.add_option("-j","--jump",action="store",type="int",dest="increment",help="jump a difference")
        self.parser_interactive.add_option("-p","--play",action="store_true",dest="play",help="play till maximum frame plotted till now")
        self.parser_interactive.add_option("-q","--quit",action="store_true",dest="quit",help="quit plotter")
        self.parser_interactive.add_option("--fontsize",action="store",type="int",dest="fontsize",help="font size for writing into file")
        self.parser_interactive.add_option("-w","--write",action="store_true",dest="write",help="write plot to file")
        self.parser_interactive.add_option("--wf",action="store",type="string",dest="filename_prefix_for_write",help="write plot to specified file")
        self.parser_interactive.add_option("--ls",action="store_true",dest="ls",help="list directories being plotted")
        self.parser_interactive.add_option("--st",action="store",type="string",dest="plot_style",help="set style of plot for the chosen directory")
        self.parser_interactive.add_option("--pw",action="store",type="int",dest="plot_width",help="style width.")
    
    def Set_Title(self):
        if(self.verbose_plot):
            self.gnuplot_plotter.title("Frame %d"%self.frame_number)

    def Set_Grid_Axis_Labels(self):
        if(self.verbose_plot):
            self.gnuplot_plotter.xlabel("x")

    def Set_Data_Axis_Labels(self):
        if(self.verbose_plot):
            self.gnuplot_plotter.ylabel(self.variable_name)

    def Set_Grid_Axis_Range(self):
        xrange_set_string="set xrange [%f:%f]"%(self.minxrange,self.maxxrange)
        self.gnuplot_plotter(xrange_set_string)

    def Set_Data_Axis_Range(self):
        data_range_paddding=(self.data_range_max-self.data_range_min)*self.padding_factor
        if(data_range_paddding!=0):
            data_range_set_string="set yrange [%f:%f]"%(self.data_range_min-data_range_paddding,self.data_range_max+data_range_paddding)
            self.gnuplot_plotter(data_range_set_string)

    def Plot(self,gnuplot_data_list):
        self.gnuplot_plotter.plot(*gnuplot_data_list)

    def Load_Additional_Options(self,options):
        """ This function will be called at the end of each loop in Main_Loop. Overwrite this function in an inherited class for additional functionality."""
        pass

    def Main_Loop(self):
        plot=True
        play=False
        selected_directory_index=-1
        plot_styles=[self.plot_style for i in range(0,len(self.directories))]
        plot_widths=[self.plot_width for i in range(0,len(self.directories))]
        min_frame=self.frame_number
        max_frame=self.frame_number
        increment=0
        command_string=""

        while(plot):
            if(self.frame_number<0):
                print "Frame number can not be negative. Clamping to 0"
                self.frame_number=0
            min_frame,max_frame=self.Set_Max_Min(self.frame_number,min_frame,max_frame)

            (gnuplot_data_list,min_plot_value,max_plot_value)=self.Get_Gnuplot_Data_And_MinMax_Bounds_For_All_Directories(self.directories,self.variable_name,self.frame_number,plot_styles,plot_widths,self.stream_type,self.verbose_plot)
            self.data_range_min,=self.Set_Max_Min(value=min_plot_value,min_value=self.data_range_min)
            self.data_range_max,=self.Set_Max_Min(value=max_plot_value,max_value=self.data_range_max)     

            self.Set_Title()
            self.Set_Grid_Axis_Labels()
            self.Set_Data_Axis_Labels()
            self.Set_Grid_Axis_Range()
            self.Set_Data_Axis_Range()
            self.Plot(gnuplot_data_list)

            if(play):
                self.frame_number+=increment
                if((self.frame_number<=max_frame) & (self.frame_number>=min_frame)):
                    print "Playing: Current frame=%d, min_frame=%d, max_frame=%d"%(self.frame_number,min_frame,max_frame)
                    continue
                else:
                    play=False
                    self.frame_number-=increment

            if(self.command==""):
                input_string=raw_input('Enter command string...\n')
            else:
                input_string=self.command
                self.command=""

            if(input_string==""):
                self.frame_number+=increment
                continue

            match_shortoption=self.shortoption_re.search(input_string)
            match_longoption=self.longoption_re.search(input_string)

            command_string=input_string
            if(match_shortoption):
                command_string=self.shortoption_re.sub(r'-\g<sopt>',command_string)
            if(match_longoption):
                command_string=self.longoption_re.sub(r'--\g<lopt>',command_string)
            if((not match_shortoption) & (not match_longoption)):
                print "Can't parse input"

            (options,args)=self.parser_interactive.parse_args(command_string.split(" "))
            
            if(options.selected_directory_index!=None):
                selected_directory_index=options.selected_directory_index
                if(selected_directory_index==-1):
                    print "selecting all directories"
                else:
                    print "selecting %d th directory = %s"%(selected_directory_index,self.directories[selected_directory_index])
            if(options.fontsize!=None):
                self.fontsize=options.fontsize
            if(bool(options.write) | (options.filename_prefix_for_write!=None)):
                output_directory="%s/output/test_%d"%(self.base_directory_path,self.test_number)
                if(not os.path.exists(output_directory)):
                    try:
                        os.makedirs(output_directory,0755)
                        print "Directory %s created succesfully!"%output_directory
                    except:
                        print "Can't create directory %s"%output_directory
                        continue
                if(options.filename_prefix_for_write!=None): filename_base="%s/%s"%(output_directory,self.Remove_Preceding_Dashes(options.filename_prefix_for_write))
                else: filename_base="%s/plot_%s_%05d"%(output_directory,self.variable_name,self.frame_number)
                filename_ps=filename_base+".ps"
                filename_pdf=filename_base+".pdf"
                self.gnuplot_plotter.hardcopy(filename_ps,enhanced=0,color=1,fontsize=self.fontsize,solid=True)
                command_string_ps2pdf="ps2pdf %s %s"%(filename_ps,filename_pdf)

                done=False
                trial_number=0
                # The loop is to get around the problem when ps2pdf is not able to see the ps file written above. Probably Gnuplot.hardcopy returns before the file write is complete.
                while(not done):
                    process=Popen(command_string_ps2pdf.split(" "));
                    process.wait()
                    if(process.returncode==None):
                        print "Error: ps2pdf hasn't terminated"
                    elif(process.returncode==0):
                        print "ps2pdf completed succesfully"
                        done=True
                    else:
                        print "ps2pdf did not completed succesfully. Retrying"
                    trial_number+=1
                    if(trial_number>3):
                        done=True

            if(options.variable_name!=None):
                self.variable_name=self.Remove_Preceding_Dashes(options.variable_name)
            if(options.frame_number!=None):
                self.frame_number=options.frame_number
            if(options.increment!=None):
                self.frame_number+=options.increment
                increment=options.increment
            if(options.play):
                play=True
                if(increment==0): increment=1
                print "Playing from currentframe=%d to max_frame=%d"%(self.frame_number,max_frame)
            if(options.ls):
                print "-------------------"
                print "Directories:"
                dir_index=0
                for directory in self.directories:
                    print "[%d] %s"%(dir_index,self.Get_Directory_Base_Name(directory))
                    dir_index+=1
                print "-------------------"
            if(options.plot_style!=None):
                self.plot_style=self.Remove_Preceding_Dashes(options.plot_style)
                if(selected_directory_index!=-1):
                    plot_styles[selected_directory_index]=self.plot_style
                else:
                    for i in range(0,len(plot_styles)):plot_styles[i]=self.plot_style
            if(options.plot_width!=None):
                self.plot_width=options.plot_width
                if(selected_directory_index!=-1):
                    plot_widths[selected_directory_index]=self.plot_width
                else:
                    for i in range(0,len(plot_widths)):plot_widths[i]=self.plot_width
            if(options.quit):
                plot=False
            self.Load_Additional_Options(options)
            #else:
            #    print "Can't recognize option"

def main():
    p=PLOT(sys.argv[1:])
    p.Main_Loop()

if __name__=="__main__":
    main()
