#!/usr/bin/python

from plot import *

class PLOT_2D(PLOT):

    # Slice Mode Functions
    def PhysBAM_Array_2D_Into_List_Slice(self,array_physbam):
        """Takes a physbam array and returns domain indices and a python list containing array data"""

        domain_indices_box=array_physbam.Domain_Indices()

        domain_indices=[]
        if(self.plot_axis==1): domain_indices=range(domain_indices_box.min_corner.x,domain_indices_box.max_corner.x+1)
        else: domain_indices=range(domain_indices_box.min_corner.y,domain_indices_box.max_corner.y+1)

        list_array=[0 for i in domain_indices]
        array_index=Vi2(self.slice_index,self.slice_index)
        list_index=0
        for i in domain_indices:
            if(self.plot_axis==1): array_index.x=i
            else: array_index.y=i
            list_array[list_index]=array_physbam.__getitem__(array_index)
            list_index+=1
        return (domain_indices,list_array)

    def Get_Node_Locations_from_PhysBAM_Grid_2D_Slice(self,grid_physbam,node_indices):
        """Reads a physbam grid and returns a python list of the node positions"""

        node_locations=[0 for i in node_indices]
        array_index=Vi2(self.slice_index,self.slice_index)
        list_index=0
        for i in node_indices:
            if(self.plot_axis==1):
                array_index.x=i
                node_locations[list_index]=grid_physbam.X(array_index).x
            else:
                array_index.y=i
                node_locations[list_index]=grid_physbam.X(array_index).y
            list_index+=1
        return node_locations
  
    def Read_Grid_And_Array_Data_Single_Variable_Slice(self,directory,variable_name,frame_number,stream_type):
        grid_file=directory+"/grid.%d.gz"%frame_number
        array_file=directory+"/"+variable_name+".%d.gz"%frame_number

        if(stream_type=='float'): grid_physbam=GRID_Vf2()
        else: grid_physbam=GRID_Vd2()
        Read_From_File(stream_type,grid_file,grid_physbam)

        if(stream_type=='float'): array_physbam=ARRAYS_2D_f()
        else: array_physbam=ARRAYS_2D_d()
        Read_From_File(stream_type,array_file,array_physbam)

        (domain_indices,list_array)=self.PhysBAM_Array_2D_Into_List_Slice(array_physbam)
        node_locations=self.Get_Node_Locations_from_PhysBAM_Grid_2D_Slice(grid_physbam,domain_indices)

        return (node_locations,list_array)

    def Read_Grid_And_Array_Data_Slice(self,directory,variable_name,frame_number,stream_type):
        if(variable_name!="adiabatic"): # plot adiabatic constant
            return self.Read_Grid_And_Array_Data_Single_Variable_Slice(directory,variable_name,frame_number,stream_type)

        (node_locations,density_array)=self.Read_Grid_And_Array_Data_Single_Variable_Slice(directory,"density",frame_number,stream_type)
        (node_locations,pressure_array)=self.Read_Grid_And_Array_Data_Single_Variable_Slice(directory,"pressure",frame_number,stream_type)

        enthalpy_array=density_array # initialization
        for index in range(len(enthalpy_array)):
            enthalpy_array[index]=pressure_array[index]/pow(density_array[index],1.4)
        return (node_locations,enthalpy_array)



    # Non-Slice Mode Functions
#    def PhysBAM_Array_2D_Into_List(self,array_physbam):
#        """Takes a physbam array and returns domain indices and a python list containing array data"""
#
#        domain_indices_box=array_physbam.Domain_Indices()
#
#        domain_indices_x=range(domain_indices_box.min_corner.x,domain_indices_box.max_corner.x+1)
#        domain_indices_y=range(domain_indices_box.min_corner.y,domain_indices_box.max_corner.y+1)
#
#        list_array=[[0 for j in domain_indices_y] for i in domain_indices_x]
#        array_index=Vi2()
#        list_index_x=0
#        for i in domain_indices_x:
#            array_index.x=i
#            list_index_y=0
#            for j in domain_indices_y:
#                array_index.y=j
#                list_array[list_index_x][list_index_y]=array_physbam.__getitem__(array_index)
#                list_index_y+=1
#            list_index_x+=1
#        return (domain_indices_x,domain_indices_y,list_array)
#
#    def Get_Node_Locations_from_PhysBAM_Grid_2D(self,grid_physbam,node_indices_x,node_indices_y):
#        """Reads a physbam grid and returns a python list of the node positions"""
#
#        node_locations_x=[0 for i in node_indices_x]
#        node_locations_y=[0 for j in node_indices_y]
#        array_index=Vi2(node_indices_x[0],node_indices_y[0]) #initialize to have a valid index
#
#        # Get the x,y locations. Assumes cartesian grid, so can iterate independently
#        # Iterate and get all the x locations
#        list_index=0
#        for i in node_indices_x:
#            array_index.x=i
#            node_locations_x[list_index]=grid_physbam.X(array_index).x
#            list_index+=1
#
#        # Iterate and get all the y locations
#        list_index=0
#        for i in node_indices_y:
#            array_index.y=i
#            node_locations_y[list_index]=grid_physbam.X(array_index).y
#            list_index+=1
#
#        return (node_locations_x,node_locations_y)
#
#    def Read_Grid_And_Array_Data(self,directory,variable_name,frame_number,stream_type):
#        grid_file=directory+"/grid.%d.gz"%frame_number
#        array_file=directory+"/"+variable_name+".%d.gz"%frame_number
#
#        if(stream_type=='float'): grid_physbam=GRID_Vf2()
#        else: grid_physbam=GRID_Vd2()
#        Read_From_File(stream_type,grid_file,grid_physbam)
#
#        if(stream_type=='float'): array_physbam=ARRAYS_2D_f()
#        else: array_physbam=ARRAYS_2D_d()
#        Read_From_File(stream_type,array_file,array_physbam)
#
#        (domain_indices_x,domain_indices_y,list_array)=self.PhysBAM_Array_2D_Into_List(array_physbam)
#        (node_locations_x,node_locations_y)=self.Get_Node_Locations_from_PhysBAM_Grid_2D(grid_physbam,domain_indices_x,domain_indices_y)
#
#        return (node_locations_x,node_locations_y,list_array)

    def Read_Grid_And_Array_Data_From_Gnuplot_File(self,gnuplot_file_name):
        # assumes that file writes data in x,y,data format in a column by column order
        gnuplot_file=open(gnuplot_file_name,'r')
        node_locations_x=[]
        node_locations_y=[]
        list_array=[]
        list_array_1d=[]
        node_location_x=0
        node_location_y=0
        first_column=True
        for line in gnuplot_file:
            line=line.strip()
            if(not line):
                #end of column
                list_array.append(list_array_1d[:]) #used [:] to copy value rather than reference
                node_locations_x.append(float(node_location_x)) # new x value in each column
                del list_array_1d[:]
                first_column=False
                continue
            (node_location_x,node_location_y,value)=line.split()
            if(first_column): node_locations_y.append(float(node_location_y))#Each column has all y values, use the first one.
            list_array_1d.append(float(value))

        gnuplot_file.close()
        return (node_locations_x,node_locations_y,list_array)

    # Check's whether to get slice mode or full data
    def Get_Gnuplot_Data_And_MinMax_Bounds(self,directory,variable_name,frame_number,plot_style,plot_width,stream_type,verbose_plot):
        plot_title=self.Get_Plot_Title(directory,verbose_plot)
        if(self.mode==1): # slice mode
            (node_locations,list_array)=self.Read_Grid_And_Array_Data_Slice(directory,variable_name,frame_number,stream_type)
            gnuplot_data=Gnuplot.Data(node_locations,list_array,title=plot_title,with="%s lw %d"%(plot_style,plot_width))
        else: #full mode
            (node_locations_x,node_locations_y,list_array)=self.Read_Grid_And_Array_Data(directory,variable_name,frame_number,stream_type)
            gnuplot_data=Gnuplot.GridData(list_array,node_locations_x,node_locations_y,title=plot_title,with="%s lw %d"%(plot_style,plot_width))

        min_array_value=self.min_recursive(list_array)
        max_array_value=self.max_recursive(list_array)
        if(not(self.set_input_minmax_xrange)): # If not given as input set the xrange here
            if(self.mode==1):
                self.minxrange=node_locations[0]
                self.maxxrange=node_locations[-1]
            else:
                self.minxrange=node_locations_x[0]
                self.maxxrange=node_locations_x[-1]
        if(not(self.set_input_minmax_yrange)): # If not given as input set the yrange here
            if(self.mode==1):
                pass
            else:
                self.minyrange=node_locations_y[0]
                self.maxyrange=node_locations_y[-1]

        return (gnuplot_data,min_array_value,max_array_value)

    def min_recursive(self,complex_list):
        if(not isinstance(complex_list,list)): # This is a single element. Just return that
            return complex_list

        minlist=[self.min_recursive(l) for l in complex_list]
        return min(minlist)

    def max_recursive(self,complex_list):
        if(not isinstance(complex_list,list)): # This is a single element. Just return that
            return complex_list

        maxlist=[self.max_recursive(l) for l in complex_list]
        return max(maxlist)

    def Set_Grid_Axis_Labels(self):
        if(self.verbose_plot):
            if(self.mode==1):
                if(self.plot_axis==1): xlabel="x"
                else: xlabel="y"
                self.gnuplot_plotter.xlabel(xlabel+" (slice %d)"%self.slice_index)
            else:
                self.gnuplot_plotter.xlabel("x")
                self.gnuplot_plotter.ylabel("y")

    def Set_Data_Axis_Labels(self):
        if(self.mode==1):
            return PLOT.Set_Data_Axis_Labels(self)

        if(self.verbose_plot):
            self.gnuplot_plotter("set zlabel \"%s\""%self.variable_name)

    def Set_Grid_Axis_Range(self):
        PLOT.Set_Grid_Axis_Range(self)

        if(self.mode!=1):
            yrange_set_string="set yrange [%f:%f]"%(self.minyrange,self.maxyrange)
            self.gnuplot_plotter(yrange_set_string)

    def Set_Data_Axis_Range(self):
        if(self.mode==1): # slice mode
            return PLOT.Set_Data_Axis_Range(self)

        data_range_paddding=(self.data_range_max-self.data_range_min)*self.padding_factor
        if(data_range_paddding!=0):
            data_range_set_string="set zrange [%f:%f]"%(self.data_range_min-data_range_paddding,self.data_range_max+data_range_paddding)
            self.gnuplot_plotter(data_range_set_string)

    def Plot(self,gnuplot_data_list):
        if(self.mode==1): # slice mode
            return PLOT.Plot(self,gnuplot_data_list)

        self.gnuplot_plotter.splot(*gnuplot_data_list)

    def Load_Additional_Options(self,options):
        if(options.mode!=None):
            self.mode=options.mode
            self.Setup_Gnuplot()
        if(options.plot_axis!=None):
            self.plot_axis=options.plot_axis
        if(options.slice_index!=None):
            self.slice_index=options.slice_index

    def Setup_Parser(self,usage):
        PLOT.Setup_Parser(self,usage)

        self.parser.add_option("-m","--mode",action="store",type="int",dest="mode",help="plotting mode: (1:slice,2:contour,3:heightfield)")
        self.parser.add_option("--minyrange",action="store",type="float",dest="minyrange",help="set minimum y range. If set, maximum y range will default to 1")
        self.parser.add_option("--maxyrange",action="store",type="float",dest="maxyrange",help="set maximum y range. If set, minimum y range will default to 0")
        self.parser.add_option("-x",action="store_const",const=1,dest="plot_axis",help="plot x axis")
        self.parser.add_option("-y",action="store_const",const=2,dest="plot_axis",help="plot y axis")
        self.parser.add_option("--slice",action="store",type="int",dest="slice_index",default=1,help="slice index of the non plotted axis. Default=%default")

    def Parse_Command_Line(self,argument_list):
        (options,args)=PLOT.Parse_Command_Line(self,argument_list)

        self.mode=2
        if(options.mode!=None):
            self.mode=options.mode

        self.minyrange=0
        self.maxyrange=1
        if((options.minyrange!=None) | (options.maxyrange!=None)): #Set default values
            self.set_input_minmax_yrange=True
        else:
            self.set_input_minmax_yrange=False
        if(options.minyrange!=None):
            self.minyrange=options.minyrange
        if(options.maxyrange!=None):
            self.maxyrange=options.maxyrange

        if(options.plot_axis!=None): self.plot_axis=options.plot_axis
        else: self.plot_axis=1
        self.slice_index=options.slice_index


        return (options,args)

    def Print_Arguments(self):
        PLOT.Print_Arguments(self)
        print "------------ PLOT_2D PARAMS ------------"
        print "mode=",self.mode
        print "minyrange=",self.minyrange
        print "maxyrange=",self.maxyrange
        print "plot_axis =",self.plot_axis
        print "slice_index =",self.slice_index
        print "--------------------------------------------"
 
    def Setup_Gnuplot_For_Contours(self):
        # prams for contours
        #self.gnuplot_plotter("set size ratio -1")
        #self.gnuplot_plotter("set grid")
        self.gnuplot_plotter("set cntrparam levels incremental 2.568e-1,.19367,6.067") # flow past a step
        #self.gnuplot_plotter("set cntrparam levels incremental 1.731,.6617,20.92") # double mac
        self.gnuplot_plotter("set contour")
        self.gnuplot_plotter("unset clabel")
        self.gnuplot_plotter("set nosurface")
        self.gnuplot_plotter("set view 0,0,1.5,1")
    
    def Setup_Gnuplot_For_Slice(self):
        self.gnuplot_plotter("unset contour")

    def Setup_Gnuplot_For_Height_Field(self):
        self.gnuplot_plotter("unset contour")
        self.gnuplot_plotter("set surface")
        self.gnuplot_plotter("set view 60,30,1,1")

    def Setup_Gnuplot(self):
        if(self.mode==1):
            self.Setup_Gnuplot_For_Slice()
        elif(self.mode==2):
            self.Setup_Gnuplot_For_Contours()
        elif(self.mode==3):
            self.Setup_Gnuplot_For_Height_Field()
        else:
            print "UNKNOWN MODE %s. Switching to slice mode"%self.mode
            self.mode=1
            self.Setup_Gnuplot_For_Slice()

    def __init__(self,argument_list):
        PLOT.__init__(self,argument_list)

        self.dimension=2

        self.parser_interactive.add_option("-m","--mode",action="store",type="int",dest="mode",help="plotting mode: (1:slice,2:contour,3:heightfield)")
        self.parser_interactive.add_option("-x",action="store_const",const=1,dest="plot_axis",help="plot x axis")
        self.parser_interactive.add_option("-y",action="store_const",const=2,dest="plot_axis",help="plot y axis")
        self.parser_interactive.add_option("--slice",action="store",type="int",dest="slice_index",help="slice index of the non plotted axis")
        self.Setup_Gnuplot()

def main():
    plot_2d=PLOT_2D(sys.argv[1:])
    plot_2d.Main_Loop()

if __name__=="__main__":
    main()
