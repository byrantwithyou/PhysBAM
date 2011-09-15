#!/usr/bin/python

from plot import *

class PLOT_2D_SLICE(PLOT):
    def PhysBAM_Array_2D_Into_List(self,array_physbam):
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

    def Get_Node_Locations_from_PhysBAM_Grid_2D(self,grid_physbam,node_indices):
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
  
    def Read_Grid_And_Array_Data_Single_Variable(self,directory,variable_name,frame_number,stream_type):
        grid_file=directory+"/grid.%d.gz"%frame_number
        array_file=directory+"/"+variable_name+".%d.gz"%frame_number

        if(stream_type=='float'): grid_physbam=GRID_Vf2()
        else: grid_physbam=GRID_Vd2()
        Read_From_File(stream_type,grid_file,grid_physbam)

        if(stream_type=='float'): array_physbam=ARRAYS_2D_f()
        else: array_physbam=ARRAYS_2D_d()
        Read_From_File(stream_type,array_file,array_physbam)

        (domain_indices,list_array)=self.PhysBAM_Array_2D_Into_List(array_physbam)
        node_locations=self.Get_Node_Locations_from_PhysBAM_Grid_2D(grid_physbam,domain_indices)

        return (node_locations,list_array)

    def Read_Grid_And_Array_Data(self,directory,variable_name,frame_number,stream_type):
        if(variable_name!="adiabatic"): # plot adiabatic constant
            return self.Read_Grid_And_Array_Data_Single_Variable(directory,variable_name,frame_number,stream_type)

        (node_locations,density_array)=self.Read_Grid_And_Array_Data_Single_Variable(directory,"density",frame_number,stream_type)
        (node_locations,pressure_array)=self.Read_Grid_And_Array_Data_Single_Variable(directory,"pressure",frame_number,stream_type)

        enthalpy_array=density_array # initialization
        for index in range(len(enthalpy_array)):
            enthalpy_array[index]=pressure_array[index]/pow(density_array[index],1.4)
        return (node_locations,enthalpy_array)

    def Setup_Parser(self,usage):
        PLOT.Setup_Parser(self,usage)
        self.parser.add_option("-x",action="store_const",const=1,dest="plot_axis",help="plot x axis")
        self.parser.add_option("-y",action="store_const",const=2,dest="plot_axis",help="plot y axis")
        self.parser.add_option("--slice",action="store",type="int",dest="slice_index",default=1,help="slice index of the non plotted axis. Default=%default")

    def Parse_Command_Line(self,argument_list):
        (options,args)=PLOT.Parse_Command_Line(self,argument_list)
        if(options.plot_axis!=None): self.plot_axis=options.plot_axis
        else: self.plot_axis=1
        self.slice_index=options.slice_index

        return (options,args)

    def Print_Arguments(self):
        PLOT.Print_Arguments(self)
        print "------------ PLOT_2D_SLICE PARAMS ------------"
        print "plot_axis =",self.plot_axis
        print "slice_index =",self.slice_index
        print "--------------------------------------------"

    def Set_Labels(self):
        PLOT.Set_Labels(self)
        if(self.verbose_plot):
            if(self.plot_axis==1): xlabel="x"
            else: xlabel="y"
            self.gnuplot_plotter.xlabel(xlabel+" (slice %d)"%self.slice_index)

    def __init__(self,argument_list):
        PLOT.__init__(self,argument_list)

        self.parser_interactive.add_option("-x",action="store_const",const=1,dest="plot_axis",help="plot x axis")
        self.parser_interactive.add_option("-y",action="store_const",const=2,dest="plot_axis",help="plot y axis")
        self.parser_interactive.add_option("--slice",action="store",type="int",dest="slice_index",help="slice index of the non plotted axis")

    def Load_Additional_Options(self,options):
        if(options.plot_axis!=None):
            self.plot_axis=options.plot_axis
        if(options.slice_index!=None):
            self.slice_index=options.slice_index

def main():
    plot_2d_slice=PLOT_2D_SLICE(sys.argv[1:])
    plot_2d_slice.Main_Loop()

if __name__=="__main__":
    main()

