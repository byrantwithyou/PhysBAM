#!/usr/bin/python

from physbam import *
from math import *
import Gnuplot, Gnuplot.funcutils
from plot import *

class CONVERGENCE(PLOT):
    def Get_Accuracy_Order(self,node_locations_all_resolutions,domain_indices_all_resolutions,list_array_all_resolutions,number_of_ghost_cells):
        number_of_real_cells_at_lowest_resolution=len(node_locations_all_resolutions[1])-2*number_of_ghost_cells
        print "number_of_real_cells_at_lowest_resolution=%d"%number_of_real_cells_at_lowest_resolution
        index_range=range(0,number_of_real_cells_at_lowest_resolution)
        accuracy_order_list=[0 for i in index_range]
        node_index_at_resolution=dict();node_location_at_resolution=dict();array_value_at_resolution=dict()
        for index in index_range:
            for resolution in [1,2,4]:
                node_index=resolution*index+number_of_ghost_cells
                node_index_at_resolution[resolution]=node_index
                node_location_at_resolution[resolution]=node_locations_all_resolutions[resolution][node_index]
                array_value_at_resolution[resolution]=list_array_all_resolutions[resolution][node_index]

            f_1x_minus_2x=array_value_at_resolution[1]-array_value_at_resolution[2]
            f_2x_minus_4x=array_value_at_resolution[2]-array_value_at_resolution[4]
            if((f_1x_minus_2x*f_2x_minus_4x>0)):
                accuracy_order=log(f_1x_minus_2x/f_2x_minus_4x,2)
                accuracy_order_list[index]=accuracy_order

        return accuracy_order_list

    def Get_Gnuplot_Data_And_MinMax_Bounds_For_All_Directories(self,directories,variable_name,frame_number,plot_styles,plot_widths,stream_type,verbose_plot):
        # Read data from all directories
        node_locations_all_resolutions=dict();domain_indices_all_resolutions=dict();list_array_all_resolutions=dict()
        for resolution in [1,2,4]:
            (node_locations_all_resolutions[resolution],domain_indices_all_resolutions[resolution],list_array_all_resolutions[resolution])=self.Read_Grid_And_Array_Data(self.directory_at_resolution[resolution],variable_name,frame_number,stream_type)
        # Compute accuracy order
        accuracy_order_list=self.Get_Accuracy_Order(node_locations_all_resolutions,domain_indices_all_resolutions,list_array_all_resolutions,self.number_of_ghost_cells)

        node_locations=node_locations_all_resolutions[1][self.number_of_ghost_cells:-self.number_of_ghost_cells]
        min_array_value=min(accuracy_order_list)
        max_array_value=max(accuracy_order_list)
        gnuplot_data=Gnuplot.Data(node_locations,accuracy_order_list,title="",with="%s lw %d"%(plot_styles[0],plot_widths[0]))

        if(self.plot_over_orig):# Get raw data from directories
            (gnuplot_data_list,min_plot_value,max_plot_value)=PLOT.Get_Gnuplot_Data_And_MinMax_Bounds_For_All_Directories(self,directories,variable_name,frame_number,plot_styles,plot_widths,stream_type,verbose_plot)
            min_plot_value,=self.Set_Max_Min(value=min_array_value,min_value=min_plot_value)
            max_plot_value,=self.Set_Max_Min(value=max_array_value,max_value=max_plot_value)
        else:
            gnuplot_data_list=[]
            min_plot_value=min_array_value
            max_plot_value=max_array_value
        gnuplot_data_list.append(gnuplot_data)

        average_order=sum(accuracy_order_list)/len(accuracy_order_list)
        print "average_order=",average_order
        return (gnuplot_data_list,min_plot_value,max_plot_value)

    def Setup_Parser(self,usage):
        usage="usage: %prog [options] <directory_1x> <directory_2x> <directory_4x>"
        PLOT.Setup_Parser(self,usage)
        self.parser.add_option("-g","--ghost",action="store",type="int",dest="number_of_ghost_cells",default=3,help="Number of ghost cells. Default=%default")
        self.parser.add_option("--plot_over_orig",action="store_true",dest="plot_over_orig",default=False,help="Plot convergence data over original data")

    def Parse_Command_Line(self,argument_list):
        (options,args)=PLOT.Parse_Command_Line(self,argument_list)
        self.number_of_ghost_cells=options.number_of_ghost_cells
        self.plot_over_orig=options.plot_over_orig

        self.directory_at_resolution=dict()
        if(len(args)==3):
            (self.directory_at_resolution[1],self.directory_at_resolution[2],self.directory_at_resolution[4])=args
        else:
            self.parser.error("incorrect number of arguments")
        return (options,args)
            
    def Print_Arguments(self):
        PLOT.Print_Arguments(self)
        print "------------ CONVERGENCE PARAMS ------------"
        print "number_of_ghost_cells =",self.number_of_ghost_cells
        print "plot_over_orig =",self.plot_over_orig
        print "directory_at_resolution=",self.directory_at_resolution
        print "--------------------------------------------"

    def __init__(self,argument_list):
        PLOT.__init__(self,argument_list)
        if(self.verbose_plot):
            self.gnuplot_plotter.ylabel("Accuracy order for "+self.variable_name)
        if(not(self.plot_over_orig)):
            self.padding_factor=0

def main():
    convergence_plotter=CONVERGENCE(sys.argv[1:])
    convergence_plotter.Main_Loop()

if __name__=="__main__":
    main()
