#!/usr/bin/python
import sys
import os
import optparse
import physbam

parser = optparse.OptionParser("usage: %prog [-s start_frame] [-e end_frame] [-o output_file] mtn_file")
parser.add_option('-s','--start_frame',default=0,type='int',help='first frame of the new motion')
parser.add_option('-e','--end_frame',default=0,type='int',help='last frame of the new motion')
parser.add_option('-o','--output',default="",type='string',help='output filename')
(options,args)=parser.parse_args()
if len(args)<1: parser.error("invalid number of arguments")

if options.output is "": options.output=args[0]+"_crop.mtn"
body_motion=physbam.BODY_MOTION_SEQUENCE_f()
physbam.Read_From_File("float",args[0]+".mtn",body_motion)
start_frame=options.start_frame
end_frame=options.end_frame
if start_frame is 0: start_frame=body_motion.saved_frame
if start_frame is 0: start_frame=1
body_motion.Resize(start_frame,end_frame)
body_motion.saved_frame=0
physbam.Write_To_File("float",options.output,body_motion)
print "Wrote to "+options.output
