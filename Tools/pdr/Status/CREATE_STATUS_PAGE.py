#!/usr/bin/python
import os
import sys
from stat import *
from datetime import datetime
from time import *
import traceback

class CREATE_STATUS_PAGE:
    def find_percent_string(self,job_directory,frame_no,in_progress):
        try:
            cout=" ".join(open(os.path.join(job_directory,"Status",frame_no+".out")).readlines())
            end=cout.rfind("</SLAVE>")
            last_command="N/A"
            if end!=-1:
                start=cout.rfind("<SLAVE>",0,end)+7
                if start!=-1:
                    last_command=cout[start:end]
            end=cout.rfind("%")
            percentage="N/A"
            if end!=-1:
                start=cout.rfind(" ",0,end)+1
                if start!=-1:
                    percentage=cout[start:end]
            return last_command,percentage
        except:
            if in_progress==1:
                return "Contacting slave to start script","N/A"
            else:
                return "Timed out (check if slavename is correct and there is enough disk space)","N/A"

    def process_frames(self,job_directory,frame_directory,in_progress):
       job_page_string="<table cellpadding=\"0\" cellspacing=\"0\">\n"
       job_page_string+="<tr>"
       job_page_string+="<th width=\"70\" valign=\"top\"><b>Frame</b></td>"
       job_page_string+="<th width=\"258\" valign=\"top\"><b>Slave Machine</b></td>"
       job_page_string+="<th width=\"258\" valign=\"top\"><b>Script Status</b></td>"
       job_page_string+="<th width=\"70\" valign=\"top\"><b>% Done</b></td>"
       job_page_string+="<th width=\"200\" valign=\"top\">&nbsp;</td>"
       if in_progress==1:
           job_page_string+="<th width=\"150\" valign=\"top\"><b>Time Elapsed</b></td>"
       job_page_string+="</tr>\n"
       full_frame_directory=os.path.join(job_directory,frame_directory)
       current_time=time()

       
       frame_item_list=[]
       for frame in os.listdir(full_frame_directory):
           start_time=0
           try:
               start_time=open(os.path.join(job_directory,"Status",frame+".start")).readlines()[0]
           except:
               pass
           slave_computer=open(os.path.join(full_frame_directory,frame)).readlines()[0]
           frame_item_list.append([int(frame),slave_computer,start_time])
       frame_item_list.sort()
                
       for frame_item in frame_item_list:
           frame=str(frame_item[0])
           slave_computer=frame_item[1]
           start_time=frame_item[2]
           elapsed_time="N/A"
           try:
               if start_time!=0:
                   elapsed_time=strftime("%H:%M:%S",gmtime(current_time-float(start_time)))
           except:
               print "Failed to get time for frame '%s'"%frame
               pass
           last_command,percentage=self.find_percent_string(job_directory,frame,in_progress)
           percentage_bar=""
           try:
               percentage_bar=("<table class=progress cellpadding=\"0\" cellspacing=\"0\" border=0 cellspacing=0 cellpadding=0 width=200><tr><td class=progress bgcolor=blue height=10 width=%d></td><td bgcolor=white width=%d class=progress></td></tr></table>"%(int(percentage)*2,(100-int(percentage))*2))
           except:
               pass
           job_page_string+="<tr><td><font size=-1>"+frame+"</font></td>"
           job_page_string+="<td><font size=-1>"+slave_computer+"</font></td>"
           job_page_string+="<td><font size=-1>"+last_command+"</font></td>"
           job_page_string+="<td><font size=-1>"+percentage+"</font></td>"
           job_page_string+="<td><font size=-1>"+percentage_bar+"</font></td>"
           if in_progress==1:
               job_page_string+="<td><font size=-1>"+elapsed_time+"</font></td></tr>\n"
       job_page_string+="</table>\n"
       return job_page_string,len(frame_item_list)


    def __init__(self,server_directory,destination_directory):
        print "Running..."
        self.server_directory=server_directory
        self.destination_directory=destination_directory
        self.jobs_directory=os.path.join(self.server_directory,"State","Jobs")
        # go through each job
        job_info_list=[]
        for i in os.listdir(self.jobs_directory):
            job_directory=os.path.join(self.jobs_directory,i)            
            infofilename=os.path.join(job_directory,"Job_Info")
            
            number_in_progress=0
            number_failed=0
            number_done=0
            number_pending=0
            
            try:
                infolines=[s.replace("\n","") for s in open(infofilename).readlines()]
                jobname,username,frames=infolines[:3]
                job_page_string="<html><head>\n"
                job_page_string+="<title>PhysBAM Distributed Renderer - Server Status</title>\n"
                job_page_string+="<META HTTP-EQUIV=\"Refresh\" CONTENT=\"5\">\n"
                job_page_string+="<link href=\"style.css\" rel=\"stylesheet\" type=\"text/css\">"
                job_page_string+="</head><body bgcolor=\"ffffff\">"
                job_page_string+="<h1>PhysBAM Distributed Renderer - Server Status</h1>\n"
                job_page_string+="<table width=700><tr>"
                job_page_string+="<th width=100>Job Name</th><td width=100>"+jobname+"</td><th width=100>Job ID</th><td width=100>"+i+"</td>"
                try:
                    job_page_string+="<th width=100>Last Update</th><td width=200>%s</td>"%strftime("%a, %d %b %Y %H:%M:%S", localtime())
                except:
                    print "Failed to get time for job '%s'"%i
                job_page_string+="</tr></table>\n"
                

                # in progress frames
                job_page_string+="<h3>In Progress Frames:</h3>\n"
                new_string,number_in_progress=self.process_frames(job_directory,"In_Progress_Frames",1)
                job_page_string+=new_string
                # failed frames
                job_page_string+="\n<h3>Failed Frames:</h3>\n"
                new_string,number_failed=self.process_frames(job_directory,"Failed_Frames",0)
                job_page_string+=new_string
                # done
                job_page_string+="\n<h3>Done Frames:</h3>\n"
                done_frames=[]
                for frame in os.listdir(os.path.join(job_directory,"Done_Frames")):
                    done_frames.append(int(frame))
                done_frames.sort()
                for frame in done_frames:
                    job_page_string+="  "+str(frame)
                number_done=len(done_frames)
                # pending
                job_page_string+="\n<h3>Pending Frames:</h3>\n"
                pending_frames=[]
                for frame in os.listdir(os.path.join(job_directory,"Pending_Frames")):
                    pending_frames.append(int(frame))
                pending_frames.sort()
                for frame in pending_frames:
                    job_page_string+="  "+str(frame)
                number_pending=len(pending_frames)
                
                # finish up the page string
                job_page_string+="</body></html>\n"

                # write the string to disk
                open(os.path.join(self.destination_directory,i+".html"),"w").write(job_page_string)
            except:
                traceback.print_exc(file=sys.stdout)
                print "Failed creating page for job '%s'"%i
                
            # remember the job and it's creation time so that it can be sorted later
            job_info=[os.stat(job_directory)[ST_MTIME],infofilename,i,number_in_progress,number_failed,number_done,number_pending]
            job_info_list.append(job_info)


        #sort the jobs so that they appear in sorted order (based on last time modified)
        self.main_page_string="<html><head>"
        self.main_page_string+="<title>PhysBAM Distributed Renderer - Server Status</title>"
        self.main_page_string+="<link href=\"style.css\" rel=\"stylesheet\" type=\"text/css\">"

        #self.main_page_string+="<META HTTP-EQUIV=\"Refresh\" CONTENT=\"5\">\n"
        self.main_page_string+="</head><body bgcolor=\"ffffff\">"
        self.main_page_string+="<h1>PhysBAM Distributed Renderer - Server Status</h1>\n"
        
        self.main_page_string+="<b>Status last updated: "+strftime("%a, %d %b %Y %H:%M:%S", localtime())+"</b><br>\n"
        self.main_page_string+="<h3>Current Jobs:</h3>\n"
        self.main_page_string+="<table cellpadding=\"0\" cellspacing=\"0\">\n"
        self.main_page_string+="<tr>"
        self.main_page_string+="<th width=\"70\" valign=\"top\"><b>Job ID</b></td>"
        self.main_page_string+="<th width=\"258\" valign=\"top\"><b>Job Name</b></td>"
        self.main_page_string+="<th width=\"150\" valign=\"top\"><b>User Name</b></td>"
        self.main_page_string+="<th width=\"50\" valign=\"top\"><b>Priority</b></td>"
        self.main_page_string+="<th width=\"100\" valign=\"top\"><b>Status</b></td>"
        self.main_page_string+="<th width=\"70\" valign=\"top\"><b># Pending</b></td>"
        self.main_page_string+="<th width=\"100\" valign=\"top\"><b># In Progress</b></td>"
        self.main_page_string+="<th width=\"70\" valign=\"top\"><b># Failed</b></td>"
        self.main_page_string+="<th width=\"70\" valign=\"top\"><b># Done</b></td>"
        self.main_page_string+="<th width=\"150\" valign=\"top\"><b>Last Activity</b></td>"
        self.main_page_string+="</tr>\n"
        job_info_list.sort()
        job_info_list.reverse()
        for i in job_info_list:
            self.main_page_string+="<tr>\n"
            modified_date,infofilename,job_web_page,number_in_progress,number_failed,number_done,number_pending=i[:7]
            def Make_Percentage_Bar_String(numbers,colors):
                total=reduce(lambda x,y:x+y,numbers)
                if total!=0: percentages=map(lambda x: int(float(x)/float(total)*100),numbers)
                else: percentages=map(lambda x:0,numbers)
                table_cells="".join(map(lambda percent,color: ("<td class=\"progress\" bgcolor=%s height=15 width=%d></td>"%(color,percent)), percentages, colors))
                return ("<table class=\"progress\"><tr>%s</tr></table>"%table_cells)
            bar_graph=Make_Percentage_Bar_String( (number_pending,number_in_progress,number_failed,number_done), ("#cccccc","#8888ee","#ee8888","#88ee88") )
            infolines=[s.replace("\n","") for s in open(infofilename).readlines()]
            jobname,username,priority=infolines[:3]
            time_modified=datetime.fromtimestamp(modified_date).strftime('%H:%M:%S %m/%d/%y')
            self.main_page_string+="<td>"+job_web_page+"</td>"
            self.main_page_string+="<td><a href="+job_web_page+".html>"+jobname+"</a></td>"
            self.main_page_string+="<td>"+username+"</td>"
            self.main_page_string+="<td>"+priority+"</td>"
            self.main_page_string+="<td>"+bar_graph+"</td>"
            self.main_page_string+="<td bgcolor=#cccccc>"+str(number_pending)+"</td>"
            self.main_page_string+="<td bgcolor=#8888ee>"+str(number_in_progress)+"</td>"
            self.main_page_string+="<td bgcolor=#ee8888>"+str(number_failed)+"</td>"
            self.main_page_string+="<td bgcolor=#88ee88>"+str(number_done)+"</td>"
            self.main_page_string+="<td>"+time_modified+"</td>\n"
            self.main_page_string+="</tr>\n"
        self.main_page_string+="</table>\n"

        # add all the free slaves to the main page
        available_slaves=[s.replace("\n","") for s in open(os.path.join(self.server_directory,"State","Free_Slaves")).readlines()]
        number_of_available_slaves=len(available_slaves)
        if available_slaves[0]=='':
            number_of_available_slaves-=1
        self.main_page_string+="<h3>Currently Free Slaves (Count: %s"%number_of_available_slaves+"):</h3>\n"
        available_slaves.sort()
        for available_slave in available_slaves:
            self.main_page_string+=available_slave+"<br>\n"
        
        # add all the linux slaves to the main page
        linux_machines=os.listdir(os.path.join(self.server_directory,"State","Slaves"))
        self.main_page_string+="<h3>Registered Linux Slaves (Count: %s"%len(linux_machines)+"):</h3>\n"
        linux_machines.sort()
        for linux_machine in linux_machines:
            platform,priority,memory,host,instance=(linux_machine).split("__")
            self.main_page_string+=host+"--"+instance+"<br>\n"

        # write out the main page
        open(os.path.join(self.destination_directory,"index.html"),"w").write(self.main_page_string)



if __name__=="__main__":
    try:
        executable,server_directory,web_directory=sys.argv
    except:
        print "Usage: %s server_directory web_directory"%sys.argv[0]
    status=CREATE_STATUS_PAGE(server_directory,web_directory)   
       
