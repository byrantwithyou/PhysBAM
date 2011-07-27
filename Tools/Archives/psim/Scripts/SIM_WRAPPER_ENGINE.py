#!/usr/bin/python
######################################################################
# Copyright 2005, Andrew Selle.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
#######################################################################
# Classes that help create a wrapper script
#######################################################################
import os
import sys
import re
import time
import threading
import popen2
from psim.Base import SOCKET

#######################################################################
# Class DUMMY_RPC for testing
#######################################################################
class DUMMY_RPC:
    def Update_Status(self,session,key,value):
        print "Updating on session=%d  %s=%s"%(session,key,repr(value))

#######################################################################
# Class PRODUCER_CONSUMER_QUEUE
#######################################################################
class PRODUCER_CONSUMER_QUEUE:
    def __init__(self):
        self.events=[]
        self.mutex=threading.Lock()

    def Produce(self,data):
        self.mutex.acquire()
        self.events.append(data)
        self.mutex.release()

    def Consume(self):
        self.mutex.acquire()
        element=self.events.pop()
        self.mutex.release()
        return element

    def Consume_All(self):
        self.mutex.acquire()
        elements=self.events
        self.events=[]
        self.mutex.release()
        return elements

#######################################################################
# Exception EXIT_REPORTING_EXCEPTION
#######################################################################
class EXIT_REPORTING_EXCEPTION(Exception):
    def __init__(self):
        pass
#######################################################################
# Class SIM_WRAPPER
#######################################################################
class SIM_WRAPPER:
    def __init__(self,cmd):
        self.cmd=cmd
        self.simulation_events=[]
        self.polled_events=[]
        self.simulation_event_queue=PRODUCER_CONSUMER_QUEUE()
        self.start_time=time.time()
        self.poll_interval=5 # 5 seconds
        self.session=int(os.environ["PSIM_SESSION_ID"])

    def run(self):
        # make the reporting thread
        self.reporting_thread=threading.Thread()
        self.reporting_thread.run=self.reporting
        self.reporting_thread.start()
        
        # launch the process and capture the pid
        pipe=popen2.Popen4(self.cmd)
        self.input,self.pid=pipe.fromchild,pipe.pid

        while 1:
            line=self.input.readline()
            if line=="": break # cmd died, break to find out why
            print line,
            line=line.strip()
            for regexp,match_func,act_func in self.simulation_events:
                match=regexp.match(line)
                try:
                    event_data=match_func(match)
                except:
                    pass
                else:
                   self.simulation_event_queue.Produce((act_func,event_data))
        # end of process detect exit code
        exit_code=pipe.wait()
        print "got exit code %d"%exit_code
        self.simulation_event_queue.Produce((self.report_exit,(exit_code,time.time())))
        self.reporting_thread.join()

    def report_exit(self,pid,rpc,session,data):
        exit_code,exit_time=data
        rpc.Update_Status(session,"exit_code",exit_code)
        rpc.Update_Status(session,"exit_time",exit_time)
        rpc.Update_Status(session,"total_duration",exit_time-self.start_time)
        raise EXIT_REPORTING_EXCEPTION
    
    def reporting(self):
        try:
            while 1:
                time.sleep(self.poll_interval)
                try:
                    rpc_connection=SOCKET.CLIENT(os.environ["PSIM_SERVER_HOST"],int(os.environ["PSIM_SERVER_PORT"]),
                                                (os.environ["PSIM_CLIENT_KEY"],os.environ["PSIM_CLIENT_CERT"],os.environ["PSIM_CA_CERT"]))
                    #rpc_connection=DUMMY_RPC()
                except:
                    print "Failed to get rpc_connection for reporting"
                    pass
                else:
                    events=self.simulation_event_queue.Consume_All()
                    for event_function,event_data in events:
                        print "events func=%s data=%s"%(repr(event_function),repr(event_data))
                        try:
                            ret=event_function(self.pid,rpc_connection,self.session,event_data)
                            print "ret=%s"%ret
                        except EXIT_REPORTING_EXCEPTION:
                            raise EXIT_REPORTING_EXCEPTION
                        except:
                            pass
                    for poll_function in self.polled_events:
                        #try:
                        poll_function(self.pid,rpc_connection,self.session)
                        #except: pass
        except EXIT_REPORTING_EXCEPTION:
            pass
