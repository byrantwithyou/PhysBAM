#!/usr/bin/python
import CONFIG
from psim.Base import SOCKET
import os
import mutex
import time
import threading
import socket

# Host representation: dictionary having hostname, locked, max_cpus, max_memory, claims
# Claims representation dictionary user,identifier,cpus,memory

class HOSTS_SERVER:
    def __init__(self):
        self.hosts_filename="/data/psim/data/hosts.py"
        self.commands=["Host_List","Add_Host","Delete_Host","Lock_Host","Unlock_Host","Claim_Host","Release_Host"]
        self.mutex=threading.Lock()
        self.next_claim_id=1
        self.Read_Host_List()

    # Public
    def Read_Host_List(self):
        self.hosts=eval(open(self.hosts_filename).read())
        for host in self.hosts.keys():
            for claim in self.hosts[host]["claims"].keys():
                if claim>self.next_claim_id: self.next_claim_id=claim
        self.next_claim_id+=1

    def Write_Host_List(self):
        open(self.hosts_filename,"w").write(repr(self.hosts))

    # Public
    def Host_List(self):
        return self.hosts

    def Add_Host(self,host,max_cpus,max_memory):
        if self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Host already added")
        self.hosts[host]={"locked":locked,"max_cpus":max_cpus,"max_memory":max_memory,"claims":{}}
        self.Write_Host_List()
    
    def Delete_Host(self,host,max_cpus,max_memory):
        if not self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Invalid host")
        del self.hosts[host]
        self.Write_Host_List()
    
    def Lock_Host(self,host):
        if not self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Invalid host")
        if self.hosts[host]["locked"]: raise SOCKED.COMMAND_EXCEPTION("Host already locked")
        self.hosts[host]["locked"]=1
        self.Write_Host_List()
    
    def Unlock_Host(self,host):
        if not self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Invalid host")
        if not self.hosts[host]["locked"]: raise SOCKED.COMMAND_EXCEPTION("Host already unlocked")
        self.hosts[host]["locked"]=0
        self.Write_Host_List()
    
    def Claim_Host(self,host,identifier,user,cpus,memory):
        if not self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Invalid host")
        host_info=self.hosts[host]
        if host_info["locked"]: raise SOCKET.COMMAND_EXCEPTION("Machine is locked")
        claim_info=host_info["claims"]
        claimed_memory,claimed_cpus=0,0
        for claim in claim_info.keys():
            claim_data=claim_info[claim]
            claimed_memory+=claim_data["memory"]
            claimed_cpus+=claim_data["cpus"]
        resource_errors=[]
        if claimed_memory+memory > host_info["max_memory"]: resource_errors.append("memory exceeds available")
        if claimed_cpus+cpus > host_info["max_cpus"]: resource_errors.append("cpus exceeds available")
        if len(resource_errors): raise SOCKET.COMMAND_EXCEPTION("Resource exceeded: "+", ".join(resource_errors))
        claim={"identifier":identifier,"user":user,"memory":memory,"cpus":cpus}
        claim_id=self.next_claim_id
        self.hosts[host]["claims"][claim_id]=claim
        self.next_claim_id+=1
        self.Write_Host_List()
        return claim_id

    def Release_Host(self,host,claim_id):
        if not self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Invalid host")
        if not self.hosts[host]["claims"].has_key(claim_id): raise SOCKET.COMMAND_EXCEPTION("Invalid claim ID")
        del self.hosts[host]["claims"][claim_id]
        self.Write_Host_List()
    
if __name__ == "__main__":
    server=HOSTS_SERVER()
    SOCKET.SERVER(socket.gethostbyname(os.environ['HOSTNAME']),CONFIG.hosts_port,server,
                  (CONFIG.server_private_key_file,CONFIG.server_certificate_file,CONFIG.ca_certificate_file))
