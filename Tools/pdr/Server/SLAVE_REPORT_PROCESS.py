import os
from RENDER_COMMAND import RENDER_COMMAND

class SLAVE_REPORT_PROCESS(RENDER_COMMAND):
    def __init__(self,server):
        RENDER_COMMAND.__init__(self,"Slave_Report_Process",server)

    def Process(self,submit_id):
        pass
        
