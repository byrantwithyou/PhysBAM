import os
from RENDER_COMMAND import RENDER_COMMAND

class ADD_SLAVE_COMMAND(RENDER_COMMAND):
    def __init__(self,server):
        RENDER_COMMAND.__init__(self,"Add_Slave",server)

    def Process(self,submit_id):
        slave_filename=os.path.basename(submit_id)
        platform,priority,memory,host,instance=slave_filename.split("__")
        if platform != "Linux" and platform != "Windows":
            raise NameError,("Slave was unacceptable platform '%s'"%platform)
        open(os.path.join(self.server.slaves_directory,slave_filename),"w")
        
