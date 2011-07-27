import os
from RENDER_COMMAND import RENDER_COMMAND

class REMOVE_SLAVE_COMMAND(RENDER_COMMAND):
    def __init__(self,server):
        RENDER_COMMAND.__init__(self,"Remove_Slave",server)

    def Process(self,submit_id):
        slavefile=os.path.basename(submit_id)
        delete_host,delete_id=slavefile.split("__")
        for i in os.listdir(self.server.slaves_directory):
            platform,priority,memory,host,id=i.split("__")
            if delete_host==host and delete_id==id:
                os.unlink(os.path.join(self.server.slaves_directory,i))
                
        
