import sys
import os
import os.path
import traceback
import PDR_UTILITY

class RENDER_COMMAND:
    def __init__(self,name,server):
        self.name=name
        self.server=server
        self.directory=os.path.join(self.server.server_root,"Commands",self.name)

    def Handle_Pending_Commands(self):
        for i in os.listdir(self.directory):
            self.server.Log(self.server.command_log,"[%-20s] '%s'" % (self.name,i))
            try:
                os.unlink(os.path.join(self.directory,i))
            except:
                print "<----- Unlink of file '%s' failed ----->" % i
            try:
                self.Process(i)
            except:
                print "<----- Command '%s' failed with '%s' ----->" % (self.name,i)
                traceback.print_exc(file=sys.stdout)

    def Initialize_Directory_Structure(self,server):
        PDR_UTILITY.Make_Directory_If_Does_Not_Exist(self.directory)
        
    
