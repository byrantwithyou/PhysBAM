import sys
from physbam import *



tri_filename=None
try:
    exe,tri_filename=sys.argv
except:
    print "Usage: %s <tri filename>"%sys.argv[0]
    sys.exit(1)


for i in ["Initialize_Incident_Old","Initialize_Incident_New","Initialize_Incident_Hash"]:
    scope=LOG_SCOPE(str(i))
    LOG.Time("Reading")
    tri=TRIANGULATED_SURFACE_f.Create()
    Read_From_File("float",tri_filename,tri)
    globals()[i](tri.mesh)
    del scope

LOG.Finish_Logging()


