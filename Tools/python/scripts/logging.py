import physbam

# Make example and driver and run
physbam.LOG.Initialize_Logging(False,False,99999,False)
physbam.LOG.Print("FUCK")
for i in range(2):
    scope=physbam.LOG_SCOPE("tester %d"%i)
    physbam.LOG.Stat("tester",100)
    for j in range(4):
        interior_scope=physbam.LOG_SCOPE("interior %d"%j)
        physbam.LOG.Print("I WEIRD\n")
        physbam.LOG.Stat("ab",100.212)
        interior_scope=None
    physbam.LOG.Print("HI")
    scope=None
physbam.LOG.Finish_Logging()
