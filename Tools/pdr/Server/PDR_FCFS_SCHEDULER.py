#!/usr/bin/python
# This class schedules jobs first come first serve, first coming is defined as job id because
# job ids are monotonically increasing

class PDR_FCFS_SCHEDULER:
    def Schedule(self,pending_jobs,free_slaves,maximum_number_of_jobs):
        class FINISHED_SCHEDULING(Exception):pass
        pending_jobs.sort(lambda x,y:int(x[0])-int(y[0]))
#        print pending_jobs
        scheduled_pairs=[]
        try:
            for job,priority,mem_req,frames in pending_jobs:
                frames.sort(lambda x,y:int(x)-int(y))  # want to do frames in order
                for frame in frames:
                    if len(free_slaves)<=0 or len(scheduled_pairs)>=maximum_number_of_jobs:
                        raise FINISHED_SCHEDULING
                    scheduled_pairs.append((job,frame,free_slaves.pop()[0]))
        except FINISHED_SCHEDULING,e: pass
#        print scheduled_pairs
        return scheduled_pairs
      
# Test code
if __name__=="__main__":
    sched=PDR_FCFS_SCHEDULER()
    
    pending=[(1,2,128,[1,2,3]),(2,3,1024,[10,20,30,40])]
    free=[("server_1",2,2000),("server_2",2,2000),("server_3",2,500),("server_4",2,500),("server_5",3,4000),("server_6",3,1024),("server_7",3,500)]
    print repr(sched.Schedule(pending,free,10))

            
            
