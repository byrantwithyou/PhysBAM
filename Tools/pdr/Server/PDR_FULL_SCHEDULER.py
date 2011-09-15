#!/usr/bin/python
# This class schedules jobs first come first serve, first coming is defined as job id because
# job ids are monotonically increasing
# It also respects memory and priority restrictions.
# Higher priority jobs starve lower priority ones

# format for pending jobs list:
#    (job_id,priority,mem_req,[frames])
# format for free slaves list:
#    (slave_name,min_priority,mem_available)
# format for return list:
#    (job_id,frame,slave_name)

class PDR_FULL_SCHEDULER:
    
    def Schedule(self,pending_jobs,free_slaves,maximum_number_of_jobs):
        class FINISHED_SCHEDULING(Exception):pass
        class NEXT_JOB(Exception):pass
        def Sort_Job_Criteria(x,y):
            if x[1]==y[1]:
                return x[0]-y[0]  # forward sort on job_id in the event of a priority tie
            return -(x[1]-y[1])   # reverse sort on priorities
        scheduled_pairs=[]
        def Sort_Slaves_Criteria(x,y):
            if x[2]==y[2]:
                return -(x[1]-y[1]) # reverse sort on priority in the event of a RAM tie
            return x[2]-y[2]        # forward sort on RAM
        scheduled_pairs=[]

        # sort the machines based on memory available (least first)
        free_slaves.sort(Sort_Slaves_Criteria)        
        print "Scheduling on the following slaves:\n "+str(free_slaves)

        # sort based on priorities, then job id
        pending_jobs.sort(Sort_Job_Criteria)
        print "Scheduling in this order:\n "+str(pending_jobs)
        try:
            for job_id,priority,mem_req,frames in pending_jobs:
                print "Scheduling job %d"%job_id
                try:
                    # sort the frames so that they are done in order                
                    frames.sort(lambda x,y:int(x)-int(y))
                    for frame in frames:
                        print "  frame %d"%frame
                        # no more slaves
                        if len(free_slaves)<=0 or len(scheduled_pairs)>=maximum_number_of_jobs:
                            raise FINISHED_SCHEDULING
                        # try to find a slave that matches the requirements
                        selected_slave_name=""
                        for slave_name,slave_priority,slave_memory in free_slaves:
                            print "    Trying slave: %s"%slave_name
                            if slave_priority<=priority:  # job priority must be high enough
                                if slave_memory>=mem_req: # slave must have enough memory
                                    selected_slave_name=slave_name
                                    free_slaves.remove((slave_name,slave_priority,slave_memory))
                                    break
                        if selected_slave_name=="": # no slave was found
                            print "  No slave match found"
                            raise NEXT_JOB # may be able to schedule other jobs
                        print "  Scheduling %d frame %d on %s"%(job_id,frame,selected_slave_name)
                        scheduled_pairs.append((job_id,frame,selected_slave_name))
                except NEXT_JOB,e: pass
        except FINISHED_SCHEDULING,e: pass
        return scheduled_pairs



# Test code
if __name__=="__main__":
    sched=PDR_FULL_SCHEDULER()
    
    pending=[(1,2,128,[1,2,3]),(2,3,1024,[10,20,30,40])]
    free=[("server_1",2,2000),("server_2",2,2000),("server_3",2,500),("server_4",2,500),("server_5",3,4000),("server_6",3,1024),("server_7",3,500)]
    print "Schedule is"
    print repr(sched.Schedule(pending,free,10))
            
            
