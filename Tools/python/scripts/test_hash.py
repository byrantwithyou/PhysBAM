import physbam
import random

pb_hash=physbam.HASHTABLE_i_i()
python_hash={}


r=random.Random()
operation=0

def FATAL_ERROR(x):
    print "FATAL ERROR %s"%x
    sys.exit(1)


def display():
    #for i in pb_hash:
    #    print "         %-15s %20d %20d"%(i.state,i.key,i.data)
    #print "           Entries %d Next Resize %d  Max Size %d Load Factor %f"%(pb_hash.Size(),pb_hash.Next_Resize(),pb_hash.Max_Size(),float(pb_hash.Size())/float(pb_hash.Max_Size()))
    pass

def insert():
    key=0
    while 1:
        key=random.randint(0,10000000)
        if not python_hash.has_key(key): break
    #print "    key=%d value=%d"%(key,operation)
    if key in pb_hash: FATAL_ERROR("python does not have but physbam does")
    pb_hash[key]=operation
    python_hash[key]=operation
    if python_hash[key]!=pb_hash[key]: FATAL_ERROR("dont match after insert")
    display()
        
def delete():
    #print "DELETE"
    if len(python_hash)<=0: return
    key=r.choice(python_hash.keys())
    if not key in pb_hash: FATAL_ERROR("key not in pb on delete")
    if python_hash[key]!=pb_hash[key]: FATAL_ERROR("keys dont match before delete")
    del python_hash[key]
    del pb_hash[key]
    if key in pb_hash: FATAL_ERROR("key in pb after delete")
    

def update():
    if len(python_hash)<=0: return
    key=r.choice(python_hash.keys())
    if not key in pb_hash: FATAL_ERROR("key not in pb on delete")
    if python_hash[key]!=pb_hash[key]: FATAL_ERROR("keys dont match before update")
    python_hash[key]=operation
    pb_hash[key]=operation
    if python_hash[key]!=pb_hash[key]: FATAL_ERROR("keys dont match after update")

operations=[insert,insert,delete,update]

for i in range(1000000):
    operation+=1
    r.choice(operations)()
    if i%100000:
        print "Did %d operations"%operation


    
#print pb_hash
print dir(pb_hash)
    
