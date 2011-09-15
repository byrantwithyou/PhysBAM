import sys,os
import signal
print "Hi"

def killhandler(signum,frame):
    print "Got signal"
    print "Killing pid=%d"%pid
    os.kill(pid,signal.SIGTERM)
    print "Exiting..."
    sys.exit(0)
signal.signal(signal.SIGTERM, killhandler)
    

pid=os.fork()
if pid>0:
    print "Forked process pid=%d\n"%pid
    killed_pid,exitcode=os.waitpid(pid,0)
    print "Quit with exitcode "+str(exitcode)
else:
    err_log=out_log=file('out.txt','a+')
    dev_null=file('/dev/null','r')
    os.dup2(out_log.fileno(),sys.stdout.fileno())
    os.dup2(err_log.fileno(),sys.stderr.fileno())
    os.dup2(dev_null.fileno(),sys.stdin.fileno())
    os.execvp("watch",["watch","w"])

