#!/usr/bin/python
import threading
import sys
import pickle
import socket
import os
import traceback
from OpenSSL import SSL
from OpenSSL import crypto

def daemonize():
    sys.stdin=None
    log=open("log.txt","w")
    sys.stdout=log
    sys.stderr=log
    PID = os.fork()
    if PID != 0:
        print "shutting down parent process:",PID
        sys.exit()
    os.setsid()
    PID = os.fork()
    if PID != 0:
        print "shutting down parent process:",PID
        sys.exit()

class CONNECTION_LOST_ERROR(Exception):
    def __init__(self,value=""):
        self.value=str(value)
        
    def __str__(self):
        return self.value

class COMMAND_EXCEPTION(Exception):
    def __init__(self,value,detail=""):
        self.value=value
        self.detail=detail

    def __str__(self):
        return "%s %s"%(self.value,self.detail)

class SOCKET_READER:
    def __init__(self,delimiter="\n"):
        #self.socket=socket
        self.buffer=""
        self.delimiter=delimiter
        
    def get(self,sock):
        while 1:
            index=self.buffer.find(self.delimiter)
            if index != -1:
                line=self.buffer[:index+len(self.delimiter)]
                self.buffer=self.buffer[index+len(self.delimiter):]
                return line
            else:
                try:
                    chunk=sock.recv(1024)
                except (SSL.WantReadError, SSL.WantWriteError, SSL.WantX509LookupError):
                    pass
                except SSL.ZeroReturnError:
                    raise CONNECTION_LOST_ERROR()
                except SSL.Error, errors:
                    raise CONNECTION_LOST_ERROR(errors)
                else:
                    if chunk=="": return ""
                    self.buffer+=chunk


class SERVER_THREAD(threading.Thread):
    def __init__(self,client_socket,client_host,commands,mutex):
        threading.Thread.__init__(self)
        self.socket=client_socket
        self.host=client_host
        self.commands=commands
        self.mutex=mutex
        self.reader=SOCKET_READER(".\n")
        self.start()


    def run(self):
        while 1:
            # Get the python command object
            pickle_chunk=0
            try:
                pickle_chunk=self.reader.get(self.socket)
            except CONNECTION_LOST_ERROR:
                break
            #print "pickle_chunk=%s"%pickle_chunk
            if pickle_chunk=='': break
            thing=None
            result=(None,None) # (exception, result)
            try:
                try: thing=pickle.loads(pickle_chunk[:-1])
                except: raise COMMAND_EXCEPTION("Non-pickled encountered from %s"%repr(self.host))
                # make sure command is of form (command, args)
                command,args=None,None
                try:
                    command,args=thing
                    print "Client %s:%d Command %s%s"%(self.host[0],self.host[1],command,repr(args))
                except: raise COMMAND_EXCEPTION("Invalid command packing")
                # Make sure command exists
                if not self.commands.has_key(command): raise COMMAND_EXCEPTION("Command \"%s\" not found"%repr(command))
                try:
                    self.mutex.acquire()
                    result=(None,self.commands[command](*args))
                finally:
                    self.mutex.release()
            except COMMAND_EXCEPTION,e:
                result=((e.value,e.detail),None)
            except:
                result=(("Other Exception",str(sys.exc_info()[1])+"\n"+traceback.format_exc()),None)

            picked_result=pickle.dumps(result)
            self.socket.sendall(picked_result+"\n")                
        print "Client %s disconnected"%repr(self.host)
        self.socket.shutdown()
        self.socket.close()

# check to make sure certificate is good
# NOTE THIS IS OUTSIDE THE CLASSES otherwise circular reference count hell happens
def verify_client_certificate(connection,certificate,errnum,depth,ok):
    print "Got certificate: %s ok=%d"%(certificate.get_subject(),ok)
    return ok

class SERVER:
    def __init__(self,bindhost,port,rpc_server,secure=None):
        self.commands=dict(map(lambda x: (x,getattr(rpc_server,x)),rpc_server.commands))
        self.mutex=rpc_server.mutex

        self.raw_listen_socket=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        # Make listen socket
        if secure:
            server_pkey_file,server_cert_file,ca_cert_file=secure
            ctx=SSL.Context(SSL.SSLv23_METHOD)
            ctx.set_verify(SSL.VERIFY_PEER,verify_client_certificate)
            ctx.use_privatekey_file(server_pkey_file)
            ctx.use_certificate_file(server_cert_file)
            ctx.load_verify_locations(ca_cert_file)
            self.listen_socket=SSL.Connection(ctx,self.raw_listen_socket)
        else:
            self.listen_socket=self.raw_listen_socket
        self.listen_socket.bind((bindhost,port))
        self.listen_socket.listen(5)

        # Server loop
        try:
            while 1:
                (client_socket,addr)=self.listen_socket.accept()
                print "Connection on %s"%repr(addr)
                SERVER_THREAD(client_socket,addr,self.commands,self.mutex)
        finally:
            self.listen_socket.close()


class CLIENT:
    # usage is (hostname, port, (privatekey_text,certificate_text,ca_certificate_text))
    def __init__(self,hostname,port,secure=None):
        # Initialize context
        if secure:
            server_ptext,server_cert_text,ca_cert_text=secure
            server_pkey=crypto.load_privatekey(crypto.FILETYPE_PEM,server_ptext)
            server_cert=crypto.load_certificate(crypto.FILETYPE_PEM,server_cert_text)
            #ca_cert=crypto.load_certificate(crypto.FILETYPE_PEM,ca_cert_text)
            ctx = SSL.Context(SSL.SSLv23_METHOD)
            ctx.set_verify(SSL.VERIFY_PEER,verify_client_certificate) # Demand a certificate
            ctx.use_privatekey(server_pkey)
            ctx.use_certificate(server_cert)
            #ctx.use_verify_locations(ca_cert)
            #ctx.use_privatekey_file("/n/curvature/data/psim/keys/client.pkey")
            #ctx.use_certificate_file("/n/curvature/data/psim/keys/client.cert")
            ctx.load_verify_locations("/n/curvature/data/psim/keys/CA.cert")
            self.socket=SSL.Connection(ctx,socket.socket(socket.AF_INET,socket.SOCK_STREAM))
        else:
            self.socket=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        self.socket.connect((hostname,port))
        self.reader=SOCKET_READER(".\n")

    def Function_Call(self,command,args):
        #return 
        pickled_command=pickle.dumps((command,args))
        self.socket.sendall(pickled_command+"\n")
        #print "sent command"
        pickled_data=self.reader.get(self.socket)
        #print "in call(%d) %s"%(os.getpid(),pickled_data)
        exception,value=None,None
        try:
            data=pickle.loads(pickled_data)
            #print repr(data)
            exception,value=data
        except:
            raise COMMAND_EXCEPTION("Failed to read return value")
        if exception!=None:
            raise COMMAND_EXCEPTION(exception[0],exception[1])
            try:
                error,error_detail=exception
            except COMMAND_EXCEPTION:
                raise
            raise COMMAND_EXCEPTION(exception[0],exception[1])

        #print "got something"
        if len(self.reader.buffer) != 0: print "WEIRD! Data left in buffer."
        return value

    def __getattr__(self,method_name):
        return lambda *x: self.Function_Call(method_name,x)
    
    def __del__(self):
        print "Destructing SOCKET.client"
