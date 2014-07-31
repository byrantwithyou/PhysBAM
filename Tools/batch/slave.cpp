#include <iostream>
#include <vector>
#include "pack.h"
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/sync/named_condition.hpp>
#include <boost/interprocess/sync/named_mutex.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/program_options.hpp>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

namespace ip = boost::interprocess;
namespace po = boost::program_options;
using namespace std;

int main(int argc,char* argv[])
{
    po::options_description desc("Options");

    string name="/tmp/batch_queue_default",in="/dev/null",out="/dev/null",err="/dev/null";
    int priority=50;
    vector<int> deps;
    vector<string> tokens;
    desc.add_options()
        ("help,h","usage info")
        ("priority,p",po::value<int>(&priority)->default_value(priority),"job priority")
        ("depends,d",po::value<vector<int> >(&deps),"depends on this job id")
        ("token,t",po::value<vector<string> >(&tokens),"command to run")
        ("name,n",po::value<string>(&name)->default_value(name),"batch queue name")
        ("in,i",po::value<string>(&in)->default_value(in),"read stdin from this file")
        ("out,o",po::value<string>(&out)->default_value(out),"write stdout to this file")
        ("err,e",po::value<string>(&err)->default_value(err),"write stderr to this file")
        ("append-stdout,a","append stdout")
        ("append-stderr,A","append stderr")
        ("kill,k","kill master")
        ("cwd,c","run command from the current working directory")
        ("env,E","set environment before running command")
        ("kill-now,K","kill master immediately");

    po::positional_options_description pod;
    pod.add("token", -1);

    po::variables_map vm;
    try
    {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
    }
    catch(...)
    {
        cout<<desc<<endl;
        return 1;
    }
    po::notify(vm);

    if(vm.count("help"))
    {
        cout<<desc<<endl;
        return 1;
    }

    int sock = socket(PF_LOCAL, SOCK_STREAM, 0);
    if(sock < 0)
    {
        cout<<"failed to create socket"<<endl;
        return 1;
    }

    struct sockaddr_un sock_name;
    sock_name.sun_family = AF_LOCAL;
    strcpy(sock_name.sun_path, name.c_str());
    if(connect(sock, (sockaddr *)&sock_name, SUN_LEN(&sock_name)) < 0)
    {
        cout<<"failed to connect socket"<<endl;
        return 1;
    }

    if(vm.count("kill") || vm.count("kill-now"))
    {
        unsigned char buff[100];
        size_t k = pack(buff, vm.count("kill") ? (size_t)-1 : (size_t)-2);
        write(sock, buff, k);
        close(sock);
        return 0;
    }

    int append_out=vm.count("append-stdout")>0;
    int append_err=vm.count("append-stderr")>0;
    bool save_env=vm.count("env")>0;
    bool save_cwd=vm.count("cwd")>0 || save_env;

    const char* cwd = "";
    if(save_cwd) cwd = get_current_dir_name();

    vector<string> env;
    if(save_env) for(int i=0;environ[i];i++) env.push_back(environ[i]);

    size_t size = 0;
    size += pack_size(size);
    size += pack_size(tokens);
    size += pack_size(deps);
    size += pack_size(in);
    size += pack_size(out);
    size += pack_size(err);
    size += pack_size(append_out);
    size += pack_size(append_err);
    size += pack_size(priority);
    size += pack_size(cwd);
    size += pack_size(env);

    unsigned char * buff = new unsigned char[size];
    size_t k = 0;
    k += pack(buff + k, size);
    k += pack(buff + k, tokens);
    k += pack(buff + k, deps);
    k += pack(buff + k, in);
    k += pack(buff + k, out);
    k += pack(buff + k, err);
    k += pack(buff + k, append_out);
    k += pack(buff + k, append_err);
    k += pack(buff + k, priority);
    k += pack(buff + k, cwd);
    k += pack(buff + k, env);
    assert(k == size);

    int w = write(sock, buff, k);

    int job_id = -1;
    if(read(sock, &job_id, sizeof(job_id)) == sizeof(job_id))
        printf("%i\n", job_id);

    if(save_cwd) free((char*)cwd);

    close(sock);

    delete [] buff;

    return 0;
}
