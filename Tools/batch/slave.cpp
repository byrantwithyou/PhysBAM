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
        ("priority,p",po::value<int>(&priority)->default_value(priority),"job priority [0, 100]")
        ("depends,d",po::value<vector<int> >(&deps),"depends on this job id")
        ("token,t",po::value<vector<string> >(&tokens),"command to run")
        ("name,n",po::value<string>(&name)->default_value(name),"batch queue name")
        ("in,i",po::value<string>(&in)->default_value(in),"read stdin from this file")
        ("out,o",po::value<string>(&out)->default_value(out),"write stdout to this file")
        ("err,e",po::value<string>(&err)->default_value(err),"write stderr to this file")
        ("append-stdout,a","append stdout")
        ("append-stderr,A","append stderr")
        ("kill,k","kill master")
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

    size_t size = 2 * sizeof(size);
    for(size_t i = 0; i < tokens.size(); i++)
        size += sizeof(size_t) + tokens[i].size();
    size += sizeof(size_t) + sizeof(deps[0]) * deps.size();
    size += sizeof(size_t) + in.size();
    size += sizeof(size_t) + out.size();
    size += sizeof(size_t) + err.size();
    size += sizeof(append_out);
    size += sizeof(append_err);
    size += sizeof(priority);

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
    assert(k == size);

    int w = write(sock, buff, k);

    int job_id = -1;
    if(read(sock, &job_id, sizeof(job_id)) == sizeof(job_id))
        printf("%i\n", job_id);

    close(sock);

    delete [] buff;

    return 0;
}
