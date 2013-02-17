#include <iostream>
#include <vector>
#include "pack.h"
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/program_options.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>

namespace ip = boost::interprocess;
namespace po = boost::program_options;
using namespace std;

int main(int argc,char* argv[])
{
    po::options_description desc("Options");

    string name="batch_queue_default",in="/dev/null",out="/dev/null",err="/dev/null";
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
        ("print-job-id,j","prints job id");

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

    ip::shared_memory_object shm(ip::open_only, name.c_str(), ip::read_write);
    ip::mapped_region region(shm, ip::read_write);
    memory_layout* layout = reinterpret_cast<memory_layout*>(region.get_address());

    if(vm.count("kill"))
    {
        try
        {
            {
                ip::scoped_lock<ip::interprocess_mutex> lock(layout->mutex);
                if(layout->filled) layout->cond_no_job.wait(lock);
                layout->filled = 1;
                layout->message_length = 0;
                layout->cond_job.notify_one();
                layout->next_job_id = -1;
            }
        }
        catch(...)
        {
            printf("failed to kill master\n");
            return 3;
        }
        return 0;
    }

    size_t size = (tokens.size()+2)*sizeof(size_t) + sizeof(int)*(deps.size()+2);
    for(size_t i = 0; i < tokens.size(); i++)
        size += tokens[i].size();
    size += sizeof(size_t) + in.size();
    size += sizeof(size_t) + out.size();
    size += sizeof(size_t) + err.size();

    if(size > layout->message_max_length)
    {
        printf("maximum job size: %zi.  This job size: %zi\n", layout->message_max_length, size);
        return 2;
    }
    int job_id = -1;

    int append_out=vm.count("append-stdout")>0;
    int append_err=vm.count("append-stderr")>0;

    unsigned char * buff = new unsigned char[size];
    size_t k = 0;
    k += pack(buff + k, tokens);
    k += pack(buff + k, deps);
    k += pack(buff + k, in);
    k += pack(buff + k, out);
    k += pack(buff + k, err);
    k += pack(buff + k, append_out);
    k += pack(buff + k, append_err);
    assert(k == size);

    try
    {
        {
            ip::scoped_lock<ip::interprocess_mutex> lock(layout->mutex);
            if(layout->filled) layout->cond_no_job.wait(lock);
            layout->filled = 1;
            layout->message_length = 0;
            layout->cond_job.notify_one();

            job_id = layout->next_job_id;
            layout->priority = priority;
            memcpy((unsigned char*)region.get_address() + layout->message_offset, buff, size);
            layout->message_length = size;
        }

        if(vm.count("print-job-id"))
            printf("%i\n", job_id);

    }
    catch(...)
    {
        printf("failed to send job\n");
        return 3;
    }

    delete [] buff;

    return 0;
}
