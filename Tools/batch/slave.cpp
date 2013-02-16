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

    string name="batch_queue_default";
    int priority=50;
    vector<int> deps;
    vector<string> tokens;
    desc.add_options()
        ("help,h","usage info")
        ("priority,p",po::value<int>(&priority)->default_value(priority),"job priority [0, 100]")
        ("depends,d",po::value<vector<int> >(&deps),"depends on this job id")
        ("token,t",po::value<vector<string> >(&tokens),"command to run")
        ("name,n",po::value<string>(&name)->default_value(name),"batch queue name")
        ("kill,k","kill master")
        ("print-job-id,i","prints job id");

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

    size_t size = (tokens.size()+2)*sizeof(size_t) + sizeof(int)*deps.size();
    for(size_t i = 0; i < tokens.size(); i++)
        size += tokens[i].size();

    if(size > layout->message_max_length)
    {
        printf("maximum job size: %zi.  This job size: %zi\n", layout->message_max_length, size);
        return 2;
    }
    int job_id = -1;

    unsigned char * buff = new unsigned char[size];
    size_t k = 0;
    k += pack(buff + k, tokens);
    k += pack(buff + k, deps);
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
