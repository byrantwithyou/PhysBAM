#include <iostream>
#include <map>
#include <vector>
#include "pack.h"
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/thread/mutex.hpp>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <signal.h>

namespace ip = boost::interprocess;
namespace po = boost::program_options;

static const size_t buffer_size = 1<<16;

struct job
{
    int id;
    int priority;
    vector<string> argv_data;
    string in,out,err;
    int append_out,append_err;
};

multimap<int,job*> jobs;
map<int,vector<int> > dependents;
map<int,pair<int,job*> > waiting_jobs;
set<int> finished_jobs;
int first_uncompleted_job;
set<pid_t> pids;

boost::condition_variable local_cond_job;
boost::mutex local_mutex;
int num_waiting_threads;
int num_threads_alive;
int die;
memory_layout* layout;

void add_job(job* j)
{
    cout<<"Add:";
    for(size_t i = 0; i < j->argv_data.size(); i++)
        cout<<" "<<j->argv_data[i];
    cout<<"   "<<j->out;
    cout<<endl;

    if(num_waiting_threads)
    {
        num_waiting_threads--;
        local_cond_job.notify_one();
    }
    jobs.insert(pair<int,job*>(j->priority,j));
}

void run_jobs()
{
    job * j = 0;
    int last_j_id = -1;
    {
        boost::unique_lock<boost::mutex> lock(local_mutex);
        num_threads_alive++;
    }
    while(!die)
    {
        boost::unique_lock<boost::mutex> lock(local_mutex);
        if(last_j_id >= 0)
        {
            finished_jobs.insert(last_j_id);
            while(!finished_jobs.empty() && *finished_jobs.begin()==first_uncompleted_job)
            {
                finished_jobs.erase(finished_jobs.begin());
                first_uncompleted_job++;
            }

            map<int,vector<int> >::iterator it = dependents.find(last_j_id);
            if(it != dependents.end())
            {
                const vector<int>& v = it->second;
                for(size_t i = 0; i < v.size(); i++)
                {
                    map<int,pair<int,job*> >::iterator it2 = waiting_jobs.find(v[i]);
                    if(it2 == waiting_jobs.end()) continue;
                    if(!--it2->second.first)
                    {
                        add_job(it2->second.second);
                        waiting_jobs.erase(it2);
                    }
                }
                dependents.erase(it);
            }
            last_j_id = -1;
        }

        if(jobs.empty())
        {
            num_waiting_threads++;
            local_cond_job.wait(lock);
        }

        multimap<int,job*>::iterator last=jobs.end();
        last--;
        j = last->second;
        jobs.erase(last);
        lock.unlock();

        if(!j)
        {
            cout<<"kill "<<jobs.size()<<endl;
            break;
        }

        last_j_id = j->id;

        vector<char*> argv;
        for(size_t i = 0; i < j->argv_data.size(); i++)
            argv.push_back(const_cast<char*>(j->argv_data[i].c_str()));
        argv.push_back(0);

        cout<<"Run:";
        for(size_t i = 0; i < j->argv_data.size(); i++)
            cout<<" "<<j->argv_data[i];
        cout<<endl;
        pid_t pid = fork();
        if(!pid)
        {
            close(0);
            close(1);
            close(2);
            int in = open(j->in.c_str(),O_RDONLY);
            int out = open(j->out.c_str(),O_WRONLY|O_CREAT|(j->append_out?O_APPEND:0),0755);
            int err = open(j->err.c_str(),O_WRONLY|O_CREAT|(j->append_err?O_APPEND:0),0755);
            if(in<0 || out<0 || err<0) exit(1);
            execvp(argv[0], &argv[0]);
            exit(1);
        }

        {
            boost::unique_lock<boost::mutex> lock(local_mutex);
            pids.insert(pid);
        }
        
        waitpid(pid, 0, 0);

        {
            boost::unique_lock<boost::mutex> lock(local_mutex);
            pids.erase(pid);
        }

        delete j;
    }

    {
        boost::unique_lock<boost::mutex> lock(local_mutex);
        num_threads_alive--;
        if(!num_threads_alive) 
        {
            die = 1;
            layout->cond_job.notify_one();
        }
    }
}

int main(int argc,char* argv[])
{
    po::options_description desc("Options");

    string name="batch_queue_default",pid_file="";
    int num_processes = boost::thread::hardware_concurrency();
    vector<int> deps;
    vector<string> tokens;
    desc.add_options()
        ("help,h","usage info")
        ("processes,p",po::value<int>(&num_processes)->default_value(num_processes),"number of processes to run")
        ("name,n",po::value<string>(&name)->default_value(name),"batch queue name")
        ("pid,P",po::value<string>(&pid_file)->default_value(name),"write pid to this file")
        ("background,b","put seld in the background when it is safe");

    po::variables_map vm;
    try
    {
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
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

    int next_job_id = -1;
    ip::shared_memory_object::remove(name.c_str());
    ip::shared_memory_object shm(ip::create_only, name.c_str(), ip::read_write);
    shm.truncate(buffer_size);
    ip::mapped_region region(shm, ip::read_write);
    memset(region.get_address(), 0, region.get_size());
    layout = new (region.get_address()) memory_layout;
    layout->message_offset = sizeof(memory_layout);
    size_t message_max_length = region.get_size() - layout->message_offset;
    layout->message_max_length = message_max_length;
    layout->message_length = 0;
    layout->next_job_id = ++next_job_id;
    layout->priority = 0;
    unsigned char *message = (unsigned char*)region.get_address() + layout->message_offset;

    if(vm.count("background"))
    {
        pid_t pid = fork();
        if(pid)
        {
            if(vm.count("pid"))
            {
                ofstream fout(pid_file.c_str());
                fout<<pid;
            }
            return 0;
        }
    }
    else if(vm.count("pid"))
    {
        ofstream fout(pid_file.c_str());
        fout<<getpid();
    }

    vector<boost::thread*> threads;
    for(int i = 0; i < num_processes; i++)
        threads.push_back(new boost::thread(run_jobs));

    {
        ip::scoped_lock<ip::interprocess_mutex> lock(layout->mutex);
        layout->cond_no_job.notify_one();
    }

    unsigned char buffer[buffer_size];

    while(!die)
    {
        cout<<"loop"<<endl;
        int job_id = -1, priority = 0;
        size_t message_size = 0;
        try
        {
            ip::scoped_lock<ip::interprocess_mutex> lock(layout->mutex);
            if(!layout->filled) layout->cond_job.wait(lock);
            if(die) break;
            layout->filled = 0;
            layout->cond_no_job.notify_one();

            cout<<"job id "<<layout->next_job_id<<std::endl;
            if(layout->next_job_id == -2)
            {
                die = 1;
                break;
            }

            if(layout->message_length <= message_max_length)
            {
                job_id = next_job_id;
                priority = layout->priority;
                memcpy(buffer, message, layout->message_length);
                message_size = layout->message_length;
            }

            layout->message_length = 0;
            layout->message_offset = sizeof(memory_layout);
            layout->message_max_length = message_max_length;
            layout->next_job_id = ++next_job_id;
            layout->priority = 0;
        }
        catch(...)
        {
            continue;
        }

        job* j = new job;
        j->id = job_id;
        j->priority = priority;
        vector<int> deps;
        try
        {
            int k = 0;
            k += unpack(buffer + k, message_size - k, j->argv_data);
            k += unpack(buffer + k, message_size - k, deps);
            k += unpack(buffer + k, message_size - k, j->in);
            k += unpack(buffer + k, message_size - k, j->out);
            k += unpack(buffer + k, message_size - k, j->err);
            k += unpack(buffer + k, message_size - k, j->append_out);
            k += unpack(buffer + k, message_size - k, j->append_err);

            if(k != message_size) cout<<"only used "<<k<<" of "<<message_size<<std::endl;
        }
        catch(...)
        {
            delete j;
            continue;
        }

        try
        {
            boost::unique_lock<boost::mutex> lock(local_mutex);

            if(job_id < 0)
            {
                delete j;
                while(num_waiting_threads--) local_cond_job.notify_one();
                for(size_t i = 0; i < threads.size(); i++)
                    jobs.insert(pair<int,job*>(INT_MIN,0));
                continue;
            }

            size_t num_deps = deps.size();
            for(size_t i = 0; i < deps.size(); i++)
            {
                if(deps[i] >= first_uncompleted_job && deps[i] < job_id && finished_jobs.find(deps[i])==finished_jobs.end())
                    dependents[deps[i]].push_back(job_id);
                else num_deps--;
            }

            if(num_deps) waiting_jobs[job_id] = pair<int,job*>(num_deps, j);
            else add_job(j);
        }
        catch(...)
        {
        }
    }

    {
        boost::unique_lock<boost::mutex> lock(local_mutex);
        for(set<pid_t>::iterator it = pids.begin(); it != pids.end(); it++)
            kill(*it, SIGTERM);
        pids.clear();

        while(num_waiting_threads--) local_cond_job.notify_one();
        for(size_t i = 0; i < threads.size(); i++)
            jobs.insert(pair<int,job*>(INT_MIN,0));
    }

    try
    {
        boost::unique_lock<boost::mutex> lock(local_mutex);

        while(num_waiting_threads--) local_cond_job.notify_one();
        for(size_t i = 0; i < threads.size(); i++)
            jobs.insert(pair<int,job*>(INT_MIN,0));
    }
    catch(...)
    {
    }

    for(size_t i = 0; i < threads.size(); i++)
    {
        threads[i]->join();
        delete threads[i];
    }

    layout->~memory_layout();
    ip::shared_memory_object::remove(name.c_str());

    return 0;
}
