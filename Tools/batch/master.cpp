#include <iostream>
#include <map>
#include <vector>
#include "pack.h"
#include <boost/program_options.hpp>
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <signal.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/un.h>
#include <sys/wait.h>
#include <unistd.h>

namespace po = boost::program_options;

static const size_t buffer_size = 1<<16;

struct job
{
    int id;
    int priority;
    vector<string> argv_data;
    string in,out,err;
    int append_out,append_err;
    string cwd;
    vector<string> env_data;
};

multimap<int,job*> jobs;
map<int,vector<int> > dependents;
map<int,pair<int,job*> > waiting_jobs;
set<int> finished_jobs;
int first_uncompleted_job;
set<pid_t> pids;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond_job = PTHREAD_COND_INITIALIZER;

int num_threads_alive;
int die;
vector<pthread_t> threads;
int sock;

void add_job(job* j)
{
    pthread_cond_signal(&cond_job);
    jobs.insert(pair<int,job*>(j->priority,j));
}

void add_kill_jobs()
{
    pthread_cond_broadcast(&cond_job);
    for(size_t i = 0; i < threads.size(); i++)
        jobs.insert(pair<int,job*>(INT_MIN,0));
}

int get_number_cpus()
{
    return sysconf(_SC_NPROCESSORS_ONLN);
}

void* run_jobs(void*)
{
    job * j = 0;
    int last_j_id = -1;
    {
        pthread_mutex_lock(&mutex);
        num_threads_alive++;
        pthread_mutex_unlock(&mutex);
    }
    while(!die)
    {
        pthread_mutex_lock(&mutex);
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

        while(jobs.empty())
            pthread_cond_wait(&cond_job,&mutex);

        multimap<int,job*>::iterator last=jobs.end();
        last--;
        j = last->second;
        jobs.erase(last);
        pthread_mutex_unlock(&mutex);

        if(!j)
        {
            break;
        }

        last_j_id = j->id;

        vector<char*> argv;
        for(size_t i = 0; i < j->argv_data.size(); i++)
            argv.push_back(const_cast<char*>(j->argv_data[i].c_str()));
        std::cout<<"command: ";
        for(size_t i = 0; i < argv.size(); i++) std::cout<<argv[i]<<" ";
        std::cout<<std::endl;
        argv.push_back(0);
        vector<char*> env;
        for(size_t i=0;i<j->env_data.size();i++)
            env.push_back(const_cast<char*>(j->env_data[i].c_str()));
        env.push_back(0);

        pid_t pid = fork();
        if(!pid)
        {
            if(j->cwd.size()) chdir(j->cwd.c_str());
            close(0);
            close(1);
            close(2);
            int in = open(j->in.c_str(),O_RDONLY);
            int out = open(j->out.c_str(),O_WRONLY|O_CREAT|(j->append_out?O_APPEND:0),0755);
            int err = open(j->err.c_str(),O_WRONLY|O_CREAT|(j->append_err?O_APPEND:0),0755);
            if(in<0 || out<0 || err<0) exit(1);
            if(env.size()>1) execvpe(argv[0], &argv[0], &env[0]);
            else execvp(argv[0], &argv[0]);
            exit(1);
        }

        pthread_mutex_lock(&mutex);
        pids.insert(pid);
        pthread_mutex_unlock(&mutex);

        waitpid(pid, 0, 0);

        pthread_mutex_lock(&mutex);
        pids.erase(pid);
        pthread_mutex_unlock(&mutex);

        delete j;
    }

    pthread_mutex_lock(&mutex);
    num_threads_alive--;
    if(!num_threads_alive) die = 1;
    pthread_mutex_unlock(&mutex);
    return 0;
}

void* listen_for_jobs(void*)
{
    std::vector<unsigned char> buffer(buffer_size);

    int next_job_id = 0;
    while(!die)
    {
        sockaddr_un addr;
        socklen_t addrlen = sizeof(addr);
        int new_sock = accept(sock, (sockaddr *)&addr, &addrlen);
        if(new_sock<0) continue;

        size_t message_size;
        if(read(new_sock, &message_size, sizeof(size_t)) == sizeof(size_t))
        {
            if(message_size == (size_t)-1)
            {
                pthread_mutex_lock(&mutex);
                add_kill_jobs();
                pthread_mutex_unlock(&mutex);
            }
            else if(message_size == (size_t)-2)
            {
                pthread_mutex_lock(&mutex);
                die = 1;
                for(set<pid_t>::iterator it = pids.begin(); it != pids.end(); it++)
                    kill(*it, SIGTERM);
                pids.clear();
                add_kill_jobs();
                pthread_mutex_unlock(&mutex);
            }
            else
            {
                message_size -= sizeof(size_t);
                if(buffer.size() < message_size)
                    buffer.resize(message_size);

                size_t got = 0;
                while(got < message_size)
                {
                    ssize_t r = read(new_sock, &buffer[0], message_size - got);
                    if(r <= 0) break;
                    got += r;
                }
                if(got == message_size)
                {
                    job* j = new job;
                    j->id = next_job_id++;
                    vector<int> deps;
                    int k = 0;
                    k += unpack(&buffer[k], message_size - k, j->argv_data);
                    k += unpack(&buffer[k], message_size - k, deps);
                    k += unpack(&buffer[k], message_size - k, j->in);
                    k += unpack(&buffer[k], message_size - k, j->out);
                    k += unpack(&buffer[k], message_size - k, j->err);
                    k += unpack(&buffer[k], message_size - k, j->append_out);
                    k += unpack(&buffer[k], message_size - k, j->append_err);
                    k += unpack(&buffer[k], message_size - k, j->priority);
                    k += unpack(&buffer[k], message_size - k, j->cwd);
                    k += unpack(&buffer[k], message_size - k, j->env_data);
                    if(k != message_size) cout<<"only used "<<k<<" of "<<message_size<<std::endl;

                    pthread_mutex_lock(&mutex);
                    size_t num_deps = deps.size();
                    for(size_t i = 0; i < deps.size(); i++)
                    {
                        if(deps[i] >= first_uncompleted_job && deps[i] < j->id && finished_jobs.find(deps[i]) == finished_jobs.end())
                            dependents[deps[i]].push_back(j->id);
                        else num_deps--;
                    }

                    if(num_deps) waiting_jobs[j->id] = pair<int,job*>(num_deps, j);
                    else add_job(j);
                    pthread_mutex_unlock(&mutex);

                    write(new_sock, &j->id, sizeof(j->id));
                }
            }
        }
        close(new_sock);
    }
}

int main(int argc,char* argv[])
{
    po::options_description desc("Options");

    string name="/tmp/batch_queue_default",pid_file="";
    int num_processes = get_number_cpus();
    vector<int> deps;
    vector<string> tokens;
    desc.add_options()
        ("help,h","usage info")
        ("processes,p",po::value<int>(&num_processes)->default_value(num_processes),"number of processes to run")
        ("name,n",po::value<string>(&name)->default_value(name),"batch queue name");

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

    struct sockaddr_un sock_name;
    if(name.length()>=sizeof(sock_name.sun_path))
    {
        cout<<"name too long"<<endl;
        return 1;
    }

    sock = socket(PF_LOCAL, SOCK_STREAM, 0);
    if(sock < 0)
    {
        cout<<"failed to create socket"<<endl;
        return 1;
    }

    sock_name.sun_family = AF_LOCAL;
    strcpy(sock_name.sun_path, name.c_str());
    unlink(name.c_str());

    if(bind(sock, (sockaddr *)&sock_name, SUN_LEN(&sock_name)) < 0)
    {
        cout<<"failed to bind socket   "<<errno<<endl;
        return 1;
    }

    if(listen(sock, 1000) < 0)
    {
        cout<<"failed to listen"<<endl;
        return 1;
    }

    pthread_t listen_thread;
    threads.resize(num_processes);
    for(int i = 0; i < num_processes; i++)
        if(pthread_create(&threads[i], NULL, run_jobs, 0))
            {
                cout<<"failed to create thread"<<endl;
                return 1;
            }

    if(pthread_create(&listen_thread, NULL, listen_for_jobs, 0))
    {
        cout<<"failed to create thread"<<endl;
        return 1;
    }

    void* ret = 0;
    for(size_t i = 0; i < threads.size(); i++)
        pthread_join(threads[i], &ret);
    pthread_cancel(listen_thread);

    unlink(name.c_str());

    pthread_mutex_destroy(&mutex);

    return 0;
}
