#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include "fixed_vector.h"
#include "request_logging.h"

request_log_info * log_info = 0;

void enable_logging(const char* file)
{
    if(!log_info)
        log_info = new request_log_info(file);
}

void finish_logging()
{
    delete log_info;
    log_info = 0;
}

request_log_info::request_log_info(const char* file)
    :request_log_file(file),out(file),next_simulation_id(0)
{
}

request_log_info::~request_log_info()
{
}

request_log_simulation_info::request_log_simulation_info()
    :id(-1),next_object(1)
{
}
request_log_simulation_info::~request_log_simulation_info()
{
}

