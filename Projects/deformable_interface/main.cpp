#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc,char* argv[])
{
    printf("hello\n");
    void* handle = dlopen("libInterface.so", RTLD_LAZY);
//    void* handle = dlopen("libPhysBAM.so", RTLD_LAZY);
    if(!handle){
        const char *p = dlerror();
        printf("error %s", p);
        exit(1);
    }

    printf("%p\n", handle);
    void (*init)() = (void (*)()) dlsym(handle, "init");
    void (*quit)() = (void (*)()) dlsym(handle, "quit");
    void (*run)(int argc, char* argv[]) = (void (*)(int argc, char* argv[])) dlsym(handle, "run");

    printf("%p %p %p\n", init, quit, run);

    init();

    run(argc, argv);

    quit();

    return 0;
}
//#####################################################################
