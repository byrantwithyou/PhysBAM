#include "DOUBLE_SLIDER.h"
#include "USER_INTERFACE.h"

int main(int argc, char **argv) 
{
#if 1
    Fl::scheme("plastic");
    USER_INTERFACE *user_interface = new USER_INTERFACE(argc,argv);
    user_interface->window->show();
#else
    Fl_Window win(0,0,800,600,"Test");
    DOUBLE_SLIDER<float> foo(40,10,300,60,"test");
    win.show();
#endif

    return Fl::run();
}
