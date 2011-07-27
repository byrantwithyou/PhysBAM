ADVECTION_NOCONA(1)

NAME
    advection_nocona - make a velocity field divergence free

SYNOPSIS
    advection_nocona [OPTION]...

DESCRIPTION
    Takes a velocity field V at position x specified by V(x)=1 and advects a passive density field D(x)=0.5*sin(4*pi*(x-0.25)-0.5*pi)+1 for x between 0.25 and 0.75 and 0 elsewhere. The output
    will be dumped into a directory named output.

    -scale RESOLUTION
        Specifies the resolution of the grid used in number of cells along each axis.

    -restart FRAME
        Specifies a frame to restart the simulation from. The intended usage is when the simulation is killed or needs to be changed after a certain frame this parameter should be used to avoid
        running the entire simulation again.

    -e FRAME
        Specifies the end frame of the simulation.

    -substep LEVEL
        Specifies the amount of output written. The higher the level the more output is written and the slower the simulation is.

    -2d
        Runs the simulation in two dimensions instead of one
    
    -3d
        Runs the simulation in three dimensions instead of one

    -threads NUMBER
        EXPERIMENTAL This parameter specifies the number of threads to use but is still in the development stages and may not give the desired performance.
