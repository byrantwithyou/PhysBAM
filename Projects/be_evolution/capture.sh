#!/bin/bash

$PHYSBAM/Projects/opengl_3d/opengl_3d$SUFFIX $1 -offscreen -so frames-%04d.png -keys $3
ffmpeg -i frames-%04d.png $2
rm frames-[0-9][0-9][0-9][0-9].png
