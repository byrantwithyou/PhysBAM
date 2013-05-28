#include "colors.inc"
#include "textures.inc"

#default{finish{ambient 0.2 diffuse 0.6}}

camera{
  location <3,5,3>
  look_at <0,0,0>
  up <0,1,0>
  angle 15   
}

#declare sphere_radius=.002;

union{

#emit sim_particles

  texture{
    pigment{color rgbt <1,0,0,0>}
    finish{
      ambient 0.2
      diffuse .75
      phong 0.25
      phong_size .5
      }}
}

union{

#emit sim_surface

  texture{
    pigment{color rgbt <0,.5,1,0>}
    finish{
      ambient 0.2
      diffuse .75
      phong 0.25
      phong_size .5
      }}
}

light_source{<100,100,50> color White*0.4*0.3
    area_light <10, 0, 0>, <0, 0, 10>, 10, 10
    adaptive 1
    jitter}

light_source{<2800,2500,700> color White*0.5 parallel point_at 0 shadowless}

#declare intensity=.06;

light_source{ <100,100,0>   color rgb <0.4,0.4,0.4>*0.9
              spotlight
              point_at<3,0,7>
              radius 0  // hotspot
              tightness 0
              falloff 0
              translate< 0, 0, 0>
              shadowless
}

light_source{ <100,100,25>   color rgb <0.2,0.2,0.2>*0.9
              spotlight
              point_at<-50,0,-50>
              radius 0  // hotspot
              tightness 1
              falloff 1
              translate< 0.1, 0.1, 0.1>
              shadowless
}