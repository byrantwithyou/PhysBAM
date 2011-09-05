#merges the face animations with the blinking eyelid animations and probably
#the neck and ears in the near future

import os

eyelid_frame = 0
blink_direction = 0

for i in range(0,1):
    
    os.system('c:\\Personal_Libraries\\Mike_T_Library\\merge_tris\\Release\\model_transformations.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Eyelids\\eyelid_1_'
              + str(eyelid_frame) + '.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Eyelids\\eyelid_2_'
              + str(eyelid_frame) + '.tri '
              + 'temp_surface.tri')

    os.system('c:\\Personal_Libraries\\Mike_T_Library\\merge_tris\\Release\\model_transformations.exe '
              + 'temp_surface.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\output_frame'
              + str(i) + '.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Tri\\render_surface'
              + str(i) + '.tri')

    if blink_direction == 0:
        eyelid_frame = eyelid_frame + 2
    elif blink_direction == 1:
        eyelid_frame = eyelid_frame - 2

    if eyelid_frame >= 7:
        blink_direction = 1
    elif eyelid_frame < 0:
        blink_direction = 2

    if eyelid_frame == 6:
        eyelid_frame = 7
    elif eyelid_frame == 1:
        eyelid_frame =0
