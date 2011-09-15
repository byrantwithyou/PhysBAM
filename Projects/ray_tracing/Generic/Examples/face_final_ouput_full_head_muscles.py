import os
import random

first_frame = 1
end_frame = 2

eyelid_frame = 0
eye_reference_frame = 1
eye_perturbation_x = 0.0
eye_perturbation_y = 0.0
eye_perturbation_z = 0.0
next_eye_movement_frame = first_frame + random.randint(20,40)
next_eye_movement_frame = 30

for i in range(first_frame, end_frame):

    if i == next_eye_movement_frame:
        eye_perturbation_x = 2.0*random.random() - 1.0
        if eye_perturbation_x < 0:
            eye_perturbation_x = -0.25 + 0.75*eye_perturbation_x
        else:
            eye_perturbation_x = 0.25 + 0.75*eye_perturbation_x
        eye_perturbation_y = 2.0*random.random() - 1.0
        if eye_perturbation_y < 0:
            eye_perturbation_y = -0.25 + 0.75*eye_perturbation_y
        else:
            eye_perturbation_y = 0.25 + 0.75*eye_perturbation_y
        eye_perturbation_z = 2.0*random.random() - 1.0
        if eye_perturbation_z < 0:
            eye_perturbation_z = -0.25 + 0.75*eye_perturbation_z
        else:
            eye_perturbation_z = 0.25 + 0.75*eye_perturbation_z
        next_eye_movement_frame = i + random.randint(20,40)
        
    os.system('c:\\Personal_Libraries\\Mike_T_Library\\morph_using_correspondences\\Release\\morph_using_correspondences.exe '
              + 'c:\\Public_Data\\Face_Data\\Eftychis_840k\\eftychis_neck.corr '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\output_frame' + str(i) + '.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\head_face_full_no_collar_840k.tri '
              + 'final_temp_surface.tri')
   
    # create left eye
    os.system('c:\\Personal_Libraries\\Mike_T_Library\\create_eye\\Release\\create_eye.exe left '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(eye_reference_frame) + ' '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' ' 
              + str(eye_perturbation_x) + ' ' + str(eye_perturbation_y) + ' ' + str(eye_perturbation_z) + ' '
              + 'left_eye_temp.tri')
    # create right eye
    os.system('c:\\Personal_Libraries\\Mike_T_Library\\create_eye\\Release\\create_eye.exe right '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(eye_reference_frame) + ' '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' ' 
              + str(eye_perturbation_x) + ' ' + str(eye_perturbation_y) + ' ' + str(eye_perturbation_z) + ' '
              + 'right_eye_temp.tri')


    # do global transformation for all objects

    #global transform face
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'final_temp_surface.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_render_surface.' + str(i))
    #global transform left eye
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'left_eye_temp.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_left_eye.' + str(i))
    #global transform right eye
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'right_eye_temp.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_right_eye.' + str(i))
    #global transform neck and shoulders
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\neck_shoulder.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_neck_shoulders.' + str(i))
    #global transform lower jaw
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test2\\lowerjaw_frame'+ str(i)+ '.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_lower_jaw.' + str(i))
    #global transform cranium
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\cranium_surface.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_cranium.' + str(i))

    # finally (as the last step!), flip all of the global-tranformed models

    #flip face
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_render_surface.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform_Flipped\\global_render_surface.' + str(i))
    #flip muscles
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\muscle_surfaces.'+ str(i) +'.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform_Flipped\\global_muscles_surface.' + str(i))
    #flip left eye
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_left_eye.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform_Flipped\\global_left_eye.' + str(i))
    #flip right eye
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_right_eye.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform_Flipped\\global_right_eye.' + str(i))
    #flip neck and shoulders
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_neck_shoulders.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform_Flipped\\global_neck_shoulders.' + str(i))
    #flip lower jaw
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_lower_jaw.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform_Flipped\\global_lower_jaw.' + str(i))
    #flip cranium
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform\\global_cranium.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Muscles_Tri_Global_Transform_Flipped\\global_cranium.' + str(i))
