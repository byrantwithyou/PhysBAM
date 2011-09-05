import os
import random
# start with same seed each time so eye twitching is deterministic.
# IMPORTANT NOTE: because we are using this method, first_frame must ALWAYS
# be 1 so that we get the same sequence of random numbers
random.seed(0)

#your input options
first_frame = 1
end_frame = 48
output_directory = '\\\\subset\\Face_Examples\\Synthetic_Sadness_2\\'
destination_directory = 'jychong@collision:/local/jychong/render_sadness2/Input'


#program options
eyelid_frame = 0
eye_reference_frame = 1
eye_perturbation_x = 0.0
eye_perturbation_y = 0.0
eye_perturbation_z = 0.0
#next_eye_movement_frame = first_frame + random.randint(20,40)
next_eye_movement_frame = first_frame + random.randint(10,20)

#extract raw frames
os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\SOLIDS_OPENGL_CONVERT\\Release\\SOLIDS_OPENGL_CONVERT.exe '
           + output_directory)              

#extract markers
os.system('c:\\Personal_Libraries\\Sergey_Library\\marker_creator\\Release\\marker_creator.exe '
           + str(first_frame) + ' ' + str(end_frame) + ' ' + output_directory)              



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
        #next_eye_movement_frame = i + random.randint(20,40)
        next_eye_movement_frame = i + random.randint(10,20)

        
    os.system('c:\\Personal_Libraries\\Mike_T_Library\\morph_using_correspondences\\Release\\morph_using_correspondences.exe '
              + 'c:\\Public_Data\\Face_Data\\Eftychis_840k\\eftychis_neck.corr '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\output_frame' + str(i) + '.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\head_face_full_no_collar_840k.tri '
              + 'temp_morphed_surface.tri')
   
    os.system('c:\\Personal_Libraries\\Mike_T_Library\\merge_tris\\Release\\model_transformations.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Eyelids\\eyelid_1_'
              + str(eyelid_frame) + '.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Eyelids\\eyelid_2_'
              + str(eyelid_frame) + '.tri '
              + 'merged_upper_eyelids.tri')
    os.system('c:\\Personal_Libraries\\Mike_T_Library\\merge_tris\\Release\\model_transformations.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\lower_eyelid_left.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\lower_eyelid_right.tri '
              + 'merged_lower_eyelids.tri')
    os.system('c:\\Personal_Libraries\\Mike_T_Library\\merge_tris\\Release\\model_transformations.exe '
              + 'merged_upper_eyelids.tri '
              + 'merged_lower_eyelids.tri '
              + 'merged_eyelids.tri')

    os.system('c:\\Personal_Libraries\\Mike_T_Library\\merge_tris\\Release\\model_transformations.exe '
              + 'temp_morphed_surface.tri '  
              + 'merged_eyelids.tri '
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
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_render_surface.' + str(i))
    #global transform lower teeth
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\lowerteeth_frame' + str(i) + '.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_lowerteeth.' + str(i))
    #global transform upper teeth
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\upperteeth_face_fitted.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_upperteeth.' + str(i))
    #global transform lower mouth
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\lowermouth_frame' + str(i) + '.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_lowermouth.' + str(i))
    #global transform left eye
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'left_eye_temp.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_left_eye.' + str(i))
    #global transform right eye
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'right_eye_temp.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_right_eye.' + str(i))
    #global transform neck and shoulders
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Global_Tri_Transform\\Release\\Global_Tri_Transform.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\transforms_for_frame' + str(i) + ' '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\neck_shoulder.tri '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_neck_shoulders.' + str(i))


    # finally (as the last step!), flip all of the global-tranformed models

    #flip face
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_render_surface.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_render_surface.' + str(i))
    #flip lower teeth
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_lowerteeth.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_lowerteeth.' + str(i))
    #flip upper teeth
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_upperteeth.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_upperteeth.' + str(i))
    #flip lower mouth
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_lowermouth.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_lowermouth.' + str(i))
    #flip left eye
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_left_eye.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_left_eye.' + str(i))
    #flip right eye
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_right_eye.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_right_eye.' + str(i))
    #flip neck and shoulders
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform\\global_neck_shoulders.' + str(i) + ' ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_neck_shoulders.' + str(i))
    #flip markers
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\markers.' + str(i) + '.tri ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_markers.' + str(i))
    #flip embedded markers
    os.system('c:\\Personal_Libraries\\"Jiayi Folder"\\Flip_Tri\\Release\\Flip_Tri.exe '
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\test\\embedded_markers.' + str(i) + '.tri ' 
              + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_embedded_markers.' + str(i))

    #upload face
    os.system('c:\\pscp.exe ' + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_render_surface.' + str(i)
               + ' ' + destination_directory)
    #upload lower teeth
    os.system('c:\\pscp.exe ' + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_lowerteeth.' + str(i)
               + ' ' + destination_directory)
    #upload upper teeth
    os.system('c:\\pscp.exe ' + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_upperteeth.' + str(i)
               + ' ' + destination_directory)
    #upload lower mouth
    os.system('c:\\pscp.exe ' + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_lowermouth.' + str(i)
               + ' ' + destination_directory)
    #upload left eye
    os.system('c:\\pscp.exe ' + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_left_eye.' + str(i)
               + ' ' + destination_directory)
    #upload right eye
    os.system('c:\\pscp.exe ' + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_right_eye.' + str(i)
               + ' ' + destination_directory)
    #upload neck and shoulders
    os.system('c:\\pscp.exe ' + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_neck_shoulders.' + str(i)
               + ' ' + destination_directory)
    #upload markers
    os.system('c:\\pscp.exe ' + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_markers.' + str(i)
               + ' ' + destination_directory)
    #upload embedded markers
    os.system('c:\\pscp.exe ' + 'c:\\Public_Data\\Face_Data\\Face_Meshes\\Animations\\Render_Full_Face_Tri_Global_Transform_Flipped\\global_embedded_markers.' + str(i)
               + ' ' + destination_directory)

