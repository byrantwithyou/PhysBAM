import os

f = open('reference_scene.scene', 'r')
file_string = f.read()

for i in range(27,28):
    # just the face
#    modified_file_string = file_string.replace('<face_model_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Tri/render_surface' + str(i) + '.tri')
#    modified_file_string = modified_file_string.replace('<face_uv_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Tri/render_surface0.uv')
#    modified_file_string = modified_file_string.replace('<face_samples_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Tri/render_surface5_s4.5_100k.samples')  

    # the whole head
#    modified_file_string = file_string.replace('<face_model_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri/render_surface' + str(i) + '.tri')
#    modified_file_string = modified_file_string.replace('<face_uv_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri/render_surface0_closest_point.uv')
#    modified_file_string = modified_file_string.replace('<face_samples_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri/render_surface5_s4.5_500k.samples')  
#    modified_file_string = modified_file_string.replace('<lowerteeth_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri/lowerteeth_frame_flipped' + str(i) + '.tri')  
#    modified_file_string = modified_file_string.replace('<left_eye_filename>','../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri/left_eye' + str(i) + '.tri')
#    modified_file_string = modified_file_string.replace('<right_eye_filename>','../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri/right_eye' + str(i) + '.tri')

    modified_file_string = file_string.replace('<face_model_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri_Global_Transform_Flipped/global_render_surface.' + str(i))
    modified_file_string = modified_file_string.replace('<face_uv_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri/render_surface0_closest_point.uv')
    modified_file_string = modified_file_string.replace('<face_samples_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri/render_surface0_s4.5_500k.samples')  
    modified_file_string = modified_file_string.replace('<lowerteeth_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri_Global_Transform_Flipped/global_lowerteeth.' + str(i))  
    modified_file_string = modified_file_string.replace('<upperteeth_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri_Global_Transform_Flipped/global_upperteeth.' + str(i))  
    modified_file_string = modified_file_string.replace('<lowermouth_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri_Global_Transform_Flipped/global_lowermouth.' + str(i))  
    modified_file_string = modified_file_string.replace('<lower_mouth_uv_filename>', '../../../../Public_Data/Face_Data/Face_Meshes/lower_mouth.uv')
    modified_file_string = modified_file_string.replace('<left_eye_filename>','../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri_Global_Transform_Flipped/global_left_eye.' + str(i))
    modified_file_string = modified_file_string.replace('<right_eye_filename>','../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri_Global_Transform_Flipped/global_right_eye.' + str(i))
    modified_file_string = modified_file_string.replace('<neck_and_shoulders_filename>','../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri_Global_Transform_Flipped/global_neck_shoulders.' + str(i))
    modified_file_string = modified_file_string.replace('<target_markers_model_filename>','../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri_Global_Transform_Flipped/global_markers.' + str(i))
    modified_file_string = modified_file_string.replace('<embedded_markers_model_filename>','../../../../Public_Data/Face_Data/Face_Meshes/Animations/Render_Full_Face_Tri_Global_Transform_Flipped/global_embedded_markers.' + str(i))


    f_out = open('temp_scene.scene', 'w')
    f_out.write(modified_file_string)
    f_out.close()
    
    os.system('c:\\Projects\\ray_tracing\\Release\\ray_tracing.exe temp_scene.scene ' + str(i))
