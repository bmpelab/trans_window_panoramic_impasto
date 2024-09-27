import cv2
import os
from natsort import natsorted
import argparse
import open3d
import numpy as np
import time
import glob
import matplotlib.pyplot as plt
import keyboard
from tqdm import tqdm

def load_ply_file(filename):
    # Load the PLY point cloud
    point_cloud = open3d.io.read_point_cloud(filename)
    return point_cloud

def load_stl_file(filename):
    # Load the PLY point cloud
    stl = open3d.io.read_triangle_mesh(filename)
    return stl

def vis_dict(ply1_dir,ply2_dir,ply3_dir,ply4_dir):
    ply1_names = natsorted(glob.glob('{}/*.ply'.format(ply1_dir)))
    ply2_names = natsorted(glob.glob('{}/*.ply'.format(ply2_dir)))
    ply3_names = natsorted(glob.glob('{}/*.ply'.format(ply3_dir)))
    # ply4_names = natsorted(glob.glob('{}/*.ply'.format(ply4_dir)))
    vis = open3d.visualization.VisualizerWithKeyCallback()
    vis.create_window()
    vis.get_render_option().background_color = [0.3, 0.3, 0.3]  # [R, G, B] values range from 0 to 1
    vis.get_render_option().mesh_show_back_face = True
    vis.get_render_option().point_size = 5
    vis.get_render_option().point_color_option = open3d.visualization.PointColorOption.Color
    idx = 0
    mode_flag = 0
    pause_flag = False
    ply1 = open3d.io.read_point_cloud(ply1_names[idx])
    vis.add_geometry(ply1)

    def mode0(event):
        nonlocal mode_flag
        mode_flag = 0

    def mode1(event):
        nonlocal mode_flag
        mode_flag = 1

    def mode2(event):
        nonlocal mode_flag
        mode_flag = 2

    def mode3(event):
        nonlocal mode_flag
        mode_flag = 3

    # def mode4(event):
    #     nonlocal mode_flag
    #     mode_flag = 4

    def mode_pause(event):
        nonlocal pause_flag
        if pause_flag:
            pause_flag = False
        else:
            pause_flag = True

    keyboard.on_press_key('0', mode0)
    keyboard.on_press_key('1', mode1)
    keyboard.on_press_key('2', mode2)
    keyboard.on_press_key('3', mode3)
    # keyboard.on_press_key('4', mode4)
    keyboard.on_press_key('P', mode_pause)

    def auto_update_next(vis):
        nonlocal idx
        nonlocal mode_flag
        nonlocal pause_flag
        ply1s = []
        ply2s = []
        ply3s = []
        ply4s = []
        for i in tqdm(range(len(ply1_names)), desc="Loading", unit="iteration"):
            ply1s.append(open3d.io.read_point_cloud(ply1_names[i]))
            ply2s.append(open3d.io.read_point_cloud(ply2_names[i]))
            ply3s.append(open3d.io.read_point_cloud(ply3_names[i]))
            # ply4s.append(open3d.io.read_point_cloud(ply4_names[i]))

        idx = 0
        mode_flag = 0
        pause_flag = False
        while True:
            if not pause_flag:
                if idx+1 < len(ply1_names):
                    idx = idx + 1
                else:
                    idx = 0

            vis.clear_geometries()

            vis.get_render_option().point_color_option = open3d.visualization.PointColorOption.Color
            if mode_flag == 0:
                vis.get_render_option().point_size = 5
                vis.add_geometry(ply1s[idx], False)
                vis.add_geometry(ply3s[idx], False)
                time.sleep(0.02)
            elif mode_flag == 1:
                vis.get_render_option().point_size = 5
                vis.add_geometry(ply2s[idx], False)
                vis.add_geometry(ply3s[idx], False)
                time.sleep(0.02)
            elif mode_flag == 2:
                vis.get_render_option().point_size = 5
                vis.add_geometry(ply1s[idx], False)
                time.sleep(0.02)
            elif mode_flag == 3:
                vis.get_render_option().point_size = 5
                vis.add_geometry(ply2s[idx], False)
                time.sleep(0.02)
            # elif mode_flag == 4:
            #     vis.get_render_option().point_size = 4
            #     vis.add_geometry(ply4s[idx], False)
            #     time.sleep(0.02)

            print("mode: "+str(mode_flag)+", frame id: "+str(idx))

            vis.poll_events()

    def exit_key(vis):
        vis.destroy_window()

    vis.register_key_callback(ord('A'), auto_update_next)
    vis.register_key_callback(ord('E'), exit_key)
    vis.poll_events()
    vis.run()
    # vis.destroy_window()


if __name__ == '__main__':
    # vis_dict("G:/chen/Documents/DataManagement/Data/202311211734_EndoNeRF/cutting_tissues_twice/mesh_texture_4",
    #          "G:/chen/Documents/DataManagement/Data/202311211734_EndoNeRF/cutting_tissues_twice/raftstereo_ply",
    #          "G:/chen/Documents/DataManagement/Data/202311211734_EndoNeRF/cutting_tissues_twice/image_ply",[])
    # vis_dict("G:/chen/Documents/DataManagement/Data/202402081824_clip_hamlyn_porcine_slam_low_defor/mesh_texture_149-370_4",
    #          "G:/chen/Documents/DataManagement/Data/202402081824_clip_hamlyn_porcine_slam_low_defor/raftstereo_ply_149-370",
    #          "G:/chen/Documents/DataManagement/Data/202402081824_clip_hamlyn_porcine_slam_low_defor/image_ply_149-370",[])
    # vis_dict("G:/chen/Documents/DataManagement/Data/202402081824_clip_hamlyn_porcine_slam_low_defor/mesh_texture_730-900_inf",
    #     "G:/chen/Documents/DataManagement/Data/202402081824_clip_hamlyn_porcine_slam_low_defor/image_ply_730-900",
    #     "G:/chen/Documents/DataManagement/Data/202402081824_clip_hamlyn_porcine_slam_low_defor/image_ply_730-900",[])
    vis_dict("G:/chen/Documents/DataManagement/Data/202311211734_EndoNeRF/pulling_soft_tissues/mesh_texture_4",
             "G:/chen/Documents/DataManagement/Data/202311211734_EndoNeRF/pulling_soft_tissues/raftstereo_ply",
             "G:/chen/Documents/DataManagement/Data/202311211734_EndoNeRF/pulling_soft_tissues/image_ply",[])
    # vis_dict("G:/chen/Documents/DataManagement/Data/202402061620_hypermap_g1g2/g2/strain_ply",
    #          "G:/chen/Documents/DataManagement/Data/202312141727_g2_scan_and_video/raftstereo_ply",
    #          "G:/chen/Documents/DataManagement/Data/202312141727_g2_scan_and_video/image_ply",
    #          "G:/chen/Documents/DataManagement/Data/202402061620_hypermap_g1g2/g2/strain_ply")