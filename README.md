# Trans-Window Panoramic Impasto for Online Tissue Deformation Recovery

Implementation for MICCAI 2024 paper **[Trans-Window Panoramic Impasto for Online Tissue Deformation Recovery](https://)** by [Jiahe Chen](http://), Etsuko Kobayashi, Ichiro Sakuma, and Naoki Tomii.

A novel problem setting and model-based framework for online tissue deformation recovery.

## Demo
https://github.com/bmpelab/trans_window_panoramic_impasto/assets/44491917/c95251c0-c616-4243-a203-20599780d6df


## Problem formulation
![overview](https://github.com/bmpelab/trans_window_panoramic_impasto/assets/44491917/5662d76e-e410-44e6-a462-d027681dfae1)

## Dataset preparation

You can download the well-prepared ex vivo dataset from [SurgEM](https://github.com/bmpelab/SurgEM.git) and in vivo dataset (derived from Hamlyn center dataset) from [here](http://). For the EndoNeRF dataset, please download via their [website](https://github.com/med-air/EndoNeRF) and process following the instructions below.

The structure of the dataset looks like:

```sh
├── datasets
    ├── surgem_ex_vivo (the surgem ex vivo dataset)
        ├── constraint_map
        ├── point_3d_map
        ├── mask
        ├── scene_flow
        ├── rectified_left
        ├── rectified_right
        ├── evaluation
        ├── rectifiedCamera.mat (camera parameter)
    ├── hamlyn_in_vivo
        ├── point_3d_map
        ├── scene_flow
        ├── rectified_left
        ├── rectified_right
        ├── evaluation
        ├── rectifiedCamera.mat
```

### Instruction to create your own dataset 

In case you capture images by your own system, image rectification may be necessary using the calibrated camera parameters. For image rectification, please follows [OpenCV_stereoRectify](https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html#ga617b1685d4059c6040827800e72ad2b6). For 3D reconstruction and 2D tracking, please check [RAFT-Stereo](https://github.com/princeton-vl/RAFT-Stereo.git) and [LiteFlowNet3](https://github.com/twhui/LiteFlowNet3.git). After that, you can create the scene flow combining the 3D reconstructed and 2D tracked results. If there are occlusion caused by surgical instrument, instrument segmentation is necessary.

## Code preparation



<!--We recommend using Miniconda to set up an environment:-->

```bash
cd xxx
```

<!--We managed to test our code on Ubuntu 18.04 with Python 3.6 and CUDA 10.2.-->

## Run

Although the data of all frames are prepared, the code iteratively load the data frame by frame; thus, it is a online approach that making use of only the current and previous frames information. We prepare two scripts for your convenience to run the test on surgem_ex_vivo and hamlyn_in_vivo dataset.

First clone or download this repository. Then, open `main_surgem.m` or `main_hamlyn.m` in MATLAB, adjust the `data_folder` to the one of the dataset.

After that, run the code. New folders (`mesh`) containing the results will be created under the same folder of the dataset.

## Evaluation

We provide ex vivo dataset with ground truth of the deformed surface obtained by 3D scanning. Please download the dataset from [here](https://github.com/bmpelab/SurgEM.git). We compare the recovered surface with the reference surface (the scanned one) for evaluating the recovery accuracy in terms of surface distance (defined in [J. Chen, et al., IJCARS, 2023](https://doi.org/10.1007/s11548-023-02889-z)). For details about the evaluation, please check the dataset.

```bash
python run_endonerf.py --config configs/{your_config_file}.txt --render_only
```

The rendered images will be saved to `logs/{expname}/renderonly_path_fixidentity_{epoch}/estim/`. Then, you can type the command below to acquire quantitative results:

```bash
python eval_rgb.py --gt_dir /path/to/data/images --mask_dir /path/to/data/gt_masks --img_dir logs/{expname}/renderonly_path_fixidentity_{epoch}/estim/
```

Note that we only evaluate photometric errors due to the difficulties in collecting geometric ground truth. 

## Bibtex

If you find this work helpful, you can cite our paper as follows:

```
@InProceedings{chen_trans-window_2024,
 author = {Chen, Jiahe and Kobayashi, Etsuko and Sakuma, Ichiro and Tomii, Naoki},
 title = {{Trans-Window Panoramic Impasto for Online Tissue Deformation Recovery}},
 booktitle = {Medical {{Image Computing}} and {{Computer Assisted Intervention}} -- {{MICCAI}} 2024},
 pages = {xxx--xxx},
 year = {2024},
}
```

If you use our ex vivo datasets, please cite this work also:

```
@article{chen2024surgem,
  title={SurgEM: A Vision-based Surgery Environment Modeling Framework for Constructing a Digital Twin towards Autonomous Soft Tissue Manipulation},
  author={Chen, Jiahe and Kobayashi, Etsuko and Sakuma, Ichiro and Tomii, Naoki},
  journal={IEEE Robotics and Automation Letters},
  year={2024},
  publisher={IEEE}
}
```

If you use the Hamlyn center datasets, please appropriately cite their works following the instructions in [Hamlyn Centre Laparoscopic / Endoscopic Video Datasets](https://hamlyn.doc.ic.ac.uk/vision/).

## Acknowledgement

- readflow code
- natsort
