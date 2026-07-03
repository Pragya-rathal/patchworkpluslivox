from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, ExecuteProcess
from launch.conditions import IfCondition
from launch.substitutions import (LaunchConfiguration, PathJoinSubstitution,
                                  PythonExpression)
from launch_ros.actions import Node
from launch_ros.substitutions import FindPackageShare


# This configuration parameters are not exposed thorugh the launch system, meaning you can't modify
# those throw the ros launch CLI. If you need to change these values, you could write your own
# launch file and modify the 'parameters=' block from the Node class.
class config:
    # TBU. Examples are as follows:
    max_range: float = 80.0
    # deskew: bool = False


def generate_launch_description():
    algorithm = LaunchConfiguration("algorithm", default="patchworkpp")

    use_sim_time = LaunchConfiguration("use_sim_time", default="true")

    # tf tree configuration, these are the likely 3 parameters to change and nothing else
    base_frame = LaunchConfiguration("base_frame", default="base_link")

    # ROS configuration
    pointcloud_topic = LaunchConfiguration("cloud_topic")
    visualize = LaunchConfiguration("visualize", default="true")

    # Optional ros bag play
    bagfile = LaunchConfiguration("bagfile", default="")

    # Patchwork++ node
    patchworkpp_node = Node(
        package="patchworkpp",
        executable="patchworkpp_node",
        name="patchworkpp_node",
        output="screen",
        remappings=[
            ("pointcloud_topic", pointcloud_topic),
        ],
        parameters=[
            {
                # ROS node configuration
                "base_frame": base_frame,
                "use_sim_time": use_sim_time,
                # Patchwork++ configuration
                "sensor_height": 1.88,
                "num_iter": 3,  # Number of iterations for ground plane estimation using PCA.
                "num_lpr": 20,  # Maximum number of points to be selected as lowest points representative.
                "num_min_pts": 0,  # Minimum number of points to be estimated as ground plane in each patch.
                "th_seeds": 0.3,
                # threshold for lowest point representatives using in initial seeds selection of ground points.
                "th_dist": 0.125,  # threshold for thickness of ground.
                "th_seeds_v": 0.25,
                # threshold for lowest point representatives using in initial seeds selection of vertical structural points.
                "th_dist_v": 0.9,  # threshold for thickness of vertical structure.
                "max_range": 80.0,  # max_range of ground estimation area
                "min_range": 1.0,  # min_range of ground estimation area
                "uprightness_thr": 0.101,
                # threshold of uprightness using in Ground Likelihood Estimation(GLE). Please refer paper for more information about GLE.
                "verbose": True,  # display verbose info
                 "max_traversable_slope_deg": 10.0,
    "roughness_threshold": 0.06,
    "vegetation_threshold": 0.08,
    "min_patch_points": 12,
    "terrain_patch_size": 0.5,
    "publish_terrain_normals": True,
    "enable_RNR": True,
    "enable_RVPF": True,
    "enable_TGR": True,
    "num_zones": 4,
    "num_sectors_each_zone": [16, 32, 54, 32],
    "num_rings_each_zone": [2, 4, 4, 4],
    "elevation_thresholds": [0.0, 0.0, 0.0, 0.0],
    "flatness_thresholds": [0.0, 0.0, 0.0, 0.0],
    "RNR_ver_angle_thr": -15.0,
    "RNR_intensity_thr": 0.2,
    "num_iter": 3,
    "num_lpr": 20,
    "num_min_pts": 10,
    "th_seeds": 0.3,
    "th_dist": 0.125,
    "th_seeds_v": 0.25,
    "th_dist_v": 0.1,
    "max_range": 80.0,
    "min_range": 2.7,
    "uprightness_thr": 0.707,
    "adaptive_seed_selection_margin": -1.2,
    "sensor_height": 1.723,
    "target_frame": "base_link",
    "mode": "czm",
    "verbose": False,
    "visualize": True,
            },
            {"algorithm": algorithm},
        ],
    )
    rviz_node = Node(
        package="rviz2",
        executable="rviz2",
        output="screen",
        arguments=[
            "-d",
            PathJoinSubstitution(
                [FindPackageShare("patchworkpp"), "rviz", "patchworkpp.rviz"]
            ),
        ],
        condition=IfCondition(visualize),
    )

    bagfile_play = ExecuteProcess(
        cmd=["ros2", "bag", "play", bagfile],
        output="screen",
        condition=IfCondition(PythonExpression(["'", bagfile, "' != ''"])),
    )
    return LaunchDescription(
        [
            DeclareLaunchArgument(
                "algorithm",
                default_value="patchworkpp",
                description="Ground segmentation algorithm: 'patchwork' or 'patchworkpp'",
            ),
            patchworkpp_node,
            rviz_node,
            bagfile_play,
        ]
    )
