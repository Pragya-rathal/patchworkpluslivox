#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Core>

// Patchwork++-ROS
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include <rclcpp/qos.hpp>
#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <std_msgs/msg/string.hpp>
#include <type_traits>

#include "GroundSegmentationServer.hpp"
#include "Utils.hpp"

namespace patchworkpp_ros {

using utils::EigenToPointCloud2;
using utils::GetTimestamps;
using utils::PointCloud2ToEigen;

// ---------------------------------------------------------------------------
// Parameter loaders
// ---------------------------------------------------------------------------

patchwork::Params GroundSegmentationServer::loadPlusplusParamsFromROS() {
  patchwork::Params params;

  params.sensor_height = declare_parameter<double>("sensor_height", params.sensor_height);
  params.num_iter      = declare_parameter<int>("num_iter", params.num_iter);
  params.num_lpr       = declare_parameter<int>("num_lpr", params.num_lpr);
  params.num_min_pts   = declare_parameter<int>("num_min_pts", params.num_min_pts);
  params.th_seeds      = declare_parameter<double>("th_seeds", params.th_seeds);

  params.th_dist    = declare_parameter<double>("th_dist", params.th_dist);
  params.th_seeds_v = declare_parameter<double>("th_seeds_v", params.th_seeds_v);
  params.th_dist_v  = declare_parameter<double>("th_dist_v", params.th_dist_v);

  params.max_range       = declare_parameter<double>("max_range", params.max_range);
  params.min_range       = declare_parameter<double>("min_range", params.min_range);
  params.uprightness_thr = declare_parameter<double>("uprightness_thr", params.uprightness_thr);

  params.verbose = get_parameter<bool>("verbose", params.verbose);

  params.max_traversable_slope_deg = declare_parameter<double>("max_traversable_slope_deg", params.max_traversable_slope_deg);
  params.roughness_threshold       = declare_parameter<double>("roughness_threshold", params.roughness_threshold);
  params.vegetation_threshold      = declare_parameter<double>("vegetation_threshold", params.vegetation_threshold);
  params.min_patch_points          = declare_parameter<int>("min_patch_points", params.min_patch_points);
  params.terrain_patch_size        = declare_parameter<double>("terrain_patch_size", params.terrain_patch_size);
  params.publish_terrain_normals   = declare_parameter<bool>("publish_terrain_normals", params.publish_terrain_normals);

  // ToDo. Support intensity
  params.enable_RNR = false;

  return params;
}

patchwork::PatchworkParams GroundSegmentationServer::loadClassicParamsFromROS() {
  patchwork::PatchworkParams params;

  params.sensor_height = declare_parameter<double>("sensor_height", params.sensor_height);
  params.max_range     = declare_parameter<double>("max_range", params.max_range);
  params.min_range     = declare_parameter<double>("min_range", params.min_range);

  params.num_iter    = declare_parameter<int>("num_iter", params.num_iter);
  params.num_lpr     = declare_parameter<int>("num_lpr", params.num_lpr);
  params.num_min_pts = declare_parameter<int>("num_min_pts", params.num_min_pts);
  params.th_seeds    = declare_parameter<double>("th_seeds", params.th_seeds);
  params.th_dist     = declare_parameter<double>("th_dist", params.th_dist);

  params.uprightness_thr = declare_parameter<double>("uprightness_thr", params.uprightness_thr);

  params.adaptive_seed_selection_margin = declare_parameter<double>(
      "adaptive_seed_selection_margin", params.adaptive_seed_selection_margin);

  params.using_global_thr = declare_parameter<bool>("using_global_thr", params.using_global_thr);
  params.global_elevation_thr =
      declare_parameter<double>("global_elevation_thr", params.global_elevation_thr);

  params.ATAT_ON        = declare_parameter<bool>("ATAT_ON", params.ATAT_ON);
  params.max_h_for_ATAT = declare_parameter<double>("max_h_for_ATAT", params.max_h_for_ATAT);
  params.num_sectors_for_ATAT =
      declare_parameter<int>("num_sectors_for_ATAT", params.num_sectors_for_ATAT);
  params.noise_bound = declare_parameter<double>("noise_bound", params.noise_bound);

  params.verbose = declare_parameter<bool>("verbose", params.verbose);

  return params;
}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

GroundSegmentationServer::GroundSegmentationServer(const rclcpp::NodeOptions &options)
    : rclcpp::Node("patchworkpp_node", options) {
  base_frame_ = declare_parameter<std::string>("base_frame", base_frame_);
  target_frame_ = declare_parameter<std::string>("target_frame", "base_link");
  std::string cloud_topic_rear = declare_parameter<std::string>("cloud_topic_rear", "");

  tf_buffer_ = std::make_unique<tf2_ros::Buffer>(this->get_clock());
  tf_listener_ = std::make_shared<tf2_ros::TransformListener>(*tf_buffer_);

  const std::string algorithm = declare_parameter<std::string>("algorithm", "patchworkpp");

  if (algorithm == "patchwork") {
    patchwork::PatchworkParams classic_params = loadClassicParamsFromROS();
    impl_ = std::make_unique<patchwork::PatchWork>(classic_params);
    RCLCPP_INFO(get_logger(), "Algorithm: patchwork (classic)");
  } else {
    patchwork::Params plusplus_params = loadPlusplusParamsFromROS();
    impl_                             = std::make_unique<patchwork::PatchWorkpp>(plusplus_params);
    RCLCPP_INFO(get_logger(), "Algorithm: patchworkpp (default)");
  }

  // Initialize subscribers
  pointcloud_sub_ = create_subscription<sensor_msgs::msg::PointCloud2>(
      "pointcloud_topic",
      rclcpp::SensorDataQoS(),
      std::bind(&GroundSegmentationServer::EstimateGround, this, std::placeholders::_1));

  if (!cloud_topic_rear.empty()) {
    pointcloud_rear_sub_ = create_subscription<sensor_msgs::msg::PointCloud2>(
        cloud_topic_rear,
        rclcpp::SensorDataQoS(),
        std::bind(&GroundSegmentationServer::EstimateGround, this, std::placeholders::_1));
  }

  /*
   * We use the following QoS setting for reliable ground segmentation.
   * If you want to run Patchwork++ in real-time and real-world operation,
   * please change the QoS setting
   */
  //  rclcpp::QoS qos((rclcpp::SystemDefaultsQoS().keep_last(1).durability_volatile()));
  rclcpp::QoS qos(rclcpp::QoSInitialization::from_rmw(rmw_qos_profile_default));
  qos.reliability(RMW_QOS_POLICY_RELIABILITY_RELIABLE);
  qos.durability(RMW_QOS_POLICY_DURABILITY_TRANSIENT_LOCAL);

  cloud_publisher_  = create_publisher<sensor_msgs::msg::PointCloud2>("/patchworkpp/cloud", qos);
  ground_publisher_ = create_publisher<sensor_msgs::msg::PointCloud2>("/patchworkpp/ground", qos);
  nonground_publisher_ = create_publisher<sensor_msgs::msg::PointCloud2>("/patchworkpp/nonground", qos);
  
  traversable_publisher_ = create_publisher<sensor_msgs::msg::PointCloud2>("/patchworkpp/traversable_ground", qos);
  nontraversable_publisher_ = create_publisher<sensor_msgs::msg::PointCloud2>("/patchworkpp/nontraversable_ground", qos);
  normals_publisher_ = create_publisher<sensor_msgs::msg::PointCloud2>("/patchworkpp/terrain_normals", qos);

  RCLCPP_INFO(this->get_logger(), "Patchwork++ ROS 2 node initialized");
}

// ---------------------------------------------------------------------------
// Callbacks
// ---------------------------------------------------------------------------

void GroundSegmentationServer::EstimateGround(
    const sensor_msgs::msg::PointCloud2::ConstSharedPtr &msg) {
  
  sensor_msgs::msg::PointCloud2 transformed_msg;
  if (msg->header.frame_id != target_frame_) {
    try {
      auto transform = tf_buffer_->lookupTransform(target_frame_, msg->header.frame_id, tf2::TimePointZero, rclcpp::Duration::from_seconds(0.1));
      tf2::doTransform(*msg, transformed_msg, transform);
    } catch (const tf2::TransformException & ex) {
      RCLCPP_WARN(this->get_logger(), "Could not transform %s to %s: %s",
          msg->header.frame_id.c_str(), target_frame_.c_str(), ex.what());
      return;
    }
  } else {
    transformed_msg = *msg;
  }

  const auto &cloud = patchworkpp_ros::utils::PointCloud2ToEigenMat(std::make_shared<sensor_msgs::msg::PointCloud2>(transformed_msg));

  // Estimate ground
  std::visit([&](auto &impl) { impl->estimateGround(cloud); }, impl_);
  cloud_publisher_->publish(patchworkpp_ros::utils::EigenMatToPointCloud2(cloud, transformed_msg.header));

  // Traversability
  std::visit([&](auto &impl) { 
      if constexpr (std::is_same_v<std::decay_t<decltype(impl)>, std::unique_ptr<patchwork::PatchWorkpp>>) {
          impl->analyze_traversability(); 
      }
  }, impl_);

  // Get ground and nonground
  Eigen::MatrixX3f ground    = std::visit([](auto &impl) { return impl->getGround(); }, impl_);
  Eigen::MatrixX3f nonground = std::visit([](auto &impl) { return impl->getNonground(); }, impl_);
  double time_taken          = std::visit([](auto &impl) { return impl->getTimeTaken(); }, impl_);
  (void)time_taken;  // available for debug logging if needed
  PublishClouds(ground, nonground, transformed_msg.header);

  // Publish traversability if patchworkpp
  std::visit([&](auto &impl) {
      if constexpr (std::is_same_v<std::decay_t<decltype(impl)>, std::unique_ptr<patchwork::PatchWorkpp>>) {
          PublishTraversabilityClouds(impl->getTraversableGround(), impl->getNontraversableGround(), impl->getTerrainNormals(), transformed_msg.header);
      }
  }, impl_);
}

void GroundSegmentationServer::PublishClouds(const Eigen::MatrixX3f &est_ground,
                                             const Eigen::MatrixX3f &est_nonground,
                                             const std_msgs::msg::Header header_msg) {
  std_msgs::msg::Header header = header_msg;
  header.frame_id              = target_frame_;
  ground_publisher_->publish(
      std::move(patchworkpp_ros::utils::EigenMatToPointCloud2(est_ground, header)));
  nonground_publisher_->publish(
      std::move(patchworkpp_ros::utils::EigenMatToPointCloud2(est_nonground, header)));
}

void GroundSegmentationServer::PublishTraversabilityClouds(
    const Eigen::MatrixX3f &traversable_ground,
    const Eigen::MatrixX3f &nontraversable_ground,
    const Eigen::MatrixX3f &terrain_normals,
    const std_msgs::msg::Header header_msg) {
  std_msgs::msg::Header header = header_msg;
  header.frame_id              = target_frame_;
  traversable_publisher_->publish(
      std::move(patchworkpp_ros::utils::EigenMatToPointCloud2(traversable_ground, header)));
  nontraversable_publisher_->publish(
      std::move(patchworkpp_ros::utils::EigenMatToPointCloud2(nontraversable_ground, header)));
  if (terrain_normals.rows() > 0) {
      normals_publisher_->publish(
          std::move(patchworkpp_ros::utils::EigenMatToPointCloud2(terrain_normals, header)));
  }
}
}  // namespace patchworkpp_ros

#include "rclcpp_components/register_node_macro.hpp"
RCLCPP_COMPONENTS_REGISTER_NODE(patchworkpp_ros::GroundSegmentationServer)
