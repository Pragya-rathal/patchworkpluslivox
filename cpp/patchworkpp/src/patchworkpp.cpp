#include "patchwork/patchworkpp.h"

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cmath>

#include "patchwork/plane_fit.h"  // xy2theta, xy2radius, point_z_cmp

using namespace std;
using namespace patchwork;

Eigen::MatrixX3f PatchWorkpp::toEigenCloud(const vector<PointXYZ> &cloud) {
  Eigen::MatrixX3f dst(cloud.size(), 3);
  int j = 0;
  for (auto &p : cloud) {
    dst.row(j++) << p.x, p.y, p.z;
  }
  return dst;
}

Eigen::VectorXi PatchWorkpp::toIndices(const vector<PointXYZ> &cloud) {
  Eigen::VectorXi dst(cloud.size());
  int j = 0;
  for (auto &p : cloud) {
    dst.row(j++) << p.idx;
  }
  return dst;
}

void PatchWorkpp::addCloud(vector<PointXYZ> &cloud, const vector<PointXYZ> &add) {
  cloud.insert(cloud.end(), add.begin(), add.end());
}

void PatchWorkpp::flush_patches(vector<Zone> &czm) {
  for (int k = 0; k < params_.num_zones; k++) {
    for (int i = 0; i < params_.num_rings_each_zone[k]; i++) {
      for (int j = 0; j < params_.num_sectors_each_zone[k]; j++) {
        // czm[k][i][j].resize(MAX_POINTS, 3);
        czm[k][i][j].clear();
      }
    }
  }

  if (params_.verbose)
    cout << "\033[1;31m"
         << "PatchWorkpp::flush_patches() - Flushed patches successfully!"
         << "\033[0m" << endl;
}

void PatchWorkpp::estimate_plane(const vector<PointXYZ> &ground) {
  if (ground.empty()) return;

  // Single-pass accumulation of mean + cross-products. Avoids the
  // per-call Eigen::MatrixX3f / centered / adjoint-product heap
  // allocations that dominated the profile. Cov is computed from the
  // raw second moments via cov_ij = (sum p_i p_j - N * mean_i * mean_j) / (N - 1).
  const size_t n = ground.size();
  double sx = 0.0, sy = 0.0, sz = 0.0;
  double sxx = 0.0, syy = 0.0, szz = 0.0;
  double sxy = 0.0, sxz = 0.0, syz = 0.0;
  for (const auto &p : ground) {
    const double x = p.x, y = p.y, z = p.z;
    sx += x;
    sy += y;
    sz += z;
    sxx += x * x;
    syy += y * y;
    szz += z * z;
    sxy += x * y;
    sxz += x * z;
    syz += y * z;
  }
  const double inv_n = 1.0 / static_cast<double>(n);
  const double mx = sx * inv_n, my = sy * inv_n, mz = sz * inv_n;
  const double denom = (n > 1) ? static_cast<double>(n - 1) : 1.0;
  const double inv_d = 1.0 / denom;

  Eigen::Matrix3f cov;
  cov(0, 0) = static_cast<float>((sxx - n * mx * mx) * inv_d);
  cov(1, 1) = static_cast<float>((syy - n * my * my) * inv_d);
  cov(2, 2) = static_cast<float>((szz - n * mz * mz) * inv_d);
  cov(0, 1) = cov(1, 0) = static_cast<float>((sxy - n * mx * my) * inv_d);
  cov(0, 2) = cov(2, 0) = static_cast<float>((sxz - n * mx * mz) * inv_d);
  cov(1, 2) = cov(2, 1) = static_cast<float>((syz - n * my * mz) * inv_d);

  pc_mean_ << static_cast<float>(mx), static_cast<float>(my), static_cast<float>(mz);

  bool is_valid = true;
  if (n < 3) is_valid = false;
  if (!cov.allFinite()) is_valid = false;

  if (is_valid) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eig;
    eig.computeDirect(cov, Eigen::ComputeEigenvectors);
    const Eigen::Vector3f evals = eig.eigenvalues().cwiseMax(0.0f);
    
    if (!evals.allFinite()) is_valid = false;
    if (evals(2) < 1e-6) is_valid = false; // Degenerate/singular

    normal_ = eig.eigenvectors().col(0);
    if (!normal_.allFinite()) is_valid = false;

    if (is_valid) {
      if (normal_(2) < 0.0f) normal_ = -normal_;
      
      double mag = normal_.norm();
      if (mag < 0.99 || mag > 1.01 || std::isnan(mag)) is_valid = false;
      
      singular_values_ << evals(2), evals(1), evals(0);
      d_ = -normal_.dot(pc_mean_);
      if (!std::isfinite(d_)) is_valid = false;
    }
  }

  if (!is_valid) {
    normal_ << 0.0f, 0.0f, 1.0f;
    singular_values_ << 10000.0f, 10000.0f, 10000.0f; // High roughness to ensure rejection
    d_ = -pc_mean_(2);
  }
}

void PatchWorkpp::extract_initial_seeds(const int zone_idx,
                                        const vector<PointXYZ> &p_sorted,
                                        vector<PointXYZ> &init_seeds,
                                        double th_seed) const {
  init_seeds.clear();

  // LPR is the mean of low point representative
  double sum = 0;
  int cnt    = 0;

  int init_idx = 0;
  if (zone_idx == 0) {
    for (int i = 0; i < p_sorted.size(); i++) {
      if (p_sorted[i].z < params_.adaptive_seed_selection_margin * params_.sensor_height) {
        ++init_idx;
      } else {
        break;
      }
    }
  }

  int adaptive_num_lpr = params_.num_lpr;
  if (!p_sorted.empty()) {
    adaptive_num_lpr = std::min(params_.num_lpr, std::max(3, static_cast<int>(p_sorted.size() * 0.1)));
  }

  // Calculate the mean height value.
  for (int i = init_idx; i < p_sorted.size() && cnt < adaptive_num_lpr; i++) {
    sum += p_sorted[i].z;
    cnt++;
  }
  double lpr_height = cnt != 0 ? sum / cnt : 0;  // in case divide by 0

  int init_seeds_num =
      0;  // iterate pointcloud, filter those height is less than lpr.height+params_.th_seeds
  for (int i = 0; i < p_sorted.size(); i++) {
    if (p_sorted[i].z < lpr_height + th_seed) {
      init_seeds.push_back(p_sorted[i]);
    }
  }
}

void PatchWorkpp::extract_initial_seeds(const int zone_idx,
                                        const vector<PointXYZ> &p_sorted,
                                        vector<PointXYZ> &init_seeds) const {
  init_seeds.clear();

  // LPR is the mean of low point representative
  double sum = 0;
  int cnt    = 0;

  int init_idx = 0;
  if (zone_idx == 0) {
    for (int i = 0; i < p_sorted.size(); i++) {
      if (p_sorted[i].z < params_.adaptive_seed_selection_margin * params_.sensor_height) {
        ++init_idx;
      } else {
        break;
      }
    }
  }

  int adaptive_num_lpr = params_.num_lpr;
  if (!p_sorted.empty()) {
    adaptive_num_lpr = std::min(params_.num_lpr, std::max(3, static_cast<int>(p_sorted.size() * 0.1)));
  }

  // Calculate the mean height value.
  for (int i = init_idx; i < p_sorted.size() && cnt < adaptive_num_lpr; i++) {
    sum += p_sorted[i].z;
    cnt++;
  }
  double lpr_height = cnt != 0 ? sum / cnt : 0;  // in case divide by 0

  int init_seeds_num = 0;
  // iterate pointcloud, filter those height is less than lpr.height+params_.th_seeds
  for (int i = 0; i < p_sorted.size(); i++) {
    if (p_sorted[i].z < lpr_height + params_.th_seeds) {
      init_seeds.push_back(p_sorted[i]);
    }
  }
}

void PatchWorkpp::estimateGround(Eigen::MatrixXf cloud_in) {
  cloud_ground_.clear();
  cloud_nonground_.clear();

  if (params_.verbose)
    cout << "\033[1;32m"
         << "PatchWorkpp::estimateGround() - Estimation starts !"
         << "\033[0m" << endl;

  clock_t beg = clock();

  // 1. Reflected Noise Removal (RNR)
  if (params_.enable_RNR) reflected_noise_removal(cloud_in);

  clock_t t1 = clock();

  // 2. Concentric Zone Model (CZM)
  flush_patches(ConcentricZoneModel_);

  clock_t t1_1 = clock();

  pc2czm(cloud_in, ConcentricZoneModel_);

  clock_t t2 = clock();

  int concentric_idx = 0;

  centers_.clear();
  normals_.clear();

  double t_flush = t1_1 - t1, t_czm = t2 - t1_1, t_sort = 0.0, t_pca = 0.0, t_gle = 0.0,
         t_revert = 0.0, t_update = 0.0;

  std::vector<patchwork::RevertCandidate> candidates;
  std::vector<double> ringwise_flatness;

  // NOTE: TBB parallelisation was evaluated for this main loop and
  // measurably HURT throughput on KITTI (24-core / 8-core / 4-core all
  // 30-50% slower than single-thread). The per-patch work is small
  // (~14 µs avg) and dominated by short-lived `std::vector` /
  // `Eigen::Matrix` allocations inside R-VPF + R-GPF, so concurrent
  // mallocs serialise on the heap and TBB scheduler overhead exceeds
  // the parallelisation benefit. Single-threaded Patchwork++ already
  // runs ~110 Hz on KITTI HDL-64E (2× the paper's reported 55 Hz on
  // i7-7700K), so there is no real-time motivation to parallelise.
  // The classic Patchwork (see cpp/patchwork/src/patchwork.cpp) does
  // benefit from TBB because it has no R-VPF and fewer allocations
  // per patch.
  for (int zone_idx = 0; zone_idx < params_.num_zones; ++zone_idx) {
    auto &zone = ConcentricZoneModel_[zone_idx];

    for (int ring_idx = 0; ring_idx < params_.num_rings_each_zone[zone_idx]; ++ring_idx) {
      const int num_sectors = params_.num_sectors_each_zone[zone_idx];
      
      // Adaptive min points based on range/zone to support sparse Livox data
      int adaptive_min_pts = params_.num_min_pts;
      if (zone_idx == 1) adaptive_min_pts = params_.num_min_pts * 0.75;
      if (zone_idx == 2) adaptive_min_pts = params_.num_min_pts * 0.5;
      if (zone_idx == 3) adaptive_min_pts = params_.num_min_pts * 0.25;
      if (adaptive_min_pts < 3) adaptive_min_pts = 3; // Absolute minimum for PCA

      // Adaptive Sector Merging: merge underpopulated sectors with next sector in the same ring
      for (int sector_idx = 0; sector_idx < num_sectors - 1; ++sector_idx) {
        if (!zone[ring_idx][sector_idx].empty() && zone[ring_idx][sector_idx].size() < static_cast<size_t>(adaptive_min_pts)) {
          zone[ring_idx][sector_idx + 1].insert(zone[ring_idx][sector_idx + 1].end(), 
                                                zone[ring_idx][sector_idx].begin(), 
                                                zone[ring_idx][sector_idx].end());
          zone[ring_idx][sector_idx].clear();
        }
      }

      clock_t t_bef_gle = clock();
      for (int sector_idx = 0; sector_idx < num_sectors; ++sector_idx) {
        if (zone[ring_idx][sector_idx].size() < static_cast<size_t>(adaptive_min_pts)) {
          addCloud(cloud_nonground_, zone[ring_idx][sector_idx]);
          continue;
        }

        std::sort(
            zone[ring_idx][sector_idx].begin(), zone[ring_idx][sector_idx].end(), point_z_cmp);

        extract_piecewiseground(
            zone_idx, zone[ring_idx][sector_idx], regionwise_ground_, regionwise_nonground_);

        centers_.push_back(PointXYZ(pc_mean_(0), pc_mean_(1), pc_mean_(2)));
        normals_.push_back(PointXYZ(normal_(0), normal_(1), normal_(2)));

        const double ground_uprightness = normal_(2);
        const double ground_elevation   = pc_mean_(2);
        const double ground_flatness    = singular_values_.minCoeff();
        const double line_variable      = singular_values_(1) != 0
                                              ? singular_values_(0) / singular_values_(1)
                                              : std::numeric_limits<double>::max();

        double heading = 0.0;
        for (int i = 0; i < 3; ++i) heading += pc_mean_(i) * normal_(i);

        const bool is_upright         = ground_uprightness > params_.uprightness_thr;
        const bool is_near_zone       = concentric_idx < params_.num_rings_of_interest;
        const bool is_heading_outside = heading < 0.0;

        bool is_not_elevated = false;
        bool is_flat         = false;
        if (concentric_idx < params_.num_rings_of_interest) {
          is_not_elevated = ground_elevation < params_.elevation_thr[concentric_idx];
          is_flat         = ground_flatness < params_.flatness_thr[concentric_idx];
        }

        if (is_upright && is_not_elevated && is_near_zone) {
          update_elevation_[concentric_idx].push_back(ground_elevation);
          update_flatness_[concentric_idx].push_back(ground_flatness);
          ringwise_flatness.push_back(ground_flatness);
        }

        if (!is_upright) {
          addCloud(cloud_nonground_, regionwise_ground_);
        } else if (!is_near_zone) {
          addCloud(cloud_ground_, regionwise_ground_);
        } else if (!is_heading_outside) {
          addCloud(cloud_nonground_, regionwise_ground_);
        } else if (is_not_elevated || is_flat) {
          addCloud(cloud_ground_, regionwise_ground_);
        } else {
          candidates.emplace_back(concentric_idx,
                                  sector_idx,
                                  ground_flatness,
                                  line_variable,
                                  pc_mean_,
                                  regionwise_ground_);
        }
        addCloud(cloud_nonground_, regionwise_nonground_);
      }
      clock_t t_aft_gle = clock();
      t_gle += t_aft_gle - t_bef_gle;

      clock_t t_bef_revert = clock();
      if (!candidates.empty()) {
        if (params_.enable_TGR) {
          temporal_ground_revert(ringwise_flatness, candidates, concentric_idx);
        } else {
          for (auto &candidate : candidates) {
            addCloud(cloud_nonground_, candidate.regionwise_ground);
          }
        }

        candidates.clear();
      }
      // ringwise_flatness must be cleared every ring; the previous
      // placement inside `if (!candidates.empty())` leaked flatnesses
      // from no-candidate rings into the next ring's TGR statistics
      // (see issue #69).
      ringwise_flatness.clear();
      clock_t t_aft_revert = clock();

      t_revert += t_aft_revert - t_bef_revert;

      concentric_idx++;
    }
  }

  clock_t t_bef_update = clock();
  update_elevation_thr();
  update_flatness_thr();
  clock_t t_aft_update = clock();

  t_update = t_aft_update - t_bef_update;

  clock_t end = clock();
  time_taken_ = end - beg;

  if (params_.verbose) {
    cout << "Time taken : " << time_taken_ / static_cast<double>(1000000)
         << "(sec) ~ "
         //  << t_flush / double(1000000)  << "(flush) + "
         << t_czm / static_cast<double>(1000000) << "(czm) + "
         << t_sort / static_cast<double>(1000000) << "(sort) + "
         << t_pca / static_cast<double>(1000000) << "(pca) + "
         << t_gle / static_cast<double>(1000000) << "(estimate)" << endl;
    //  << t_revert / double(1000000) << "(revert) + "
    //  << t_update / double(1000000) << "(update)" << endl;
  }

  if (params_.verbose)
    cout << "\033[1;32m"
         << "PatchWorkpp::estimateGround() - Estimation is finished !"
         << "\033[0m" << endl;
}

void PatchWorkpp::update_elevation_thr(void) {
  for (int i = 0; i < params_.num_rings_of_interest; i++) {
    if (update_elevation_[i].empty()) continue;

    double update_mean = 0.0, update_stdev = 0.0;
    calc_mean_stdev(update_elevation_[i], update_mean, update_stdev);
    if (i == 0) {
      params_.elevation_thr[i] = update_mean + 3 * update_stdev;
      params_.sensor_height    = -update_mean;
    } else {
      params_.elevation_thr[i] = update_mean + 2 * update_stdev;
    }

    // if (params_.verbose) cout << "elevation threshold [" << i << "]: " <<
    // params_.elevation_thr[i] << endl;

    int exceed_num = update_elevation_[i].size() - params_.max_elevation_storage;
    if (exceed_num > 0)
      update_elevation_[i].erase(update_elevation_[i].begin(),
                                 update_elevation_[i].begin() + exceed_num);
  }
}

void PatchWorkpp::update_flatness_thr(void) {
  for (int i = 0; i < params_.num_rings_of_interest; i++) {
    if (update_flatness_[i].empty()) break;
    if (update_flatness_[i].size() <= 1) break;

    double update_mean = 0.0, update_stdev = 0.0;
    calc_mean_stdev(update_flatness_[i], update_mean, update_stdev);
    params_.flatness_thr[i] = update_mean + update_stdev;

    // if (params_.verbose) { cout << "flatness threshold [" << i << "]: " <<
    // params_.flatness_thr[i] << endl; }

    int exceed_num = update_flatness_[i].size() - params_.max_flatness_storage;
    if (exceed_num > 0)
      update_flatness_[i].erase(update_flatness_[i].begin(),
                                update_flatness_[i].begin() + exceed_num);
  }
}

void PatchWorkpp::reflected_noise_removal(Eigen::MatrixXf &cloud_in) {
  if (cloud_in.cols() < 4) {
    cout << "RNR requires intensity information !" << endl;
    return;
  }

  double min_intensity = std::numeric_limits<double>::max();
  double max_intensity = std::numeric_limits<double>::lowest();
  for (int i = 0; i < cloud_in.rows(); i++) {
     double I = cloud_in.row(i)(3);
     if (I < min_intensity) min_intensity = I;
     if (I > max_intensity) max_intensity = I;
  }
  double intensity_range = max_intensity - min_intensity;
  if (intensity_range < 1e-3) intensity_range = 1.0;

  int cnt = 0;
  for (int i = 0; i < cloud_in.rows(); i++) {
    double r =
        sqrt(cloud_in.row(i)(0) * cloud_in.row(i)(0) + cloud_in.row(i)(1) * cloud_in.row(i)(1));
    double z                = cloud_in.row(i)(2);
    double ver_angle_in_deg = atan2(z, r) * 180 / M_PI;
    
    double normalized_I = (cloud_in.row(i)(3) - min_intensity) / intensity_range;

    if (ver_angle_in_deg < params_.RNR_ver_angle_thr && z < -params_.sensor_height - 0.8 &&
        normalized_I < params_.RNR_intensity_thr) {
      cloud_nonground_.push_back(
          PointXYZ(cloud_in.row(i)(0), cloud_in.row(i)(1), cloud_in.row(i)(2), i));
      cloud_in.row(i)(2) = std::numeric_limits<float>::min();
      cnt++;
    }
  }

  if (params_.verbose)
    cout << "PatchWorkpp::reflected_noise_removal() - Number of Noises : " << cnt << endl;
}

void PatchWorkpp::temporal_ground_revert(const std::vector<double> &ring_flatness,
                                         const std::vector<patchwork::RevertCandidate> &candidates,
                                         int concentric_idx) {
  if (params_.verbose)
    std::cout << "\033[1;34m"
              << "=========== Temporal Ground Revert (TGR) ==========="
              << "\033[0m" << endl;

  double mean_flatness = 0.0, stdev_flatness = 0.0;
  calc_mean_stdev(ring_flatness, mean_flatness, stdev_flatness);

  if (params_.verbose) {
    cout << "[" << candidates[0].concentric_idx << ", " << candidates[0].sector_idx << "]"
         << " mean_flatness: " << mean_flatness << ", stdev_flatness: " << stdev_flatness
         << std::endl;
  }

  for (const auto &candidate : candidates) {
    // Debug
    if (params_.verbose) {
      cout << "\033[1;33m" << candidate.sector_idx << "th flat_sector_candidate"
           << " / flatness: " << candidate.ground_flatness
           << " / line_variable: " << candidate.line_variable
           << " / ground_num : " << candidate.regionwise_ground.size() << "\033[0m" << endl;
    }

    double mu_flatness = mean_flatness + 1.5 * stdev_flatness;
    if (ring_flatness.size() < 3 || mu_flatness < 1e-5) {
       // Graceful fallback to learned/static threshold for Livox sparse data
       mu_flatness = params_.flatness_thr[concentric_idx];
       if (mu_flatness < 1e-5) mu_flatness = 0.05; // safe default
    }
    double prob_flatness =
        1 / (1 + exp((candidate.ground_flatness - mu_flatness) / (mu_flatness / 10)));

    if (candidate.regionwise_ground.size() > 1500 &&
        candidate.ground_flatness < params_.th_dist * params_.th_dist)
      prob_flatness = 1.0;

    double prob_line = 1.0;
    if (candidate.line_variable > 8.0) {  //&& candidate.line_dir > M_PI/4)//
      // if (params_.verbose) cout << "line_dir: " << candidate.line_dir << endl;
      prob_line = 0.0;
    }

    bool revert = prob_line * prob_flatness > 0.5;

    if (concentric_idx < params_.num_rings_of_interest) {
      if (revert) {
        if (params_.verbose) {
          cout << "\033[1;32m"
               << "REVERT TRUE"
               << "\033[0m" << endl;
        }
        addCloud(cloud_ground_, candidate.regionwise_ground);
      } else {
        if (params_.verbose) {
          cout << "\033[1;31m"
               << "FINAL REJECT"
               << "\033[0m" << endl;
        }
        addCloud(cloud_nonground_, candidate.regionwise_ground);
      }
    }
  }

  if (params_.verbose)
    std::cout << "\033[1;34m"
              << "===================================================="
              << "\033[0m" << endl;
}

// For adaptive
void PatchWorkpp::extract_piecewiseground(const int zone_idx,
                                          const vector<PointXYZ> &src,
                                          vector<PointXYZ> &dst,
                                          vector<PointXYZ> &non_ground_dst) {
  // 0. Initialization
  if (!ground_pc_.empty()) ground_pc_.clear();
  if (!dst.empty()) dst.clear();
  if (!non_ground_dst.empty()) non_ground_dst.clear();

  // 1. Region-wise Vertical Plane Fitting (R-VPF)
  // : removes potential vertical plane under the ground plane
  // src_wo_verticals_ and src_tmp_ are reused instance scratch buffers
  // (see header) — `clear()` keeps capacity, so the per-patch malloc
  // pressure on the glibc heap goes away after the first few patches.
  src_wo_verticals_.clear();
  src_wo_verticals_.insert(src_wo_verticals_.end(), src.begin(), src.end());

  if (params_.enable_RVPF) {
    for (int i = 0; i < params_.num_iter; i++) {
      extract_initial_seeds(zone_idx, src_wo_verticals_, ground_pc_, params_.th_seeds_v);
      estimate_plane(ground_pc_);

      if (zone_idx == 0 && normal_(2) < params_.uprightness_thr) {
        src_tmp_.clear();
        src_tmp_.swap(src_wo_verticals_);  // src_tmp_ now holds the old src_wo_verticals_;
                                           // src_wo_verticals_ is empty (capacity retained).

        for (const auto &point : src_tmp_) {
          double distance = calc_point_to_plane_d(point, normal_, d_);

          if (abs(distance) < params_.th_dist_v) {
            non_ground_dst.push_back(point);
          } else {
            src_wo_verticals_.push_back(point);
          }
        }
      } else
        break;
    }
  }

  // 2. Region-wise Ground Plane Fitting (R-GPF)
  // : fits the ground plane

  extract_initial_seeds(zone_idx, src_wo_verticals_, ground_pc_);
  estimate_plane(ground_pc_);

  for (int i = 0; i < params_.num_iter; i++) {
    ground_pc_.clear();

    for (const auto &point : src_wo_verticals_) {
      double distance = calc_point_to_plane_d(point, normal_, d_);

      double adaptive_th_dist = params_.th_dist;
      if (zone_idx == 1) adaptive_th_dist *= 1.2;
      else if (zone_idx == 2) adaptive_th_dist *= 1.5;
      else if (zone_idx == 3) adaptive_th_dist *= 2.0;

      if (i < params_.num_iter - 1) {
        if (distance < adaptive_th_dist) {
          ground_pc_.push_back(point);
        }
      } else {
        if (distance < adaptive_th_dist) {
          dst.push_back(point);
        } else {
          non_ground_dst.push_back(point);
        }
      }
    }

    if (i < params_.num_iter - 1) {
      estimate_plane(ground_pc_);
    } else {
      estimate_plane(dst);
    }
  }

  if (dst.size() + non_ground_dst.size() != src.size()) {
    cout << "\033[1;33m"
         << "Points are Missing/Adding !!! Please Check !! "
         << "\033[0m" << endl;
    cout << "gnd size: " << dst.size() << ", non gnd size: " << non_ground_dst.size()
         << ", src: " << src.size() << endl;
  }
}

double PatchWorkpp::calc_point_to_plane_d(const PointXYZ &p,
                                          const Eigen::VectorXf &normal,
                                          double d) {
  return normal(0) * p.x + normal(1) * p.y + normal(2) * p.z + d;
}

void PatchWorkpp::calc_mean_stdev(const std::vector<double> &vec, double &mean, double &stdev) {
  if (vec.size() <= 1) return;

  mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();

  for (int i = 0; i < vec.size(); i++) {
    stdev += (vec.at(i) - mean) * (vec.at(i) - mean);
  }
  stdev /= vec.size() - 1;
  stdev = sqrt(stdev);
}

void PatchWorkpp::pc2czm(const Eigen::MatrixXf &src, std::vector<Zone> &czm) {
  double max_range = params_.max_range, min_range = params_.min_range;
  double min_range_0 = min_ranges_[0], min_range_1 = min_ranges_[1], min_range_2 = min_ranges_[2],
         min_range_3 = min_ranges_[3];
  int num_ring_0 = params_.num_rings_each_zone[0], num_sector_0 = params_.num_sectors_each_zone[0];
  int num_ring_1 = params_.num_rings_each_zone[1], num_sector_1 = params_.num_sectors_each_zone[1];
  int num_ring_2 = params_.num_rings_each_zone[2], num_sector_2 = params_.num_sectors_each_zone[2];
  int num_ring_3 = params_.num_rings_each_zone[3], num_sector_3 = params_.num_sectors_each_zone[3];

  for (int i = 0; i < src.rows(); i++) {
    float x = src.row(i)(0), y = src.row(i)(1), z = src.row(i)(2);

    if (z == std::numeric_limits<float>::min()) continue;

    double r = xy2radius(x, y);
    int ring_idx, sector_idx;
    if ((r <= max_range) && (r > min_range)) {
      // double theta = xy2theta(pt.x, pt.y);
      double theta = xy2theta(x, y);

      if (r < min_range_1) {  // In First rings
        ring_idx   = min(static_cast<int>(((r - min_range_0) / ring_sizes_[0])), num_ring_0 - 1);
        sector_idx = min(static_cast<int>((theta / sector_sizes_[0])), num_sector_0 - 1);
        czm[0][ring_idx][sector_idx].emplace_back(PointXYZ(x, y, z, i));
      } else if (r < min_range_2) {
        ring_idx   = min(static_cast<int>(((r - min_range_1) / ring_sizes_[1])), num_ring_1 - 1);
        sector_idx = min(static_cast<int>((theta / sector_sizes_[1])), num_sector_1 - 1);
        czm[1][ring_idx][sector_idx].emplace_back(PointXYZ(x, y, z, i));
      } else if (r < min_range_3) {
        ring_idx   = min(static_cast<int>(((r - min_range_2) / ring_sizes_[2])), num_ring_2 - 1);
        sector_idx = min(static_cast<int>((theta / sector_sizes_[2])), num_sector_2 - 1);
        czm[2][ring_idx][sector_idx].emplace_back(PointXYZ(x, y, z, i));
      } else {  // Far!
        ring_idx   = min(static_cast<int>(((r - min_range_3) / ring_sizes_[3])), num_ring_3 - 1);
        sector_idx = min(static_cast<int>((theta / sector_sizes_[3])), num_sector_3 - 1);
        czm[3][ring_idx][sector_idx].emplace_back(PointXYZ(x, y, z, i));
      }

    } else {
      cloud_nonground_.push_back(PointXYZ(x, y, z, i));
    }
  }
  if (params_.verbose)
    cout << "\033[1;33m"
         << "PatchWorkpp::pc2czm() - Divides pointcloud into the concentric zone model successfully"
         << "\033[0m" << endl;
}

void PatchWorkpp::analyze_traversability() {
  cloud_traversable_ground_.clear();
  cloud_nontraversable_ground_.clear();
  terrain_normals_.clear();

  if (cloud_ground_.empty()) return;

  std::unordered_map<int, std::vector<PointXYZ>> grid;
  
  double min_x = 100000, max_x = -100000;
  double min_y = 100000, max_y = -100000;
  for (const auto& p : cloud_ground_) {
    if (p.x < min_x) min_x = p.x;
    if (p.x > max_x) max_x = p.x;
    if (p.y < min_y) min_y = p.y;
    if (p.y > max_y) max_y = p.y;
  }

  int grid_width = std::ceil((max_x - min_x) / params_.terrain_patch_size) + 1;

  for (const auto& p : cloud_ground_) {
    int x_idx = std::floor((p.x - min_x) / params_.terrain_patch_size);
    int y_idx = std::floor((p.y - min_y) / params_.terrain_patch_size);
    int idx = y_idx * grid_width + x_idx;
    grid[idx].push_back(p);
  }

  for (const auto& pair : grid) {
    const auto& points = pair.second;
    if (points.size() < static_cast<size_t>(params_.min_patch_points)) {
      addCloud(cloud_traversable_ground_, points);
      continue;
    }

    const size_t n = points.size();
    double sx = 0, sy = 0, sz = 0;
    double sxx = 0, syy = 0, szz = 0;
    double sxy = 0, sxz = 0, syz = 0;
    double min_z = 100000, max_z = -100000;
    for (const auto& p : points) {
      sx += p.x; sy += p.y; sz += p.z;
      sxx += p.x * p.x; syy += p.y * p.y; szz += p.z * p.z;
      sxy += p.x * p.y; sxz += p.x * p.z; syz += p.y * p.z;
      if (p.z < min_z) min_z = p.z;
      if (p.z > max_z) max_z = p.z;
    }
    double inv_n = 1.0 / n;
    double mx = sx * inv_n, my = sy * inv_n, mz = sz * inv_n;
    double denom = (n > 1) ? static_cast<double>(n - 1) : 1.0;
    double inv_d = 1.0 / denom;

    Eigen::Matrix3f cov;
    cov(0, 0) = (sxx - n * mx * mx) * inv_d;
    cov(1, 1) = (syy - n * my * my) * inv_d;
    cov(2, 2) = (szz - n * mz * mz) * inv_d;
    cov(0, 1) = cov(1, 0) = (sxy - n * mx * my) * inv_d;
    cov(0, 2) = cov(2, 0) = (sxz - n * mx * mz) * inv_d;
    cov(1, 2) = cov(2, 1) = (syz - n * my * mz) * inv_d;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eig(cov, Eigen::ComputeEigenvectors);
    Eigen::Vector3f evals = eig.eigenvalues().cwiseMax(0.0f);
    Eigen::Vector3f normal = eig.eigenvectors().col(0);
    if (normal(2) < 0.0f) normal = -normal;

    double roughness = evals(0); 
    double slope_deg = std::acos(std::min(1.0f, std::max(-1.0f, normal(2)))) * 180.0 / M_PI;
    double height_variance = max_z - min_z;

    if (params_.publish_terrain_normals) {
      terrain_normals_.push_back(PointXYZ(normal(0), normal(1), normal(2), 0));
    }

    bool is_traversable = true;
    if (slope_deg > params_.max_traversable_slope_deg) is_traversable = false;
    if (roughness > params_.roughness_threshold) is_traversable = false;
    if (height_variance > params_.vegetation_threshold) is_traversable = false;

    if (is_traversable) {
      addCloud(cloud_traversable_ground_, points);
    } else {
      addCloud(cloud_nontraversable_ground_, points);
    }
  }
}
