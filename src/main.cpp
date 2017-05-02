#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/progressive_morphological_filter.h>

#include <boost/thread/thread.hpp>

#include <pcl/common/common_headers.h>
#include <pcl/range_image/range_image.h>
#include <pcl/visualization/range_image_visualizer.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include<pcl/io/ply_io.h>


#include <pcl/point_types.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/progressive_morphological_filter.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/visualization/point_cloud_handlers.h>

#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>

typedef pcl::PointXYZ PointType;

// --------------------
// -----Parameters-----
// --------------------
float angular_resolution_x = 0.5f,
      angular_resolution_y = angular_resolution_x;
pcl::RangeImage::CoordinateFrame coordinate_frame = pcl::RangeImage::CAMERA_FRAME;
bool live_update = false;

// --------------
// -----Help-----
// --------------
void printUsage (const char* progName){
  std::cout << "\n\nUsage: "<<progName<<" [options] <scene.pcd>\n\n"
            << "Options:\n"
            << "-------------------------------------------\n"
            << "-rx <float>  angular resolution in degrees (default "<<angular_resolution_x<<")\n"
            << "-ry <float>  angular resolution in degrees (default "<<angular_resolution_y<<")\n"
            << "-c <int>     coordinate frame (default "<< (int)coordinate_frame<<")\n"
            << "-l           live update - update the range image according to the selected view in the 3D viewer.\n"
            << "-h           this help\n"
            << "\n\n";
}

void setViewerPose (pcl::visualization::PCLVisualizer& viewer, const Eigen::Affine3f& viewer_pose){
  Eigen::Vector3f pos_vector = viewer_pose * Eigen::Vector3f(0, 0, 0);
  Eigen::Vector3f look_at_vector = viewer_pose.rotation () * Eigen::Vector3f(0, 0, 1) + pos_vector;
  Eigen::Vector3f up_vector = viewer_pose.rotation () * Eigen::Vector3f(0, -1, 0);
  viewer.setCameraPosition (pos_vector[0], pos_vector[1], pos_vector[2],
                            look_at_vector[0], look_at_vector[1], look_at_vector[2],
                            up_vector[0], up_vector[1], up_vector[2]);
}

// --------------
// -----Main-----
// --------------
int main (int argc, char** argv){
  // --------------------------------------
  // -----Parse Command Line Arguments-----
  // --------------------------------------
  if (pcl::console::find_argument (argc, argv, "-h") >= 0)
  {
    printUsage (argv[0]);
    return 0;
  }
  if (pcl::console::find_argument (argc, argv, "-l") >= 0)
  {
    live_update = true;
    std::cout << "Live update is on.\n";
  }
  if (pcl::console::parse (argc, argv, "-rx", angular_resolution_x) >= 0)
    std::cout << "Setting angular resolution in x-direction to "<<angular_resolution_x<<"deg.\n";
  if (pcl::console::parse (argc, argv, "-ry", angular_resolution_y) >= 0)
    std::cout << "Setting angular resolution in y-direction to "<<angular_resolution_y<<"deg.\n";
  int tmp_coordinate_frame;
  if (pcl::console::parse (argc, argv, "-c", tmp_coordinate_frame) >= 0)
  {
    coordinate_frame = pcl::RangeImage::CoordinateFrame (tmp_coordinate_frame);
    std::cout << "Using coordinate frame "<< (int)coordinate_frame<<".\n";
  }
  angular_resolution_x = pcl::deg2rad (angular_resolution_x);
  angular_resolution_y = pcl::deg2rad (angular_resolution_y);

  // ------------------------------------------------------------------
  // -----Read pcd file or create example point cloud if not given-----
  // ------------------------------------------------------------------
  pcl::PointCloud<PointType>::Ptr point_cloud_ptr (new pcl::PointCloud<PointType>);
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::PointCloud<PointType>& point_cloud = *point_cloud_ptr;
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_f (new pcl::PointCloud<pcl::PointXYZ>);


  Eigen::Affine3f scene_sensor_pose (Eigen::Affine3f::Identity ());
  std::vector<int> pcd_filename_indices = pcl::console::parse_file_extension_argument (argc, argv, "ply");
	if (!pcd_filename_indices.empty ()){
		std::string filename = argv[pcd_filename_indices[0]];



		pcl::PLYReader Reader;
		std::cout << filename << '\n';
		Reader.read(filename, point_cloud);

		// if (pcl::io::loadPCDFile (filename, point_cloud) == -1){
		// 	std::cout << "Was not able to open file \""<<filename<<"\".\n";
		// 	printUsage (argv[0]);
		// 	return 0;
		// }

		scene_sensor_pose = Eigen::Affine3f (Eigen::Translation3f (point_cloud.sensor_origin_[0], point_cloud.sensor_origin_[1], point_cloud.sensor_origin_[2])) *Eigen::Affine3f (point_cloud.sensor_orientation_);
	}
	else{
		std::cout << "\nNo *.pcd file given => Genarating example point cloud.\n\n";
		for (float x=-0.5f; x<=0.5f; x+=0.01f){
			for (float y=-0.5f; y<=0.5f; y+=0.01f){
				PointType point;  point.x = x;  point.y = y;  point.z = 2.0f - y;
				point_cloud.points.push_back (point);
			}
		}
		point_cloud.width = (int) point_cloud.points.size ();  point_cloud.height = 1;
	}


    std::cout << "all that filtering bullshit" << '\n';


    std::cout << "all that filtering bullshit ENDS" << '\n';

  // -----------------------------------------------
  // -----Create RangeImage from the PointCloud-----
  // -----------------------------------------------
  float noise_level = 0.0;
  float min_range = 0.0f;
  int border_size = 1;
  boost::shared_ptr<pcl::RangeImage> range_image_ptr(new pcl::RangeImage);
  pcl::RangeImage& range_image = *range_image_ptr;
  range_image.createFromPointCloud (point_cloud, angular_resolution_x, angular_resolution_y,
                                    pcl::deg2rad (360.0f), pcl::deg2rad (180.0f),
                                    scene_sensor_pose, coordinate_frame, noise_level, min_range, border_size);

  // --------------------------------------------
  // -----Open 3D viewer and add point cloud-----
  // --------------------------------------------
  pcl::visualization::PCLVisualizer viewer ("3D Viewer");
  viewer.setBackgroundColor (1, 1, 1);
  // pcl::visualization::PointCloudColorHandlerCustom<pcl::PointWithRange> range_image_color_handler (range_image_ptr, 0, 0, 0);
  // viewer.addPointCloud (range_image_ptr, range_image_color_handler, "range image");
  viewer.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "range image");

  viewer.addCoordinateSystem (1.0f, "global");

  pcl::visualization::PointCloudColorHandlerCustom<PointType> point_cloud_color_handler (point_cloud_ptr, 150, 0, 0);
  // viewer.addPointCloud (point_cloud_ptr, point_cloud_color_handler, "original point cloud");

    std::cerr << "Cloud before filtering: " << std::endl;
    std::cerr << point_cloud << std::endl;

    pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
    sor.setInputCloud (point_cloud_ptr);
    sor.setMeanK (50);
    sor.setStddevMulThresh (1.0);
    sor.filter (*cloud_filtered);

    std::cerr << "Cloud after filtering: " << std::endl;
    std::cerr << *cloud_filtered << std::endl;

    // point_cloud_color_handler (cloud_filtered, 0, 0, 150);

    // Create the filtering object: downsample the dataset using a leaf size of 1cm
    pcl::VoxelGrid<pcl::PointXYZ> vg;
    vg.setInputCloud (cloud_filtered);
    vg.setLeafSize (0.01f, 0.01f, 0.01f);
    vg.filter (*cloud_filtered);
    std::cout << "PointCloud after filtering has: " << cloud_filtered->points.size ()  << " data points." << std::endl; //*
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color(cloud_filtered, 0, 255, 0);
    viewer.addPointCloud(cloud_filtered, single_color, "filtered point cloud");

    // Create the segmentation object for the planar model and set all the parameters
    pcl::SACSegmentation<pcl::PointXYZ> seg;
    pcl::PointIndices::Ptr inliers (new pcl::PointIndices);
    pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_plane (new pcl::PointCloud<pcl::PointXYZ> ());
    pcl::PCDWriter writer;
    seg.setOptimizeCoefficients (true);
    seg.setModelType (pcl::SACMODEL_PLANE);
    seg.setMethodType (pcl::SAC_RANSAC);
    seg.setMaxIterations (100);
    seg.setDistanceThreshold (0.02);


    int i=0, nr_points = (int) cloud_filtered->points.size ();
    while (cloud_filtered->points.size () > 0.3 * nr_points){
        // Segment the largest planar component from the remaining cloud
        seg.setInputCloud (cloud_filtered);
        seg.segment (*inliers, *coefficients);
        if (inliers->indices.size () == 0){
            std::cout << "Could not estimate a planar model for the given dataset." << std::endl;
            break;
        }

        // Extract the planar inliers from the input cloud
        pcl::ExtractIndices<pcl::PointXYZ> extract;
        extract.setInputCloud (cloud_filtered);
        extract.setIndices (inliers);
        extract.setNegative (false);

        // Get the points associated with the planar surface
        extract.filter (*cloud_plane);
        std::cout << "PointCloud representing the planar component: " << cloud_plane->points.size () << " data points." << std::endl;
        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color(cloud_filtered, 0, 0, 255);
        viewer.addPointCloud(cloud_plane, single_color, "plane point cloud");

        // Remove the planar inliers, extract the rest
        extract.setNegative (true);
        extract.filter (*cloud_f);
        *cloud_filtered = *cloud_f;
    }








  viewer.initCameraParameters ();
  setViewerPose(viewer, range_image.getTransformationToWorldSystem ());

  // --------------------------
  // -----Show range image-----
  // --------------------------
  pcl::visualization::RangeImageVisualizer range_image_widget ("Range image");
  range_image_widget.showRangeImage (range_image);

  //--------------------
  // -----Main loop-----
  //--------------------
  while (!viewer.wasStopped ())
  {
    range_image_widget.spinOnce ();
    viewer.spinOnce ();
    pcl_sleep (0.01);

    if (live_update)
    {
      scene_sensor_pose = viewer.getViewerPose();
      range_image.createFromPointCloud (point_cloud, angular_resolution_x, angular_resolution_y,
                                        pcl::deg2rad (360.0f), pcl::deg2rad (180.0f),
                                        scene_sensor_pose, pcl::RangeImage::LASER_FRAME, noise_level, min_range, border_size);
      range_image_widget.showRangeImage (range_image);
    }
  }
}
