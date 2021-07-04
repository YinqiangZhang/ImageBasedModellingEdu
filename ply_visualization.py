import os
import numpy as np
import open3d as o3d

filename = "points.ply"
filepath = os.path.join(os.getcwd(), filename)
pcd = o3d.io.read_point_cloud(filepath)
np_pcd = np.asarray(pcd.points)
print(f"total number of points: {len(np_pcd)}")
o3d.visualization.draw_geometries([pcd])