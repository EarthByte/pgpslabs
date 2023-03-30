import geopandas as gpd
import cartopy.crs as ccrs
import os
import shutil
import matplotlib.pyplot as plt
import matplotlib as mpl
import pygplates
import numpy as np
import ptt  # plate tectonic tools
import shapely
from shapely.geometry import LineString
import warnings
import moviepy.editor as mpy

import slab_tracker_utils as slab
import splits_and_merges as snm

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

#######################################################
# Define Input files
#######################################################

plate_model_dir = 'C:/Users/Grace/Documents/SlabWindow/Global_Model_WD_Internal_Release_2022_v2'

RotFile_List = [
    '%s/Alps_Mesh_Rotations.rot' % plate_model_dir,
    '%s/Andes_Flat_Slabs_Rotations.rot' % plate_model_dir,
    '%s/Andes_Rotations.rot' % plate_model_dir,
    '%s/Australia_Antarctica_Mesh_Rotations.rot' % plate_model_dir,
    '%s/Australia_North_Zealandia_Rotations.rot' % plate_model_dir,
    '%s/Eurasia_Arabia_Mesh_Rotations.rot' % plate_model_dir,
    '%s/Global_250-0Ma_Rotations.rot' % plate_model_dir,
    '%s/Global_410-250Ma_Rotations.rot' % plate_model_dir,
    '%s/North_America_Flat_Slabs_Rotations.rot' % plate_model_dir,
    '%s/North_America_Mesh_Rotations.rot' % plate_model_dir,
    '%s/North_China_Mesh_Rotations.rot' % plate_model_dir,
    '%s/South_Atlantic_Rotations.rot' % plate_model_dir,
    '%s/South_China_DeformingModel.rot' % plate_model_dir,
    '%s/Southeast_Asia_Rotations.rot' % plate_model_dir]

GPML_List = [
    '%s/Alps_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Alps_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/America_Anyui_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/America_Anyui_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Andes_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Andes_Flat_Slabs_Topologies.gpml' % plate_model_dir,
    '%s/Andes_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Arctic_Eurasia_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Arctic_Eurasia_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Australia_Antarctica_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Australia_Antarctica_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Australia_North_Zealandia_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Australia_North_Zealandia_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Baja_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Coral_Sea_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Coral_Sea_Topologies.gpml' % plate_model_dir,
    '%s/East_African_Rift_Deforming_Mesh_and_Topologies.gpml' % plate_model_dir,
    '%s/Eastern_China_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/East-West_Gondwana_Deforming_Mesh_and_Topologies.gpml' % plate_model_dir,
    '%s/Ellesmere_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Eurasia_Arabia_Deforming_Mesh_and_Topologies.gpml' % plate_model_dir,
    '%s/Global_Mesozoic-Cenozoic_PlateBoundaries.gpml' % plate_model_dir,
    '%s/Global_Paleozoic_PlateBoundaries.gpml' % plate_model_dir,
    '%s/Greater_India_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Greater_India_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Inactive_Meshes_and_Topologies.gpml' % plate_model_dir,
    '%s/North_America_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/North_Atlantic_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/North_Atlantic_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/North_China_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Northern_Andes_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Northern_Andes_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Papua_New_Guinea_Deforming_Meshes.gpml' % plate_model_dir,
    '%s/Papua_New_Guinea_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Scotia_Deforming_Mesh_and_Topologies.gpml' % plate_model_dir,
    '%s/Siberia_Eurasia_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Siberia_Eurasia_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/South_Atlantic_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/South_Atlantic_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/South_China_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/South_China_DeformingElements.gpml' % plate_model_dir,
    '%s/South_Zealandia_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/South_Zealandia_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Southeast_Asia_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Southeast_Asia_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/West_Antarctic_Zealandia_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/West_Antarctica_Zealandia_Mesh_Topologies.gpml' % plate_model_dir,
    '%s/Western_North_America_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Western_Tethys_Deforming_Mesh.gpml' % plate_model_dir,
    '%s/Western_Tethys_Tectonic_Boundary_Topologies.gpml' % plate_model_dir]

coastline_filename = 'C:/Users/Grace/Documents/SlabWindow/Global_Model_WD_Internal_Release_2022_v2/StaticGeometries/Coastlines/Global_coastlines_low_res.shp'
coastlines = pygplates.FeatureCollection(coastline_filename)

# OPTIONAL
# specify the convention for agegrids
# assumes that the file name has 1 number within it that specifies the age in Myr
# the list object below must contain 2 strings, which specify the bit before and after this number
# agegrid_filename = ['/Users/Simon/Data/AgeGrids/Agegrids_30m_20151002_2015_v1_r756/agegrid_30m_','.grd']

# If you do not have age grids, use the next line instead:
agegrid_filename = None

#####################################
rotation_model = pygplates.RotationModel(RotFile_List)
topology_features = pygplates.FeatureCollection()
for file in GPML_List:
    topology_feature = pygplates.FeatureCollection(file)
    topology_features.add(topology_feature)

#############################
# INPUT PARAMETERS

start_time = 45
end_time = 0
time_step = 1.0
dip_angle_degrees = 45.0
line_tessellation_distance = np.radians(1.0)
reconstruction_time = end_time
bin_time = 15  # how many rows of slab will be displayed on the map at any time
global_map = True  # set true if producing global map, false to set region by co-ordinates
mp4 = True  # set true to produce a .mp4 file, set false for a .gif file.
handle_splits = False

point_size = 1  # size of individual points denoting slab (I find 01 good for global maps, and ~40 good for regional ones)
fps = 2  # number of frames per second in the final  animation
vmin = 0
vmax = 660  # use this to control the scale over which the colours denoting the depth of the subducting slab change

# Define the scope of your map (used only if global_map=False)
X_start = 10  # beginning of x-axis
X_end = 100  # end of x-axis
Y_start = 10  # start of y-axis
Y_end = 80  # end of y-axis

# By way of example the above values will be used by xlim ylim as ax.set_xlim(10, 100) and # ax.set_ylim(10, 80)

###  Potential sites of interest
# Central America
# ax.set_xlim(75, 105) #left/right
# ax.set_ylim(5, 25) # bottom/top

# North America
# ax.set_xlim(10, 100) #left/right
# ax.set_ylim(10, 80) # bottom/top

# Chile
# ax.set_xlim(90, 150)
# ax.set_ylim(-60, 0)

# Try to use small circle path for stage rotation to rotate along velocity dip.
# Ie, distance to stage rotation pole matches distance to original stage pole.
# use_small_circle_path = False

coord_handler = ccrs.PlateCarree()

output_filename = 'subduction_3d_geometries_time_%0.2fMa_dip_%0.2fdeg.asc' % (reconstruction_time, dip_angle_degrees)

#############################
# Create Figure

fig = plt.figure(figsize=(20, 15))
if global_map:
    ax = fig.add_subplot(111, projection=ccrs.Robinson(central_longitude=180))
    ax.set_global
else:
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree(central_longitude=180))

# array into which results will be stored

output_data = []

dip_angle_radians = np.radians(dip_angle_degrees)

# Main Warping Function
def subducting_slab():

    if handle_splits:
        plate_disappearance_time_lut = snm.get_plate_disappearance_time_lut(topology_features,
                                                                            rotation_model,
                                                                            time_list,
                                                                            verbose=True)

        print(plate_disappearance_time_lut)

    # loop over a series of times at which we want to extract trench iso-sub-chrons
    for time in time_list:

        print('time %0.2f Ma' % time)

        # call function to get subduction boundary segments
        subduction_boundary_sections = slab.getSubductionBoundarySections(topology_features,
                                                                          rotation_model,
                                                                          time)

        # Set up an age grid interpolator for this time, to be used
        # for each tessellated line segment
        if agegrid_filename is not None:
            grdfile = '%s%d%s' % (agegrid_filename[0], time, agegrid_filename[1])
            lut = slab.make_age_interpolator(grdfile)

        # Loop over each segment
        for segment_index, subduction_segment in enumerate(subduction_boundary_sections):

            # find the overriding plate id (and only continue if we find it)
            overriding_and_subducting_plates = slab.find_subducting_plate(subduction_segment, False)
            # overriding_and_subducting_plates = slab.find_overriding_and_subducting_plates(subduction_segment, time)

            if not overriding_and_subducting_plates:
                continue
            overriding_plate, subducting_plate, subduction_polarity = overriding_and_subducting_plates

            overriding_plate_id = overriding_plate.get_resolved_feature().get_reconstruction_plate_id()
            subducting_plate_id = subducting_plate.get_resolved_feature().get_reconstruction_plate_id()

            # if (opid!=224 or cpid!=909):
            # if (subducting_plate_id!=911 and subducting_plate_id!=909):
            # if subducting_plate_id<900:
            #    continue

            subducting_plate_disappearance_time = -1.
            if handle_splits:
                for plate_disappearance in plate_disappearance_time_lut:
                    if plate_disappearance[0] == subducting_plate_id:
                        subducting_plate_disappearance_time = plate_disappearance[1]

            tessellated_line = subduction_segment.get_resolved_geometry().to_tessellated(line_tessellation_distance)

            # print(len(tessellated_line.get_points()))

            if agegrid_filename is not None:
                x = tessellated_line.to_lat_lon_array()[:, 1]
                y = tessellated_line.to_lat_lon_array()[:, 0]
                subduction_ages = lut.ev(np.radians(y + 90.), np.radians(x + 180.))
            else:
                # if no age grids, just fill the ages with zero
                subduction_ages = [0. for point in tessellated_line.to_lat_lon_array()[:, 1]]

            # CALL THE MAIN WARPING FUNCTION
            (points,
             point_depths,
             polyline) = slab.warp_subduction_segment(tessellated_line,
                                                      rotation_model,
                                                      subducting_plate_id,
                                                      overriding_plate_id,
                                                      subduction_polarity,
                                                      time,
                                                      end_time,
                                                      time_step,
                                                      dip_angle_radians,
                                                      subducting_plate_disappearance_time)

            # add points for this segment to the plot

            ax.scatter(polyline.to_lat_lon_array()[:, 1],
                       polyline.to_lat_lon_array()[:, 0],
                       c=np.array(point_depths), s=point_size,
                       edgecolors=None, zorder=2, cmap=plt.cm.gnuplot_r, vmin=vmin, vmax=vmax,
                       transform=coord_handler)
            # ax.plot(polyline.to_lat_lon_array()[:, 1], polyline.to_lat_lon_array()[:, 0],
            #       '-r', markersize=1, zorder=1)

            output_data.append([time, polyline, point_depths, subduction_ages])


ax.set_aspect('auto')


##################################
#       Map and Animations       #
##################################

### Create some functions to manipulate points and topologies

def create_geodataframe_general(pygplates_recon_geom, reconstruction_time):
    """ This is a general function to convert reconstructed features (e.g.
    reconstructed coastlines) from pygplates into a GeoDataFrame. This helps
    avoid plotting artefacts. Note that the input geometry must be a polygon.

    Input:
        - pygplates.ReconstructedFeatureGeometry (i.e., output of pygplates.reconstruct)
        - recontruction time - this is just for safekeeping in the geodataframe!
    Output:
        - gpd.GeoDataFrame of the feature"""

    # create new and empty geodataframe
    recon_gpd = gpd.GeoDataFrame()
    recon_gpd['NAME'] = None
    recon_gpd['PLATEID1'] = None
    recon_gpd['PLATEID2'] = None
    recon_gpd['FROMAGE'] = None
    recon_gpd['TOAGE'] = None
    recon_gpd['geometry'] = None
    recon_gpd['reconstruction_time'] = None
    recon_gpd = recon_gpd.set_crs(epsg=4326)

    date_line_wrapper = pygplates.DateLineWrapper()

    names = []
    plateid1s = []
    plateid2s = []
    fromages = []
    toages = []
    geometrys = []
    reconstruction_times = []

    for i, seg in enumerate(pygplates_recon_geom):

        wrapped_polygons = date_line_wrapper.wrap(seg.get_reconstructed_geometry())
        for poly in wrapped_polygons:
            ring = np.array([(p.get_longitude(), p.get_latitude()) for p in poly.get_exterior_points()])
            ring[:, 1] = np.clip(ring[:, 1], -89, 89)  # anything approaching the poles creates artefacts

            name = seg.get_feature().get_name()
            plateid = seg.get_feature().get_reconstruction_plate_id()
            conjid = seg.get_feature().get_conjugate_plate_id()
            from_age, to_age = seg.get_feature().get_valid_time()
            # append things
            names.append(name)
            plateid1s.append(plateid)
            plateid2s.append(conjid)
            fromages.append(from_age)
            toages.append(to_age)
            geometrys.append(shapely.geometry.Polygon(ring))
            reconstruction_times.append(reconstruction_time)

    # write to geodataframe
    recon_gpd['NAME'] = names
    recon_gpd['PLATEID1'] = plateid1s
    recon_gpd['PLATEID2'] = plateid2s
    recon_gpd['FROMAGE'] = fromages
    recon_gpd['TOAGE'] = toages
    recon_gpd['geometry'] = geometrys
    recon_gpd['reconstruction_time'] = reconstruction_times

    return recon_gpd


def create_geodataframe_topologies(topologies, reconstruction_time):
    """ This is a function to convert topologies from pygplates into a GeoDataFrame
    This helps select the closed topological plates ('gpml:TopologicalClosedPlateBoundary',
    and also helps resolve plotting artefacts from crossing the dateline.
    This function does NOT incorporate various plate boundary types into the geodataframe!

    Input:
        - pygplates.Feature. This is designed for `topologies`, which comes from:
              resolved_topologies = ptt.resolve_topologies.resolve_topologies_into_features(
                                        rotation_model, topology_features, reconstruction_time)
              topologies, ridge_transforms, ridges, transforms, trenches, trench_left, trench_right, other = resolved_topologies
        - recontruction time - this is just for safekeeping in the geodataframe!
    Output:
        - gpd.GeoDataFrame of the feature"""

    # function for getting closed topologies only
    # i.e., the plates themselves, NOT all the features for plotting!

    # set up the empty geodataframe
    recon_gpd = gpd.GeoDataFrame()
    recon_gpd['NAME'] = None
    recon_gpd['PLATEID1'] = None
    recon_gpd['PLATEID2'] = None
    recon_gpd['FROMAGE'] = None
    recon_gpd['TOAGE'] = None
    recon_gpd['geometry'] = None
    recon_gpd['reconstruction_time'] = None
    recon_gpd['gpml_type'] = None
    recon_gpd = recon_gpd.set_crs(epsg=4326)

    # some empty things to write stuff to
    names = []
    plateid1s = []
    plateid2s = []
    fromages = []
    toages = []
    geometrys = []
    reconstruction_times = []
    gpml_types = []

    # a dateline wrapper! so that they plot nicely and do nice things in geopandas
    date_line_wrapper = pygplates.DateLineWrapper()

    for i, seg in enumerate(topologies):
        gpmltype = seg.get_feature_type()

        # polygon and wrap
        polygon = seg.get_geometry()
        wrapped_polygons = date_line_wrapper.wrap(polygon)
        for poly in wrapped_polygons:
            ring = np.array([(p.get_longitude(), p.get_latitude()) for p in poly.get_exterior_points()])
            ring[:, 1] = np.clip(ring[:, 1], -89, 89)  # anything approaching the poles creates artefacts
            for wrapped_point in poly.get_exterior_points():
                wrapped_point_lat_lon = wrapped_point.get_latitude(), wrapped_point.get_longitude()

            # might result in two polys - append to loop here (otherwise we will be missing half the pacific etc)
            name = seg.get_name()
            plateid = seg.get_reconstruction_plate_id()
            conjid = seg.get_conjugate_plate_id()
            from_age, to_age = seg.get_valid_time()

            names.append(name)
            plateid1s.append(plateid)
            plateid2s.append(conjid)
            fromages.append(from_age)
            toages.append(to_age)
            geometrys.append(shapely.geometry.Polygon(ring))
            reconstruction_times.append(reconstruction_time)
            gpml_types.append(str(gpmltype))

    # write to geodataframe
    recon_gpd['NAME'] = names
    recon_gpd['PLATEID1'] = plateid1s
    recon_gpd['PLATEID2'] = plateid2s
    recon_gpd['FROMAGE'] = fromages
    recon_gpd['TOAGE'] = toages
    recon_gpd['geometry'] = geometrys
    recon_gpd['reconstruction_time'] = reconstruction_times
    recon_gpd['gpml_type'] = gpml_types

    return recon_gpd


def create_geodataframe_lines(pygplates_feature, reconstruction_time, flip_geom=False):
    """ This is a general function to convert reconstructed LINE features (e.g.
    reconstructed mid ocean ridges) from pygplates into a GeoDataFrame. This helps
    avoid plotting artefacts. Note that the input geometry must be a line geometry.

    Input:
        - pygplates.ReconstructedFeatureGeometry (i.e., output of pygplates.reconstruct)
        - recontruction time - this is just for safekeeping in the geodataframe!
    Output:
        - gpd.GeoDataFrame of the feature

    To make life easier for trenches, there is an additional parameter 'flip_geom'.
    Use this to invert the order of either the left OR right trenches, so they can be plotted (e.g. using pygmt) with a single command.

    For example:
    gpd_trench_left = create_geodataframe_lines(trench_left, reconstruction_time, flip_geom=True)
    gpd_trench_right = create_geodataframe_lines(trench_right, reconstruction_time)
    gpd_trenches = gpd_trench_left.append(gpd_trench_right)

    can then be plotted in pygmt using:
    fig.plot(gpd_trenches, pen='2p,black', style='f0.3/0.1+l+t')

    """

    # some empty things to write stuff to
    names = []
    plateid1s = []
    plateid2s = []
    left_plates = []
    right_plates = []
    fromages = []
    toages = []
    geometrys = []
    reconstruction_times = []
    gpml_types = []

    # set up the empty geodataframe
    recon_gpd = gpd.GeoDataFrame()
    recon_gpd['NAME'] = None
    recon_gpd['PLATEID1'] = None
    recon_gpd['PLATEID2'] = None
    recon_gpd['left_plates'] = None
    recon_gpd['right_plates'] = None
    recon_gpd['FROMAGE'] = None
    recon_gpd['TOAGE'] = None
    recon_gpd['geometry'] = None
    recon_gpd['reconstruction_time'] = None
    recon_gpd['gpml_type'] = None
    recon_gpd = recon_gpd.set_crs(epsg=4326)

    #    date_line_wrapper = pygplates.DateLineWrapper()

    # loop through
    for feature in pygplates_feature:

        # Add logic to Cartopy to connect points across the date line instead of going around the back of the world.
        geometry = feature.get_all_geometries()[0].to_lat_lon_array()[::-1, ::-1]

        sanitize_longs = []
        cross_date_line = False  # Create flag to define whether LineSting goes across the date line
        for ix, ea in enumerate(geometry):
            if ix == 0:
                old_long = ea[0]  # if long/lat index = 0, then leave it alone.
            else:  # For every long/lat tuple other than 0
                difflong = old_long - ea[0]
                if abs(difflong) > 180:
                    cross_date_line = True  # Flag as crossing dateline
                old_long = ea[0]

        for ix, ea in enumerate(geometry):
            if cross_date_line:
                if ea[0] < 0:
                    ea[0] = ea[0] + 360  # If the long/lat tuple crosses the date line, add 360 to the longitude.
                sanitize_longs.append(
                    ea)  # append the longitude, repeat the process for all relevant linestring points.

        if cross_date_line:
            geometry = sanitize_longs

        if flip_geom is True:
            geometries = np.flip(geometry, axis=0)
        else:
            geometries = geometry

        # construct shapely geometry
        geom = shapely.geometry.LineString(geometries)

        names.append(feature.get_name())
        left_plates.append(feature.get_left_plate())
        right_plates.append(feature.get_right_plate())
        plateid1s.append(feature.get_reconstruction_plate_id())
        plateid2s.append(feature.get_conjugate_plate_id())
        fromages.append(feature.get_valid_time()[0])
        toages.append(feature.get_valid_time()[1])
        reconstruction_times.append(reconstruction_time)
        gpml_types.append(str(feature.get_feature_type()))

        if geom.is_valid:
            geometrys.append(geom)

    # write to geodataframe
    recon_gpd['NAME'] = names
    recon_gpd['PLATEID1'] = plateid1s
    recon_gpd['PLATEID2'] = plateid2s
    recon_gpd['FROMAGE'] = fromages
    recon_gpd['TOAGE'] = toages
    recon_gpd['geometry'] = geometrys
    recon_gpd['reconstruction_time'] = reconstruction_times
    recon_gpd['gpml_type'] = gpml_types
    recon_gpd['left_plates'] = left_plates
    recon_gpd['right_plates'] = right_plates

    return recon_gpd


def create_geodataframe_points(pygplates_feature, reconstruction_time, flip_geom=False):
    """ This is a general function to convert reconstructed point/multipoint features
    (i.e. any kind of point data, e.g. hotspots, magpicks, fossils) from pygplates into
    a GeoDataFrame. Note that the input geometry must be a multipoint/point geometry.

    Input:
        - pygplates.ReconstructedFeatureGeometry (i.e., output of pygplates.reconstruct)
        - recontruction time - this is just for safekeeping in the geodataframe!
    Output:
        - gpd.GeoDataFrame of the feature

    """

    # some empty things to write stuff to
    names = []
    plateid1s = []
    plateid2s = []
    fromages = []
    toages = []
    geometrys = []
    reconstruction_times = []
    gpml_types = []

    # set up the empty geodataframe
    recon_gpd = gpd.GeoDataFrame()
    recon_gpd['NAME'] = None
    recon_gpd['PLATEID1'] = None
    recon_gpd['PLATEID2'] = None
    recon_gpd['FROMAGE'] = None
    recon_gpd['TOAGE'] = None
    recon_gpd['geometry'] = None
    recon_gpd['reconstruction_time'] = None
    recon_gpd['gpml_type'] = None
    recon_gpd = recon_gpd.set_crs(epsg=4326)

    # loop through
    for feature in pygplates_feature:
        # feature is a 'pygplates.ReconstructedFeatureGeometry'
        rlon = feature.get_reconstructed_geometry().to_lat_lon_array()[0][1]
        rlat = feature.get_reconstructed_geometry().to_lat_lon_array()[0][0]

        # construct shapely geometry
        geom = shapely.geometry.Point(rlon, rlat)

        # need to get_feature() first before we can get some of the various properties
        names.append(feature.get_feature().get_name())
        plateid1s.append(feature.get_feature().get_reconstruction_plate_id())
        plateid2s.append(feature.get_feature().get_conjugate_plate_id())
        fromages.append(feature.get_feature().get_valid_time()[0])
        toages.append(feature.get_feature().get_valid_time()[1])
        reconstruction_times.append(reconstruction_time)
        gpml_types.append(str(feature.get_feature().get_feature_type()))

        if geom.is_valid:
            geometrys.append(geom)

    # write to geodataframe
    recon_gpd['NAME'] = names
    recon_gpd['PLATEID1'] = plateid1s
    recon_gpd['PLATEID2'] = plateid2s
    recon_gpd['FROMAGE'] = fromages
    recon_gpd['TOAGE'] = toages
    recon_gpd['geometry'] = geometrys
    recon_gpd['reconstruction_time'] = reconstruction_times
    recon_gpd['gpml_type'] = gpml_types

    return recon_gpd


### Velocity Arrows

def make_GPML_velocity_feature(Lon, Lat):
    # function to make a velocity mesh nodes at an arbitrary set of points defined in Lat
    # Lon and Lat are assumed to be 1d arrays.

    # Add points to a multipoint geometry
    multi_point = pygplates.MultiPointOnSphere([(float(lat), float(lon)) for lat, lon in zip(Lat, Lon)])

    # Create a feature containing the multipoint feature, and defined as MeshNode type
    meshnode_feature = pygplates.Feature(pygplates.FeatureType.create_from_qualified_string('gpml:MeshNode'))
    meshnode_feature.set_geometry(multi_point)
    meshnode_feature.set_name('Velocity Mesh Nodes from pygplates')

    output_feature_collection = pygplates.FeatureCollection(meshnode_feature)
    return output_feature_collection


# Generate points for regular long,lat grid
Xnodes = np.arange(-180, 180, 20)
Ynodes = np.arange(-90, 90, 20)
Xg, Yg = np.meshgrid(Xnodes, Ynodes)
Xg = Xg.flatten()
Yg = Yg.flatten()

velocity_domain_features = make_GPML_velocity_feature(Xg, Yg)
# velocity_domain_features

# Calculate velocities using a delta time interval of 1 My.
delta_time = 1.

# All domain points and associated (magnitude, azimuth, inclination) velocities for the current time.
all_domain_points = []
all_velocities = []

# Partition our velocity domain features into our topological plate polygons at the current 'time'.
plate_partitioner = pygplates.PlatePartitioner(topology_features, rotation_model, reconstruction_time)

# loop through...
for velocity_domain_feature in velocity_domain_features:

    # A velocity domain feature usually has a single geometry but we'll assume it can be any number.
    # Iterate over them all.
    for velocity_domain_geometry in velocity_domain_feature.get_geometries():
        for velocity_domain_point in velocity_domain_geometry.get_points():

            # save the velocity domain point in question
            all_domain_points.append(velocity_domain_point)

            # get the plate
            partitioning_plate = plate_partitioner.partition_point(velocity_domain_point)

            if partitioning_plate:
                # We need the newly assigned plate ID to get the equivalent stage rotation of that tectonic plate.
                partitioning_plate_id = partitioning_plate.get_feature().get_reconstruction_plate_id()

                # Get the stage rotation of partitioning plate from 'reconstruction_time + delta_time' to 'reconstruction_time'.
                equivalent_stage_rotation = rotation_model.get_rotation(reconstruction_time,
                                                                        partitioning_plate_id,
                                                                        reconstruction_time + delta_time)

                # Calculate velocity at the velocity domain point.
                # This is from 'reconstruction_time + delta_time' to 'reconstruction_time' on the partitioning plate.
                velocity_vectors = pygplates.calculate_velocities([velocity_domain_point],
                                                                  equivalent_stage_rotation,
                                                                  delta_time)

                # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
                velocities = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                    [velocity_domain_point], velocity_vectors)
                all_velocities.append(velocities[0])

            else:
                all_velocities.append((0, 0, 0))

# new lists with velocity magnitude and azimuth
pt_vel_mag = []
pt_vel_az = []

# iterate through
for velocity_vector in all_velocities:
    pt_vel_mag.append(velocity_vector[0])
    pt_vel_az.append(velocity_vector[1])

# new lists with our coordinates
pt_lon = []
pt_lat = []
for pt in all_domain_points:
    pt_lon.append(pt.to_lat_lon()[1])
    pt_lat.append(pt.to_lat_lon()[0])


def get_plate_velocities(velocity_domain_features, topology_features,
                         rotation_model, time, delta_time, rep='north_east'):
    """ Function to get velocites via pygplates. Returns the array of tuples "all_velocities",
    which by default will be north-east-down format"""

    # All domain points and associated (magnitude, azimuth, inclination) velocities for the current time.
    all_domain_points = []
    all_velocities = []

    # Partition our velocity domain features into our topological plate polygons at the current 'time'.
    plate_partitioner = pygplates.PlatePartitioner(topology_features, rotation_model, time)

    for velocity_domain_feature in velocity_domain_features:
        # A velocity domain feature usually has a single geometry but we'll assume it can be any number.
        # Iterate over them all.
        for velocity_domain_geometry in velocity_domain_feature.get_geometries():
            for velocity_domain_point in velocity_domain_geometry.get_points():

                all_domain_points.append(velocity_domain_point)

                partitioning_plate = plate_partitioner.partition_point(velocity_domain_point)

                if partitioning_plate:
                    # We need the newly assigned plate ID to get the equivalent stage rotation of that tectonic plate.
                    partitioning_plate_id = partitioning_plate.get_feature().get_reconstruction_plate_id()

                    # Get the stage rotation of partitioning plate from 'time + delta_time' to 'time'.
                    equivalent_stage_rotation = rotation_model.get_rotation(time,
                                                                            partitioning_plate_id,
                                                                            time + delta_time)

                    # Calculate velocity at the velocity domain point.
                    # This is from 'time + delta_time' to 'time' on the partitioning plate.
                    velocity_vectors = pygplates.calculate_velocities(
                        [velocity_domain_point],
                        equivalent_stage_rotation,
                        delta_time)

                    if rep == 'mag_azim':
                        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
                        velocities = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                            [velocity_domain_point],
                            velocity_vectors)
                        all_velocities.append(velocities[0])

                    elif rep == 'north_east':
                        # Convert global 3D velocity vectors to local (north, east, down) tuples (one tuple per point).
                        velocities = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                            [velocity_domain_point],
                            velocity_vectors)
                        all_velocities.append(velocities[0])
                else:
                    all_velocities.append((0, 0, 0))

    return all_velocities


### Animation Loop

outdir = 'tmp'
if not os.path.exists(outdir):
    os.makedirs(outdir)


def reconstruct_and_plot(reconstruction_time):
    if not global_map:
        ax.set_xlim(X_start, X_end)
        ax.set_ylim(Y_start, Y_end)
    # Reconstruct Velocity Arrows to Reconstruction Time

    # Call the function we created above to get the velocities
    all_velocities = get_plate_velocities(velocity_domain_features,
                                          topology_features,
                                          rotation_model,
                                          reconstruction_time,
                                          delta_time,
                                          rep='north_east')

    # prepare for plotting
    pt_vel_n = []
    pt_vel_e = []

    for vel in all_velocities:
        if not hasattr(vel, 'get_y'):
            pt_vel_e.append(vel[1])
            pt_vel_n.append(vel[0])
        else:
            pt_vel_e.append(vel.get_y())
            pt_vel_n.append(vel.get_x())

    u = np.asarray(pt_vel_e).reshape((Ynodes.shape[0], Xnodes.shape[0]))
    v = np.asarray(pt_vel_n).reshape((Ynodes.shape[0], Xnodes.shape[0]))

    # Reconstruct Topology Features to Reconstruction Time

    reconstructed_coastlines = []
    pygplates.reconstruct(coastlines, rotation_model, reconstructed_coastlines, reconstruction_time)
    gdf_coastlines = create_geodataframe_general(reconstructed_coastlines, reconstruction_time)

    # ---
    # use plate tectonic tools to get topologies
    resolved_topologies = ptt.resolve_topologies.resolve_topologies_into_features(
        rotation_model, topology_features, reconstruction_time, transform_segment_deviation_in_radians=0)

    # all the topologies at X time
    topologies, ridge_transforms, ridges, transforms, trenches, trench_left, trench_right, other = resolved_topologies

    # for getting the plates themselves - as polygons.
    gdf_topologies = create_geodataframe_topologies(topologies, reconstruction_time)
    gdf_topologies_plates = gdf_topologies[(gdf_topologies['gpml_type'] == 'gpml:TopologicalClosedPlateBoundary')]

    # for getting components of the topologies as lines - similar to if we exported topologies from GPlates
    gdf_ridge_transforms = create_geodataframe_lines(ridge_transforms, reconstruction_time)

    gdf_trench_left = create_geodataframe_lines(trench_left, reconstruction_time, flip_geom=True)
    gdf_trench_right = create_geodataframe_lines(trench_right, reconstruction_time)
    gdf_trenches = gdf_trench_left.append(gdf_trench_right)

    # Text on map
    ax.text(0.08, 0.95, '%i Ma' % reconstruction_time, transform=ax.transAxes, horizontalalignment='right',
            verticalalignment='bottom', fontsize=26)

    # Plot coastlines by some kind of colour.
    gdf_coastlines.plot(column='PLATEID1', color='grey', ax=ax,
                        transform=coord_handler, alpha=0.75)

    # Plot topologies
    gdf_ridge_transforms.plot(ax=ax, color='firebrick', transform=coord_handler)
    gdf_trenches.plot(ax=ax, color='navy', transform=coord_handler)
    # gdf_reconstructed_hotspots.plot(ax=ax, color='orange', marker='^',  transform=coord_handler)

    # Plot arrows using quiver. This time we use Xnodes and Ynodes, rather than pt_lon, pt_lat
    Q = ax.quiver(Xnodes, Ynodes, u, v, scale=2000, transform=ccrs.PlateCarree(), color='#babcbf')

    # make quiver key.
    plt.quiverkey(Q, 0.11, 0.05, 50, '50 mm/yr', labelpos='W')

    slab.write_subducted_slabs_to_xyz(output_filename, output_data)
    fig.savefig(outdir + '/Slab_Windows_time_%0.2fMa_dip_%0.2fdeg.png' % (reconstruction_time, dip_angle_degrees))

    print('... Saved figure for %s Ma' % reconstruction_time)


### Take saved images and create an animation

plot_times = np.arange(start_time, end_time - time_step, -time_step)
frame_list = []

for i in plot_times:
    ax.clear()
    time_list = np.arange(i + (bin_time * time_step), i - time_step, -time_step)
    end_time = i - time_step
    subducting_slab()
    reconstruct_and_plot(i)
    frame_list.append(outdir + '/Slab_Windows_time_%0.2fMa_dip_%0.2fdeg.png' % (i, dip_angle_degrees))

clip = mpy.ImageSequenceClip(frame_list, fps=fps)

if mp4:
    clip.write_videofile('Slab_Windows_time_%0.2fMa_dip_%0.2fdeg.mp4' % (reconstruction_time, dip_angle_degrees))
else:
    clip.write_gif('Slab_Windows_time_%0.2fMa_dip_%0.2fdeg.gif' % (reconstruction_time, dip_angle_degrees))

if os.path.exists(outdir):
    shutil.rmtree(outdir)
