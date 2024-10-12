from osgeo import ogr, osr
from shapely.geometry import LineString, Polygon, shape
from shapely.ops import polygonize, unary_union
from shapely.strtree import STRtree
import pyproj
import osmnx as ox
import simplekml
import json
import random
import os
import subprocess
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, Polygon

#Configure
extent = ([2698690,1111036],  [2706116,1116231])
input_shapefile= os.path.join("input","dissolved_wanderwege.shp")
output_shapefile = "r_w_dissolved.shp"
parking_shapefile = "r_parking.shp"  # This should be the correct path to your parking data
parking_buffer_meter = 500  # Use 500 meters for proximity check
max_length_meter = 6000  # Define the maximum length in meters
max_height_meter = 400
output_directory = "result"

# Step 0 CLIP SWISS HIKING TRAIL:

min_x, min_y = extent[0]
max_x, max_y = extent[1]

# Create a temporary polygon shapefile for the clipping extent
clip_shapefile = "clip_extent.shp"
if os.path.exists(clip_shapefile):
    os.remove(clip_shapefile)

# Create the polygon WKT for the extent
wkt_polygon = f"POLYGON(({min_x} {min_y}, {max_x} {min_y}, {max_x} {max_y}, {min_x} {max_y}, {min_x} {min_y}))"

# Use ogr2ogr to create a shapefile for the clipping polygon

subprocess.run([
    "ogr2ogr",
    "-f", "ESRI Shapefile",  # Output format
    output_shapefile,        # Output shapefile
    input_shapefile,         # Input shapefile
    "-clipsrc", wkt_polygon,  # Use the dynamically created WKT polygon
    "-a_srs", "EPSG:2056",  # Set the spatial reference
])


print(f"Clipping complete. Result saved in {output_shapefile}")


# Step 0: Get parking from OSM

# OSMnx expects coordinates in EPSG:4326 (latitude/longitude), so we need to reproject from Swiss coordinates (EPSG:2056) to EPSG:4326

proj_2056_to_4326 = pyproj.Transformer.from_crs(2056, 4326, always_xy=True).transform

# Reproject the extent from EPSG:2056 to EPSG:4326
min_lon, min_lat = proj_2056_to_4326(min_x, min_y)
max_lon, max_lat = proj_2056_to_4326(max_x, max_y)

# Create a bounding box in lat/lon
bbox = (min_lat, max_lat, min_lon, max_lon)

# Download OSM data for "amenity=parking" within the bounding box
tags = {"amenity": "parking"}
gdf = ox.geometries_from_bbox(*bbox, tags)

if not gdf.empty:
    print(f"Found {len(gdf)} parking locations.")

    # Save the result to a shapefile
    parking= "parking.gpkg"  # Change to shapefile extension
    gdf.to_file(parking, driver="GPKG")  # Use ESRI Shapefile driver
    print(f"Parking data saved to {parking}")
    def convert_geopackage_to_epsg2056(input_gpkg, output_gpkg):
        # Read the GeoPackage
        gdf = gpd.read_file(input_gpkg)

        # Convert to EPSG:2056
        gdf_2056 = gdf.to_crs(epsg=2056)

        # Write the result back to a new GeoPackage
        gdf_2056.to_file(output_gpkg, driver='GPKG')

        print(f"Conversion complete. Data is now in EPSG:2056 and stored in {output_gpkg}")

    convert_geopackage_to_epsg2056(parking, "parking_2056.gpkg")
    def convert_points_to_polygons_and_merge(input_gpkg, output_gpkg, buffer_distance=20):
        # Read the GeoPackage
        gdf = gpd.read_file(input_gpkg)

        # Separate points and polygons
        points = gdf[gdf.geometry.type == 'Point']
        polygons = gdf[gdf.geometry.type == 'Polygon']

        # Convert points to polygons using buffer
        point_polygons = points.copy()
        point_polygons['geometry'] = points.geometry.buffer(buffer_distance)

        # Merge the new polygons with existing polygons
        all_polygons = gpd.GeoDataFrame(pd.concat([polygons, point_polygons], ignore_index=True))

        # Ensure the CRS is EPSG:2056
        all_polygons.crs = "EPSG:2056"

        # Write the result back to a new GeoPackage
        all_polygons.to_file(output_gpkg, driver='ESRI Shapefile')

        print(f"Conversion complete. All data is now stored as polygons in {output_gpkg}")

    convert_points_to_polygons_and_merge("parking_2056.gpkg", parking_shapefile)
else:
    print("No parking locations found in the defined extent.")



# Step 1: Load the LineString shapefile

driver = ogr.GetDriverByName("ESRI Shapefile")
input_ds = driver.Open(output_shapefile, 0)  # 0 means read-only
input_layer = input_ds.GetLayer()

# Create output directory if it doesn't exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Extract the MultiLineString geometry
line_strings = []
for feature in input_layer:
    geom = feature.GetGeometryRef()

    # Check if it's a MultiLineString
    if geom.GetGeometryName() == "MULTILINESTRING":
        # Loop through each LineString inside the MultiLineString
        for i in range(geom.GetGeometryCount()):
            linestring = geom.GetGeometryRef(i)
            coords = linestring.GetPoints()  # Extract points for each LineString
            line_strings.append(LineString(coords))

# Check if we have extracted any valid LineStrings
if not line_strings:
    print("No valid LineStrings were found within the MultiLineString.")
else:
    print(f"Found {len(line_strings)} LineString(s).")

# Polygonize the LineStrings to get base polygons
lines = unary_union(line_strings)  # Merge all lines into one multiline
polygons = list(polygonize(lines))  # Create polygons from the closed loops

# Output the number of generated polygons
print(f"Generated {len(polygons)} polygons.")

# Function to merge polygons based on neighbors
def merge_polygons(polygons):
    # Ensure that input polygons are valid shapely Polygon objects
    if not all(isinstance(p, Polygon) for p in polygons):
        raise ValueError("All elements in the input list must be shapely Polygon objects.")

    merged_polygons = []
    combination_count = 0

    # Create a spatial index with the original polygons
    spatial_index = STRtree(polygons)

    # Iterate over each polygon
    for polygon in polygons:
        # Query the spatial index for nearby polygons (geometry objects)
        nearby_indices = spatial_index.query(polygon)

        # Use the indices to get the actual polygons
        nearby_polygons = [polygons[i] for i in nearby_indices]

        # Find neighbors that touch or intersect with the current polygon
        neighbors = [p for p in nearby_polygons if p is not polygon and (p.touches(polygon) or p.intersects(polygon))]

        # If there are neighbors, try to merge them
        if neighbors:
            possible_merge = unary_union([polygon] + neighbors)
            if isinstance(possible_merge, Polygon) and possible_merge.is_valid:
                merged_polygons.append(possible_merge)
                combination_count += 1
                print(f"Merged {combination_count} polygon(s)")

    return merged_polygons

# Get all merged polygons
all_polygons = polygons + merge_polygons(polygons)
print("Get all merged polygons DONE")

# Write each polygon to its own GeoPackage file with EPSG:2056 CRS
for idx, polygon in enumerate(all_polygons):
    output_gpkg_path = f"{output_directory}/polygon_{idx}.gpkg"
    output_driver = ogr.GetDriverByName("GPKG")
    output_ds = output_driver.CreateDataSource(output_gpkg_path)

    # Create the spatial reference for EPSG:2056
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(2056)  # EPSG:2056

    output_layer = output_ds.CreateLayer(f"polygon_{idx}", geom_type=ogr.wkbPolygon, srs=srs)

    # Create the feature with the merged polygon
    feature_def = output_layer.GetLayerDefn()
    new_feature = ogr.Feature(feature_def)

    # Check geometry validity
    if polygon.is_valid:
        new_feature.SetGeometry(ogr.CreateGeometryFromWkb(polygon.wkb))
        output_layer.CreateFeature(new_feature)
    else:
        print(f"Invalid polygon geometry for polygon_{idx}. Skipping saving.")

    # Cleanup
    new_feature = None
    output_ds = None

input_ds = None
print("Step 1 complete: Individual polygons saved as separate GeoPackage files with EPSG:2056 CRS.")

# Step 2: Process the generated GeoPackage files and filter polygons by proximity to parking areas
# Load the parking areas shapefile
parking_driver = ogr.GetDriverByName("ESRI Shapefile")
parking_ds = parking_driver.Open(parking_shapefile, 0)  # 0 means read-only
parking_layer = parking_ds.GetLayer()

# Create a list to store all parking geometries
parking_geometries = []
for feature in parking_layer:
    geom = feature.GetGeometryRef()
    # Convert to GeoJSON string and then to Shapely shape
    geojson = geom.ExportToJson()  # Export to GeoJSON
    parking_geometries.append(shape(json.loads(geojson)))  # Convert to Shapely shape

# Create a union of all parking geometries for proximity checking
all_parking = unary_union(parking_geometries)

# Iterate through GeoPackage files in the output directory
for filename in os.listdir(output_directory):
    if filename.endswith(".gpkg"):
        gpkg_path = os.path.join(output_directory, filename)

        # Open the GeoPackage
        gpkg_driver = ogr.GetDriverByName("GPKG")
        gpkg_ds = gpkg_driver.Open(gpkg_path, 0)  # 0 means read-only
        gpkg_layer = gpkg_ds.GetLayer()

        # Create an output list for valid polygons
        valid_polygons = []

        # Check each polygon in the GeoPackage
        for feature in gpkg_layer:
            geom = feature.GetGeometryRef()
            # Convert OGR geometry to GeoJSON string
            geojson = geom.ExportToJson()  # Export to GeoJSON
            shapely_geom = shape(json.loads(geojson))  # Convert to Shapely shape

            # Check if the polygon is within 500m of any parking area
            if all_parking.distance(shapely_geom) <= parking_buffer_meter:  # 500m proximity check
                valid_polygons.append(shapely_geom)

        # Check if valid polygons exist, if not delete the GeoPackage file
        if not valid_polygons:
            print(f"Deleting {filename} as it is not withing {parking_buffer_meter} meters...")
            os.remove(gpkg_path)  # Delete the GeoPackage file

        # Cleanup
        gpkg_ds = None  # Cleanup after processing the GeoPackage

print("Step 2 complete: Filtered polygons saved to new GeoPackage files based on parking proximity.")

# Step 3: Remove GeoPackage files containing polygons longer than 6000 meters


for filename in os.listdir(output_directory):
    if filename.endswith(".gpkg"):
        gpkg_path = os.path.join(output_directory, filename)

        # Open the GeoPackage
        gpkg_driver = ogr.GetDriverByName("GPKG")
        gpkg_ds = gpkg_driver.Open(gpkg_path, 0)  # 0 means read-only
        gpkg_layer = gpkg_ds.GetLayer()

        # Check each polygon in the GeoPackage
        for feature in gpkg_layer:
            geom = feature.GetGeometryRef()
            # Convert OGR geometry to GeoJSON string
            geojson = geom.ExportToJson()  # Export to GeoJSON
            shapely_geom = shape(json.loads(geojson))  # Convert to Shapely shape

            # Check the length of the polygon
            if shapely_geom.length > max_length_meter or shapely_geom.length < 1000:
                print(f"Deleting {filename} as it contains polygons longer than {max_length_meter} meters.")
                os.remove(gpkg_path)  # Delete the GeoPackage file
                break  # Exit loop once the file is marked for deletion

        # Cleanup
        gpkg_ds = None  # Cleanup after processing the GeoPackage

print("Step 3 complete: GeoPackage files deleted if they contained polygons longer than 6000 meters.")
# Step 4: Sum the heights along the polygons and delete files with total height greater than 300 meters
for filename in os.listdir(output_directory):
    if filename.endswith(".gpkg"):
        gpkg_path = os.path.join(output_directory, filename)

        # Open the GeoPackage
        gpkg_driver = ogr.GetDriverByName("GPKG")
        gpkg_ds = gpkg_driver.Open(gpkg_path, 0)  # 0 means read-only
        gpkg_layer = gpkg_ds.GetLayer()

        # Check each polygon in the GeoPackage
        total_height = 0  # Initialize total height sum
        for feature in gpkg_layer:
            geom = feature.GetGeometryRef()
            # Convert OGR geometry to GeoJSON string
            geojson = geom.ExportToJson()  # Export to GeoJSON
            shapely_geom = shape(json.loads(geojson))  # Convert to Shapely shape

            # Calculate the sum of differences in heights along the polygon
            heights = [point[2] for point in shapely_geom.exterior.coords]  # Extract heights (3rd value in the vertex)

            # Calculate the differences and sum them
            height_differences = [abs(heights[i] - heights[i - 1]) for i in range(1, len(heights))]
            total_height += sum(height_differences)

        # Check if the total height exceeds the limit
        if total_height > max_height_meter:
            print(f"Deleting {filename} as total height {total_height} meters exceeds {max_height_meter} meters.")
            os.remove(gpkg_path)  # Delete the GeoPackage file

        # Cleanup
        gpkg_ds = None  # Cleanup after processing the GeoPackage

print("Step 4 complete: GeoPackage files deleted if they contained polygons with total height greater than 300 meters.")
# Function to determine the color based on length
""" def length_to_color(length, max_length_meter):
    normalized_length = min(max(length, 1000), max_length_meter)  # Clamp length between 1000 and max_length_meter

    if normalized_length <= (1000 + max_length_meter) / 2:
        # Transition from Green (at 1000) to Yellow (at midpoint)
        r = int(255 * ((normalized_length - 1000) / (max_length_meter - 1000)))  # Red increases
        g = 255  # Green stays at maximum
    else:
        # Transition from Yellow (midpoint) to Red (at max_length_meter)
        r = 255  # Red stays at maximum
        g = int(255 * (1 - (normalized_length - (1000 + max_length_meter) / 2) / (max_length_meter - (1000 + max_length_meter) / 2)))  # Green decreases

    b = 0  # Blue remains constant
    alpha = 255  # Fully opaque

    # Format: aabbggrr
    return f'{alpha:02x}{b:02x}{g:02x}{r:02x}'  # Return hex color """

def length_to_color(length, max_length_meter):
    # Generate random values for red, green, and blue channels
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)

    alpha = 255  # Fully opaque

    # Format: aabbggrr
    return f'{alpha:02x}{b:02x}{g:02x}{r:02x}'  # Return hex color with random RGB values

# Function to convert coordinates from EPSG:2056 to WGS84
def transform_coordinates(coords):
    source = osr.SpatialReference()
    source.ImportFromEPSG(2056)  # EPSG:2056
    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)  # EPSG:4326

    transform = osr.CoordinateTransformation(source, target)
    return [(transform.TransformPoint(x, y)[:2]) for x, y, _ in coords]

# Step 5: Convert valid GeoPackages to KML using simplekml
# List to hold placemark data
placemark_data = []

# Step 5: Convert valid GeoPackages to one KML file using simplekml
for filename in os.listdir(output_directory):
    if filename.endswith(".gpkg"):
        gpkg_path = os.path.join(output_directory, filename)

        # Open the GeoPackage
        gpkg_driver = ogr.GetDriverByName("GPKG")
        gpkg_ds = gpkg_driver.Open(gpkg_path, 0)  # 0 means read-only
        gpkg_layer = gpkg_ds.GetLayer()

        # Process each polygon
        for feature in gpkg_layer:
            geom = feature.GetGeometryRef()
            geojson = geom.ExportToJson()  # Export to GeoJSON
            shapely_geom = shape(json.loads(geojson))  # Convert to Shapely shape

            # Calculate total height and length
            heights = [point[2] for point in shapely_geom.exterior.coords]
            height_differences = [abs(heights[i] - heights[i - 1]) for i in range(1, len(heights))]
            total_height = sum(height_differences)
            length = shapely_geom.length

            # Round heights and length to integers
            total_height = round(total_height)
            length = round(length)

            # Update the max values if necessary
            if total_height > max_height_meter:
                max_height_meter = total_height
            if length > max_length_meter:
                max_length_meter = length

            # Convert coordinates to WGS84
            coords_wgs84 = transform_coordinates(shapely_geom.exterior.coords)

            # Determine the color based on length
            color = length_to_color(length, max_length_meter)

            # Determine the thickness based on total height (scaled)
            thickness = round(max(min(total_height / max_height_meter * 10, 10), 1))

            # Convert length to kilometers and round to one decimal place
            length_km = round(length / 1000, 1)

            # Store placemark data
            placemark_data.append({
                'name': f"Länge:{length_km}km Höhe:{total_height}m",
                'coords': [(lon, lat, 0) for lat, lon in coords_wgs84],
                'color': color,
                'thickness': thickness
            })

            # Save each polygon as a separate KML file with {length}_{total_height}.kml
            individual_kml = simplekml.Kml()
            individual_placemark = individual_kml.newlinestring(name=f"Länge:{length_km}km_Höhe:{total_height}m")
            individual_placemark.coords = [(lon, lat, 0) for lat, lon in coords_wgs84]
            individual_placemark.style.linestyle.color = color  # Set color
            individual_placemark.style.linestyle.width = thickness  # Set thickness

            # Add extended data
            individual_placemark.extendeddata.newdata(name="type", value="linepolygon")

            # Save the individual KML file
            individual_kml_filename = f"{length}_{total_height}.kml"
            individual_kml_path = os.path.join(output_directory, individual_kml_filename)
            individual_kml.save(individual_kml_path)

        # Cleanup
        gpkg_ds = None  # Cleanup after processing the GeoPackage

# Sort the placemark data by thickness (from largest to smallest)
placemark_data.sort(key=lambda x: x['thickness'], reverse=True)

# Create a KML object for the combined file
kml = simplekml.Kml()

# Set the KML document name to be {max_height_meter}_{max_length_meter}
kml.document.name = f"{max_height_meter}_{max_length_meter}"

# Add sorted placemarks to the KML
for placemark_info in placemark_data:
    placemark = kml.newlinestring(name=placemark_info['name'])
    placemark.coords = placemark_info['coords']
    placemark.style.linestyle.color = placemark_info['color']  # Set color
    placemark.style.linestyle.width = placemark_info['thickness']  # Set thickness

    # Add extended data
    placemark.extendeddata.newdata(name="type", value="linepolygon")

# Save the combined KML file using the max height and max length in the filename
kml_filename = f"{max_height_meter}_{max_length_meter}.kml"
kml_path = os.path.join(output_directory, kml_filename)
kml.save(kml_path)

print(f"Step 5 complete: All GeoPackage files converted to {kml_filename} with the specified structure and sorted placemarks.")

files_to_delete = ["parking__merged_2056.gpkg", "parking_2056.gpkg", "parking.gpkg"]
# Loop through the list and delete each file
for file in files_to_delete:
    try:
        os.remove(file)
        print(f"Deleted: {file}")
    except:
        print(f"File not found: {file}")