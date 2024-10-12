```markdown
# Swiss Hiking Trail Planner

This project involves analyzing Swiss hiking trails using various geospatial data processing techniques. It includes operations such as clipping, proximity checks, polygon generation, and filtering based on specific criteria. The primary goal is to identify and visualize hiking routes in relation to parking areas, ensuring that the selected trails meet user-defined conditions.

## Example
[Hiking in Ticino with max Distance 6km- max climbing 400m - next parking lot within 500m](https://s.geo.admin.ch/hw0827f8188j)


## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Key Features](#key-features)
- [Functions](#functions)
- [Contributing](#contributing)
- [License](#license)

## Prerequisites

To run this project, ensure you have the following installed:

- Python 3.x
- Required libraries:
  - `osgeo`
  - `shapely`
  - `pyproj`
  - `osmnx`
  - `simplekml`
  - `geopandas`
  - `pandas`

You can install the required libraries using pip:

```bash
pip install osgeo shapely pyproj osmnx simplekml geopandas pandas
```

## Installation

Clone this repository to your local machine:

```bash
git clone https://github.com/yourusername/swiss-hiking-trail-analysis.git
cd swiss-hiking-trail-analysis
```

## Usage

1. Set your configurations in the script, including paths for input shapefiles and output directories.

2. Run the script to process the hiking trails and generate the desired outputs.
 
```bash
python hiking_multipoly_parking.py
```

### Configurations

- `extent`: Set the bounding box for your analysis.
- `input_shapefile`: Path to the input shapefile containing hiking trails.
- `parking_buffer_meter`: Buffer distance for parking proximity checks.
- `max_length_meter`: Maximum length for valid hiking trails.
- `max_height_meter`: Maximum height for valid hiking trails.
- `output_directory`: Directory where results will be saved.

## Key Features

- Clipping hiking trails based on a specified geographic extent.
- Fetching parking locations from OpenStreetMap (OSM) within a defined bounding box.
- Generating polygons from LineString geometries and merging neighboring polygons.
- Filtering polygons based on proximity to parking areas and specified length/height constraints.
- Exporting results to GeoPackage files and generating KML files for visualization.

## Functions

- `empty_directory(directory_path)`: Removes all files from the specified directory.
- `convert_geopackage_to_epsg2056(input_gpkg, output_gpkg)`: Converts GeoPackage geometries to EPSG:2056.
- `convert_points_to_polygons_and_merge(input_gpkg, output_gpkg, buffer_distance)`: Converts point geometries to polygons and merges them.
- `merge_polygons(polygons)`: Merges neighboring polygons based on proximity.
- `length_to_color(length, max_length_meter)`: Generates a color based on polygon length.
- `transform_coordinates(coords)`: Transforms coordinates from EPSG:2056 to WGS84.

## Contributing

Contributions are welcome! Please create a pull request or open an issue for any feature requests or bugs.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
```

