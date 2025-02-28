import json
import sys
from pyproj import Transformer

def convert_coordinates(input_filename, output_filename=None):
    """
    Convert GeoJSON coordinates from EPSG:26910 to EPSG:4326
    
    Args:
        input_filename (str): Path to input GeoJSON file
        output_filename (str, optional): Path to output GeoJSON file. If None, will use input_filename with "_4326" suffix
    """
    if output_filename is None:
        output_filename = input_filename.replace('.geojson', '_4326.geojson')
    
    # Create transformer from EPSG:26910 to EPSG:4326 (WGS84)
    transformer = Transformer.from_crs("EPSG:26910", "EPSG:4326", always_xy=True)
    
    # Read input GeoJSON
    with open(input_filename, 'r') as f:
        geojson_data = json.load(f)
    
    # Update CRS property
    if 'crs' in geojson_data:
        geojson_data['crs'] = {
            "type": "name",
            "properties": {
                "name": "EPSG:4326"
            }
        }
    
    # Transform coordinates in all features
    for feature in geojson_data['features']:
        geometry = feature.get('geometry', {})
        if not geometry:
            continue
            
        geom_type = geometry.get('type')
        
        if geom_type == 'Point':
            x, y = geometry['coordinates']
            lon, lat = transformer.transform(x, y)
            geometry['coordinates'] = [lon, lat]
        
        elif geom_type == 'LineString':
            new_coords = []
            for coord in geometry['coordinates']:
                x, y = coord
                lon, lat = transformer.transform(x, y)
                new_coords.append([lon, lat])
            geometry['coordinates'] = new_coords
        
        elif geom_type == 'Polygon':
            new_rings = []
            for ring in geometry['coordinates']:
                new_coords = []
                for coord in ring:
                    x, y = coord
                    lon, lat = transformer.transform(x, y)
                    new_coords.append([lon, lat])
                new_rings.append(new_coords)
            geometry['coordinates'] = new_rings
            
        elif geom_type == 'MultiPoint':
            new_points = []
            for point in geometry['coordinates']:
                x, y = point
                lon, lat = transformer.transform(x, y)
                new_points.append([lon, lat])
            geometry['coordinates'] = new_points
            
        elif geom_type == 'MultiLineString':
            new_lines = []
            for line in geometry['coordinates']:
                new_coords = []
                for coord in line:
                    x, y = coord
                    lon, lat = transformer.transform(x, y)
                    new_coords.append([lon, lat])
                new_lines.append(new_coords)
            geometry['coordinates'] = new_lines
            
        elif geom_type == 'MultiPolygon':
            new_polygons = []
            for polygon in geometry['coordinates']:
                new_rings = []
                for ring in polygon:
                    new_coords = []
                    for coord in ring:
                        x, y = coord
                        lon, lat = transformer.transform(x, y)
                        new_coords.append([lon, lat])
                    new_rings.append(new_coords)
                new_polygons.append(new_rings)
            geometry['coordinates'] = new_polygons
    
    # Write transformed GeoJSON to output file in condensed format
    with open(output_filename, 'w') as f:
        json.dump(geojson_data, f, separators=(',', ':'))
    
    print(f"Converted {input_filename} from EPSG:26910 to EPSG:4326")
    print(f"Output saved to {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python geojson_converter.py input_file.geojson [output_file.geojson]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    convert_coordinates(input_file, output_file)