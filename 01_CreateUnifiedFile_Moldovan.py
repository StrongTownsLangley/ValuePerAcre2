import json
import csv
import math
from shapely.geometry import shape, Point, Polygon, mapping
import shapely.ops
from tqdm import tqdm  # Import tqdm for progress bars
from rtree import index  # Import rtree for spatial indexing

def load_json_file(file_path):
    """Load a JSON file with UTF-8 encoding to avoid character decoding issues."""
    with open(file_path, 'r', encoding='utf-8') as file:
        return json.load(file)

def load_csv_file(file_path):
    """Load a CSV file containing parcel tax information."""
    tax_data = {}
    with open(file_path, 'r', encoding='utf-8') as file:
        reader = csv.reader(file)
        header = next(reader)  # Skip header row
        
        # Identify column indices
        parcel_id_idx = 0  # ParcelMaster.parcelnumber is the first column and matches with MSTPARCEL
        taxable_value_idx = 1  # Parcels.current_taxable is the second column
        
        print(f"CSV Header: {header}")
        
        for row in reader:
            try:
                parcel_id = row[parcel_id_idx].strip()
                taxable_value = float(row[taxable_value_idx]) if row[taxable_value_idx].strip() else 0
                tax_data[parcel_id] = taxable_value
            except (IndexError, ValueError) as e:
                print(f"Error processing row: {row}, Error: {e}")
    
    print(f"Loaded {len(tax_data)} tax records from CSV")
    return tax_data

def extract_point_from_geometry(geometry):
    """Safely extract coordinates from a Point geometry."""
    try:
        if geometry['type'] == 'Point':
            return geometry['coordinates']
        else:
            # If not a point, try to get the centroid
            geom_shape = shape(geometry)
            centroid = geom_shape.centroid
            return [centroid.x, centroid.y]
    except Exception as e:
        print(f"Error extracting point from geometry: {e}")
        return None

def fill_polygon_holes(geometry):
    """
    Fill the holes in a polygon or multipolygon geometry.
    This will create a solid polygon without interior rings.
    """
    try:
        # Convert the GeoJSON geometry to a shapely geometry
        geom = shape(geometry)
        
        if geom.geom_type == 'Polygon':
            # Create a new polygon without holes (just the exterior ring)
            exterior = geom.exterior
            filled_poly = Polygon(exterior.coords)
            return mapping(filled_poly)
            
        elif geom.geom_type == 'MultiPolygon':
            # Handle each polygon in the multipolygon
            filled_polys = []
            for poly in geom.geoms:
                exterior = poly.exterior
                filled_poly = Polygon(exterior.coords)
                filled_polys.append(filled_poly)
            
            # Create a new multipolygon without holes
            filled_multipoly = shapely.ops.unary_union(filled_polys)
            return mapping(filled_multipoly)
        
        # If not a polygon/multipolygon, return original
        return geometry
    except Exception as e:
        print(f"Error filling polygon holes: {e}")
        return geometry

def identify_complex_polygons(parcel_shapes):
    """
    Identify complex developments where multiple buildings share a common area.
    These typically have complex shapes or large parcels with holes.
    """
    complex_polygons = {}
    
    for parcel_id, parcel_data in parcel_shapes.items():
        properties = parcel_data.get("properties", {})
        geometry = parcel_data.get("geometry", {})
        
        # Check for common keywords in property fields that might indicate complex developments
        is_complex = False
        for field, value in properties.items():
            if isinstance(value, str) and any(keyword in value.upper() for keyword in ["COMPLEX", "CMPLX", "COMMON", "CONDO"]):
                is_complex = True
                break
        
        # Check if it's a polygon with holes (interior rings)
        has_holes = False
        try:
            if geometry.get('type') == 'Polygon':
                # A polygon with holes has more than one coordinate array
                has_holes = len(geometry.get('coordinates', [])) > 1
            elif geometry.get('type') == 'MultiPolygon':
                # For multipolygons, check each polygon for holes
                for poly_coords in geometry.get('coordinates', []):
                    if len(poly_coords) > 1:
                        has_holes = True
                        break
        except Exception:
            pass
        
        # Check if it's a multipolygon with many parts (potential apartment complex)
        is_multipart = False
        try:
            if geometry.get('type') == 'MultiPolygon':
                is_multipart = len(geometry.get('coordinates', [])) > 3  # Arbitrary threshold
        except Exception:
            pass
        
        # Check for large area parcels which might be complexes
        is_large_parcel = False
        try:
            area = shape(geometry).area
            # Set an arbitrary threshold for "large" parcels
            is_large_parcel = area > 10000  # Adjust based on your data's scale
        except Exception:
            pass
            
        # Mark as complex if any criteria match
        if is_complex or has_holes or is_multipart or is_large_parcel:
            print(f"Found complex development: {parcel_id}")
            
            # For complex polygons, fill the holes to ensure we catch all units
            if has_holes:
                # Create a copy of the feature with filled holes
                filled_geometry = fill_polygon_holes(geometry)
                
                # Create a new shape from the filled geometry
                filled_shape = shape(filled_geometry)
                
                # Store both the original and filled versions
                complex_polygons[parcel_id] = {
                    "original": parcel_data,
                    "filled": {
                        "geometry": filled_geometry,
                        "shape": filled_shape,
                        "properties": properties
                    }
                }
            else:
                # Just store a reference to the original data
                complex_polygons[parcel_id] = {
                    "original": parcel_data,
                    "filled": parcel_data  # No filling needed
                }
    
    return complex_polygons

def find_parcels_within_complex(complex_id, complex_data, all_parcels, spatial_idx, idx_to_parcel_id):
    """
    Find all parcels that are located within a complex development.
    Uses the filled (hole-free) version of complex polygons to ensure all units are captured.
    """
    # Use the filled (no holes) shape for detection
    complex_shape = complex_data["filled"]["shape"]
    contained_parcels = []
    
    # Get complex bounds for spatial query
    minx, miny, maxx, maxy = complex_shape.bounds
    
    # Query the spatial index
    potential_parcels = list(spatial_idx.intersection((minx, miny, maxx, maxy)))
    
    # Get property info for matching
    complex_props = complex_data["filled"]["properties"]
    
    # Track everything we find for debugging
    all_found = []
    
    for idx in potential_parcels:
        parcel_id = idx_to_parcel_id.get(idx)
        
        if parcel_id and parcel_id != complex_id and parcel_id in all_parcels:
            parcel_data = all_parcels[parcel_id]
            parcel_shape = parcel_data["shape"]
            
            # Different criteria for containment
            criteria = {
                "centroid_inside": False,
                "mostly_overlapping": False
            }
            
            # Check if parcel centroid is within the complex
            parcel_centroid = parcel_shape.centroid
            criteria["centroid_inside"] = complex_shape.contains(parcel_centroid)
            
            # If not, check if polygons overlap substantially
            if not criteria["centroid_inside"]:
                try:
                    if complex_shape.intersects(parcel_shape):
                        overlap_area = complex_shape.intersection(parcel_shape).area
                        parcel_area = parcel_shape.area
                        if parcel_area > 0 and overlap_area / parcel_area > 0.5:  # More than 50% overlap
                            criteria["mostly_overlapping"] = True
                except Exception:
                    pass
            
            # Determine if this parcel should be included
            should_include = criteria["centroid_inside"] or criteria["mostly_overlapping"]
            
            # Record this for analysis
            all_found.append({
                "parcel_id": parcel_id,
                "criteria": criteria,
                "included": should_include
            })
            
            # If it meets our criteria, include it
            if should_include:
                contained_parcels.append(parcel_id)
    
    # Debug: analyze what we found
    if len(all_found) > 0:
        print(f"Complex {complex_id}: Found {len(contained_parcels)} contained parcels out of {len(all_found)} candidates")
    
    return contained_parcels

def main():
    # Load the data files
    print("Loading data files...")
    tax_data = load_csv_file("moldovan_taxparcelinformation.csv")
    parcels_geojson = load_json_file("moldovan_parcels.json")
    
    print(f"Loaded tax data for {len(tax_data)} parcels")
    print(f"Loaded parcel geometries: {len(parcels_geojson['features'])} features")
    
    # Check for sample matches to verify data alignment
    sample_count = 0
    for feature in parcels_geojson["features"][:10]:  # Check first 10 for samples
        mst_parcel = feature["properties"].get("MSTPARCEL")
        if mst_parcel and mst_parcel in tax_data:
            print(f"Sample match found: MSTPARCEL={mst_parcel}, Taxable Value={tax_data[mst_parcel]}")
            sample_count += 1
    
    print(f"Found {sample_count} sample matches between GeoJSON and tax data")
    
    # Process polygon parcels
    parcels_with_shapes = {}
    unified_parcels = {}
    
    # Process parcels from the parcels file
    print("Processing polygon parcels...")
    for feature in tqdm(parcels_geojson["features"], desc="Processing parcels"):
        # Use MSTPARCEL field to match with the CSV data
        parcel_id = feature["properties"].get("PIN") or feature["properties"].get("PARCEL")
        mst_parcel_id = feature["properties"].get("MSTPARCEL")  # This is the field that matches the CSV
        
        if parcel_id is None or parcel_id == "":
            continue  # Skip parcels without ID
        
        try:
            # Create a shapely geometry for spatial operations
            parcel_shape = shape(feature["geometry"])
            
            if parcel_shape.is_empty:
                continue  # Skip empty geometries
                
            # Try to fix any invalid geometries
            if not parcel_shape.is_valid:
                try:
                    parcel_shape = parcel_shape.buffer(0)  # Buffer by 0 often fixes invalid geometries
                except Exception:
                    # If buffer fails, try to keep going anyway
                    pass
            
            # Get the taxable value from our CSV data using MSTPARCEL as the key
            taxable_value = tax_data.get(mst_parcel_id, 0)
            
            # Store the parcel info along with its shape
            parcels_with_shapes[parcel_id] = {
                "geometry": feature["geometry"],
                "shape": parcel_shape,
                "properties": feature["properties"],
                "combined_units": 0,
                "taxable_value": taxable_value,
                "combined_parcels": []
            }
            
            # Initialize in unified parcels too
            unified_parcels[parcel_id] = {
                "geometry": feature["geometry"],
                "properties": feature["properties"],
                "combined_units": 0,
                "combined_parcels": [],
                "taxable_value": taxable_value,
                "parent_parcel": None  # Will be set if this parcel is inside another
            }
        except Exception as e:
            print(f"Error processing parcel {parcel_id}: {e}")
    
    print(f"Processed {len(parcels_with_shapes)} valid parcels")
    
    # Create spatial index for efficient overlap checking
    print("Creating spatial index...")
    spatial_idx = index.Index()
    idx_to_parcel_id = {}
    
    # Add all parcels to the spatial index
    for idx, (parcel_id, parcel_data) in enumerate(parcels_with_shapes.items()):
        parcel_shape = parcel_data["shape"]
        minx, miny, maxx, maxy = parcel_shape.bounds
        spatial_idx.insert(idx, (minx, miny, maxx, maxy))
        # Store the mapping from index to parcel_id
        idx_to_parcel_id[idx] = parcel_id
    
    # Step 1: Identify complex multi-building developments
    print("Identifying complex developments...")
    complex_developments = identify_complex_polygons(parcels_with_shapes)
    print(f"Found {len(complex_developments)} complex developments")
    
    # Step 2: For each complex, find all parcels within it
    contained_parcels = {}
    container_parcels = {}
    
    for complex_id, complex_data in tqdm(complex_developments.items(), desc="Processing complex developments"):
        contained = find_parcels_within_complex(
            complex_id, complex_data, parcels_with_shapes, spatial_idx, idx_to_parcel_id
        )
        
        # Record the relationship
        if contained:
            container_parcels[complex_id] = contained
            for child_id in contained:
                contained_parcels[child_id] = complex_id
                
                # Update parent-child relationship in unified parcels
                if child_id in unified_parcels and complex_id in unified_parcels:
                    unified_parcels[child_id]["parent_parcel"] = complex_id
    
    # Step 3: Regular overlap detection for remaining parcels
    print("Processing remaining overlaps...")
    for parcel_id, parcel_data in tqdm(parcels_with_shapes.items(), desc="Finding additional overlaps"):
        # Skip parcels already identified as part of a complex
        if parcel_id in contained_parcels or parcel_id in container_parcels:
            continue
            
        parcel_shape = parcel_data["shape"]
        parcel_centroid = parcel_shape.centroid
        
        # Get bounds for spatial index query
        minx, miny, maxx, maxy = parcel_shape.bounds
        
        # Query the spatial index for potential containing parcels
        for idx in spatial_idx.intersection((minx, miny, maxx, maxy)):
            other_id = idx_to_parcel_id.get(idx)
            
            if (other_id and other_id != parcel_id and 
                other_id not in contained_parcels and 
                other_id in parcels_with_shapes):
                
                other_data = parcels_with_shapes[other_id]
                other_shape = other_data["shape"]
                
                # Skip parcels that are smaller
                if other_shape.area <= parcel_shape.area:
                    continue
                    
                # Check if this parcel's centroid is inside the other parcel
                if other_shape.contains(parcel_centroid):
                    # Record the relationship
                    contained_parcels[parcel_id] = other_id
                    
                    if other_id not in container_parcels:
                        container_parcels[other_id] = []
                        
                    container_parcels[other_id].append(parcel_id)
                    
                    # Update parent-child relationship in unified parcels
                    if parcel_id in unified_parcels and other_id in unified_parcels:
                        unified_parcels[parcel_id]["parent_parcel"] = other_id
                        
                    break  # Found a container, no need to check others
    
    print(f"Found {len(contained_parcels)} parcels contained within others")
    
    # For any container parcel, ensure it has the combined values from all its contained parcels
    for container_id, contained_ids in container_parcels.items():
        if container_id in unified_parcels:
            container_data = unified_parcels[container_id]
            
            # Initialize if not already done
            if "combined_units" not in container_data or container_data["combined_units"] == 0:
                container_data["combined_units"] = len(contained_ids)
                
            if "combined_parcels" not in container_data or not container_data["combined_parcels"]:
                container_data["combined_parcels"] = contained_ids
                
            # Sum up taxable values from contained parcels
            for child_id in contained_ids:
                if child_id in unified_parcels:
                    container_data["taxable_value"] += unified_parcels[child_id]["taxable_value"]
    
    # Create the output GeoJSON - exclude parcels that are contained within others
    output_features = []
    excluded_count = 0
    
    print("Creating final output file...")
    
    # Create a set of all parcels that should be excluded (all combined parcels)
    all_excluded_parcels = set()
    for container_id, contained_ids in container_parcels.items():
        if container_id in unified_parcels:
            all_excluded_parcels.update(contained_ids)
    
    # Also exclude parcels that are directly marked as contained
    for child_id, parent_id in contained_parcels.items():
        all_excluded_parcels.add(child_id)

    print(f"Total parcels to be excluded from output: {len(all_excluded_parcels)}")
    
    # Create features for output, excluding contained parcels
    for parcel_id, parcel_data in tqdm(unified_parcels.items(), desc="Creating output features"):
        # Skip parcels that are contained within others or listed as combined into another parcel
        if parcel_id in all_excluded_parcels:
            excluded_count += 1
            continue
            
        # Skip parcels with no taxable value and no combined units 
        if parcel_data["taxable_value"] == 0 and parcel_data["combined_units"] == 0:
            continue
        
        # Special handling for complex parcels
        if parcel_id in container_parcels and container_parcels[parcel_id]:
            # Ensure contained units are properly recorded
            parcel_data["combined_units"] = max(parcel_data["combined_units"], len(container_parcels[parcel_id]))
            
            # Make sure all contained parcels are in the combined_parcels list
            for contained_id in container_parcels[parcel_id]:
                if contained_id not in parcel_data["combined_parcels"]:
                    parcel_data["combined_parcels"].append(contained_id)
                    
            # Update the complex geometry to fill in the holes
            if parcel_id in complex_developments:
                # Use the filled geometry without holes
                filled_geometry = complex_developments[parcel_id]["filled"]["geometry"]
                parcel_data["geometry"] = filled_geometry
        
        # Update the properties to include the combined unit count and taxable value
        parcel_data["properties"]["TaxableValue"] = parcel_data["taxable_value"]
        parcel_data["properties"]["CombinedUnits"] = parcel_data["combined_units"]
        parcel_data["properties"]["CombinedParcels"] = ",".join(parcel_data["combined_parcels"]) if parcel_data["combined_parcels"] else ""
        
        # Create a feature for this parcel
        feature = {
            "type": "Feature",
            "id": parcel_id,
            "geometry": parcel_data["geometry"],
            "properties": parcel_data["properties"]
        }
        
        output_features.append(feature)
    
    # Create the final GeoJSON
    output_geojson = {
        "type": "FeatureCollection",
        "crs": parcels_geojson.get("crs", {
            "type": "name",
            "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}
        }),  # Use the original CRS or default to WGS 84
        "features": output_features
    }
    
    # Write to file
    with open("moldovan_unified_parcels.geojson", "w", encoding="utf-8") as outfile:
        json.dump(output_geojson, outfile)
    
    print(f"Process complete!")
    print(f"Total parcels processed: {len(parcels_with_shapes)}")
    print(f"Parcels in output file: {len(output_features)}")
    print(f"Parcels excluded (contained within others): {excluded_count}")
    print(f"Output written to moldovan_unified_parcels.geojson")

if __name__ == "__main__":
    main()