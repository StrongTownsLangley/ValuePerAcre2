import json
import math
from shapely.geometry import shape, Point, Polygon, mapping
import shapely.ops
from tqdm import tqdm  # Import tqdm for progress bars
from rtree import index  # Import rtree for spatial indexing
from jsonfileloader import JsonFileLoader
from taxablevaluecalculator import TaxableValueCalculator
from totalassessedvaluecalculator import TotalAssessedValueCalculator
from geometrypointextractor import GeometryPointExtractor
from polygonholefiller import PolygonHoleFiller

def identify_complex_polygons(parcel_shapes):
    """
    Identify complex developments where multiple buildings share a common area.
    These typically have "CMPLX" or similar in their unit field, or are large parcels
    with holes or complex shapes.
    """
    complex_polygons = {}
    
    for parcel_id, parcel_data in parcel_shapes.items():
        properties = parcel_data.get("properties", {})
        unit = properties.get("UNIT")
        geometry = parcel_data.get("feature", {}).get("geometry", {})
        
        # Check if it's marked as a complex
        is_complex = unit and (unit == "CMPLX" or "COMPLEX" in unit or "COMMON" in unit)
        
        # Or check if it's a polygon with holes (interior rings)
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
        
        # Or check if it's a multipolygon with many parts (potential apartment complex)
        is_multipart = False
        try:
            if geometry.get('type') == 'MultiPolygon':
                is_multipart = len(geometry.get('coordinates', [])) > 3  # Arbitrary threshold
        except Exception:
            pass
            
        # Also check for "common" in the street name or other fields
        common_in_fields = False
        for field, value in properties.items():
            if isinstance(value, str) and "COMMON" in value.upper():
                common_in_fields = True
                break
                
        # Mark as complex if any criteria match
        if is_complex or has_holes or is_multipart or common_in_fields:
            print(f"Found complex development: {parcel_id} with unit {unit}")
            
            # For complex polygons, fill the holes to ensure we catch all units
            if has_holes:
                # Create a copy of the feature with filled holes
                filled_feature = parcel_data["feature"].copy()
                filled_feature["geometry"] = PolygonHoleFiller.fill_polygon_holes(geometry)
                
                # Create a new shape from the filled geometry
                filled_shape = shape(filled_feature["geometry"])
                
                # Store both the original and filled versions
                complex_polygons[parcel_id] = {
                    "original": parcel_data,
                    "filled": {
                        "feature": filled_feature,
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
    
    # Get address info for matching
    complex_house = complex_data["filled"]["properties"].get("HOUSE")
    complex_street = complex_data["filled"]["properties"].get("STREET")
    
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
                "mostly_overlapping": False,
                "address_match": False
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
            
            # Check address match
            parcel_house = parcel_data["properties"].get("HOUSE")
            parcel_street = parcel_data["properties"].get("STREET")
            
            # Match address if available
            if (parcel_house and complex_house and parcel_house == complex_house and
                parcel_street and complex_street and parcel_street == complex_street):
                criteria["address_match"] = True
            
            # Check if unit ID suggests it's part of the complex
            is_unit = False
            parcel_unit = parcel_data["properties"].get("UNIT")
            if parcel_unit and parcel_unit not in ("CMPLX", "COMPLEX", "COMMON"):
                # It has a unit number, likely part of the complex
                is_unit = True
            
            # Check if it's part of an address range
            # e.g., if complex is 8068 and parcels are 8050-8099, they're probably related
            address_range_match = False
            if parcel_house and complex_house and parcel_street == complex_street:
                try:
                    ph = int(parcel_house)
                    ch = int(complex_house)
                    if abs(ph - ch) < 50:  # Arbitrary threshold for nearby addresses
                        address_range_match = True
                except (ValueError, TypeError):
                    pass
            
            # Determine if this parcel should be included
            should_include = (
                (criteria["centroid_inside"] or criteria["mostly_overlapping"]) and 
                (criteria["address_match"] or address_range_match or is_unit)
            )
            
            # Record this for analysis
            all_found.append({
                "parcel_id": parcel_id,
                "criteria": criteria,
                "is_unit": is_unit,
                "address_range_match": address_range_match,
                "included": should_include
            })
            
            # If it meets our criteria, include it
            if should_include:
                contained_parcels.append(parcel_id)
    
    # Debug: analyze what we found
    if len(all_found) > 0:
        print(f"Complex {complex_id}: Found {len(contained_parcels)} contained parcels out of {len(all_found)} candidates")
    
    return contained_parcels

def create_description_mapping(tax_rates):
    """Create a mapping from Description values to tax types."""
    mapping = {}
    for entry in tax_rates:
        # Each entry has "code" and "description" fields
        mapping[entry["description"]] = entry["code"]
    return mapping

def load_geosource_assessments(file_path):
    """Load the geosource assessments JSON file and create a lookup by FOLIO."""
    data = JsonFileLoader.load_json_file(file_path)
    
    # Process into a dictionary keyed by FOLIO
    lookup = {}
    
    # Check if it's a list structure (not wrapped in a FeatureCollection)
    if isinstance(data, list):
        items = data
    # Otherwise it might be within a structure with features or items
    elif "features" in data:
        items = data["features"]
    elif "items" in data:
        items = data["items"]
    else:
        # Try to find the actual data
        if isinstance(data, dict) and len(data) == 1:
            # Might be a single key with the actual list
            first_key = list(data.keys())[0]
            if isinstance(data[first_key], list):
                items = data[first_key]
            else:
                items = []
        else:
            items = []
    
    # Extract attributes
    for item in items:
        if "attributes" in item:
            attributes = item["attributes"]
            folio = attributes.get("FOLIO")
            if folio:
                lookup[folio] = attributes
    
    print(f"Loaded {len(lookup)} assessments from geosource file")
    return lookup

def main():
    # Load the data files
    print("Loading data files...")
    tax_data = JsonFileLoader.load_json_file("2024_tax_rates.json")
    # Assuming tax_rates is in the format described
    tax_rates = [
        {"code": "Residential", "description": "Residential", "rate": 3.17967},
        {"code": "Utilities", "description": "Utilities", "rate": 38.86440},
        {"code": "MajorIndustry", "description": "Major Industry", "rate": 6.00953},
        {"code": "LightIndustry", "description": "Light Industry", "rate": 7.44413},
        {"code": "Business", "description": "Business Other", "rate": 9.61864},
        {"code": "Rec_NonProfit", "description": "Recreation", "rate": 5.13217},
        {"code": "Farm", "description": "Farm", "rate": 17.99}
    ]
    
    assessments_geojson = JsonFileLoader.load_json_file("2024_assessments.geojson")
    parcels_geojson = JsonFileLoader.load_json_file("2024_parcels.geojson")
    
    # Load the geosource assessments file
    geosource_assessments = load_geosource_assessments("geosource_assessments_2024.json")
    
    # Create mapping from Description values to tax types
    description_mapping = create_description_mapping(tax_rates)
    
    # Process polygon parcels first
    parcels_with_shapes = {}
    unified_parcels = {}
    
    # Process parcels from the parcels file
    print("Processing polygon parcels...")
    for feature in tqdm(parcels_geojson["features"], desc="Processing parcels"):
        parcel_id = feature["properties"].get("FOLIO")
        
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
            
            # Store the parcel info along with its shape
            parcels_with_shapes[parcel_id] = {
                "feature": feature,
                "shape": parcel_shape,
                "properties": feature["properties"],
                "combined_units": 0,
                "taxable_value": 0,
                "assessed_value": 0,  # Initialize with 0
                "combined_parcels": []
            }
            
            # Initialize in unified parcels too
            unified_parcels[parcel_id] = {
                "geometry": feature["geometry"],
                "properties": feature["properties"],
                "combined_units": 0,
                "combined_parcels": [],
                "taxable_value": 0,
                "assessed_value": 0,  # Initialize with 0
                "parent_parcel": None  # Will be set if this parcel is inside another
            }
        except Exception as e:
            print(f"Error processing parcel {parcel_id}: {e}")
    
    print(f"Processed {len(parcels_with_shapes)} valid parcels")
    
    # First pass: Create spatial index for efficient overlap checking
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
    
    # Create a spatial index for assessments too
    print("Creating spatial index for assessments...")
    assessment_idx = index.Index()
    assessment_points = {}
    
    for idx, feature in enumerate(assessments_geojson["features"]):
        point_coords = GeometryPointExtractor.extract_point_from_geometry(feature["geometry"])
        if point_coords:
            assessment_points[idx] = {
                "folio": feature["properties"].get("Folio"),
                "coords": point_coords,
                "properties": feature["properties"]
            }
            # Add to spatial index with a small buffer
            minx, miny = point_coords[0] - 0.001, point_coords[1] - 0.001
            maxx, maxy = point_coords[0] + 0.001, point_coords[1] + 0.001
            assessment_idx.insert(idx, (minx, miny, maxx, maxy))
    
    # Match assessment data with parcels and calculate tax
    matched_count = 0
    unmatched_count = 0
    apartment_count = 0
    direct_match_count = 0
    geosource_match_count = 0
    
    # Process assessments in batches for better performance
    batch_size = 1000
    total_assessments = len(assessments_geojson["features"])
    num_batches = (total_assessments + batch_size - 1) // batch_size
    
    print(f"Processing {total_assessments} assessments in {num_batches} batches...")
    
    # First pass: Match direct FOLIO matches and store their IDs
    processed_folios = set()
    
    for batch_num in tqdm(range(num_batches), desc="Processing assessment batches"):
        start_idx = batch_num * batch_size
        end_idx = min(start_idx + batch_size, total_assessments)
        
        for i in range(start_idx, end_idx):
            feature = assessments_geojson["features"][i]
            properties = feature["properties"]
            folio = properties.get("Folio")
            
            if folio is None or folio == "" or folio in processed_folios:
                continue
            
            # Mark this folio as processed to avoid duplicates
            processed_folios.add(folio)
            
            # Calculate the taxable value and assessed value
            taxable_value = TaxableValueCalculator.calculate_taxable_value(properties, tax_rates, description_mapping)
            assessed_value = TotalAssessedValueCalculator.calculate_total_assessed_value(properties)
            
            # Check if there's a direct FOLIO match in our parcels (regular case)
            if folio in unified_parcels:
                unified_parcels[folio]["taxable_value"] = taxable_value
                unified_parcels[folio]["assessed_value"] = assessed_value
                
                # If this parcel is contained within another, add the value to the parent as well
                parent_id = unified_parcels[folio]["parent_parcel"]
                if parent_id:
                    if parent_id in unified_parcels:
                        unified_parcels[parent_id]["taxable_value"] += taxable_value
                        unified_parcels[parent_id]["assessed_value"] += assessed_value
                        unified_parcels[parent_id]["combined_units"] += 1
                        if folio not in unified_parcels[parent_id]["combined_parcels"]:
                            unified_parcels[parent_id]["combined_parcels"].append(folio)
                
                direct_match_count += 1
                matched_count += 1
            else:
                # Special handling for complex developments
                # Check if this is one of the units in a complex development
                for complex_id, unit_ids in container_parcels.items():
                    if folio in unit_ids and complex_id in unified_parcels:
                        # Add the value to the complex
                        unified_parcels[complex_id]["taxable_value"] += taxable_value
                        unified_parcels[complex_id]["assessed_value"] += assessed_value
                        unified_parcels[complex_id]["combined_units"] += 1
                        if folio not in unified_parcels[complex_id]["combined_parcels"]:
                            unified_parcels[complex_id]["combined_parcels"].append(folio)
                        matched_count += 1
                        apartment_count += 1
                        break
                else:
                    # This assessment doesn't match any parcel in the primary source
                    # Check the geosource assessments
                    if folio in geosource_assessments:
                        geosource_properties = geosource_assessments[folio]
                        geosource_taxable = TaxableValueCalculator.calculate_taxable_value(geosource_properties, tax_rates, description_mapping)
                        geosource_assessed = TotalAssessedValueCalculator.calculate_total_assessed_value(geosource_properties)
                        
                        # Use these values for the parcel if it exists
                        if folio in unified_parcels:
                            unified_parcels[folio]["taxable_value"] = geosource_taxable
                            unified_parcels[folio]["assessed_value"] = geosource_assessed
                            
                            # If this parcel is contained within another, add the value to the parent as well
                            parent_id = unified_parcels[folio]["parent_parcel"]
                            if parent_id and parent_id in unified_parcels:
                                unified_parcels[parent_id]["taxable_value"] += geosource_taxable
                                unified_parcels[parent_id]["assessed_value"] += geosource_assessed
                                unified_parcels[parent_id]["combined_units"] += 1
                                if folio not in unified_parcels[parent_id]["combined_parcels"]:
                                    unified_parcels[parent_id]["combined_parcels"].append(folio)
                                    
                            geosource_match_count += 1
                            matched_count += 1
                        else:
                            # Still unmatched
                            unmatched_count += 1
                            if unmatched_count <= 10:
                                print(f"Unmatched: Assessment FOLIO={folio} (not found in parcels)")
                    else:
                        # Assessment not found in either source
                        unmatched_count += 1
                        if unmatched_count <= 10:
                            print(f"Unmatched: Assessment FOLIO={folio} (not found in either source)")

    # Now process remaining parcels that had no assessment data from the first pass
    print("Processing parcels without assessment data...")
    remaining_parcel_count = 0
    
    for parcel_id, parcel_data in tqdm(unified_parcels.items(), desc="Checking geosource for remaining parcels"):
        # Skip parcels that already have values
        if parcel_data["taxable_value"] > 0 or parcel_data["assessed_value"] > 0:
            continue
            
        # Skip parcels that are contained within others (they'll get values from their parents)
        if parcel_data["parent_parcel"]:
            continue
            
        # Check if this parcel has data in the geosource file
        if parcel_id in geosource_assessments:
            geosource_properties = geosource_assessments[parcel_id]
            parcel_data["taxable_value"] = TaxableValueCalculator.calculate_taxable_value(geosource_properties, tax_rates, description_mapping)
            parcel_data["assessed_value"] = TotalAssessedValueCalculator.calculate_total_assessed_value(geosource_properties)
            geosource_match_count += 1
            matched_count += 1
            remaining_parcel_count += 1
            
            if remaining_parcel_count <= 10:
                print(f"Found in geosource: Parcel FOLIO={parcel_id}")

    # For any container parcel that didn't get assessments directly, 
    # ensure it has the combined values from all its contained parcels
    for container_id, contained_ids in container_parcels.items():
        if container_id in unified_parcels:
            container_data = unified_parcels[container_id]
            
            # Initialize if not already done
            if "combined_units" not in container_data or container_data["combined_units"] == 0:
                container_data["combined_units"] = len(contained_ids)
                
            if "combined_parcels" not in container_data or not container_data["combined_parcels"]:
                container_data["combined_parcels"] = contained_ids
                
            # Sum up taxable values from contained parcels if needed
            if container_data["taxable_value"] == 0:
                for child_id in contained_ids:
                    if child_id in unified_parcels:
                        container_data["taxable_value"] += unified_parcels[child_id]["taxable_value"]
                        container_data["assessed_value"] += unified_parcels[child_id]["assessed_value"]
    
    # Create parallel processing capabilities
    print("Second pass optimization: Creating spatial lookups...")
    
    # Create a dictionary for fast FOLIO lookup in unified_parcels
    unified_folio_lookup = set(unified_parcels.keys())
    
    # Process unmatched assessments using spatial lookup
    print("Matching unmatched assessments by location...")
    second_pass_matches = 0
    
    # Create a list of unmatched assessments to process
    unmatched_assessments = []
    
    for i, feature in enumerate(assessments_geojson["features"]):
        properties = feature["properties"]
        folio = properties.get("Folio")
        
        if folio is None or folio == "" or folio in unified_folio_lookup:
            continue  # Skip if already matched
        
        # Check if it's in geosource
        if folio in geosource_assessments:
            properties = geosource_assessments[folio]  # Use geosource properties
        
        point_coords = GeometryPointExtractor.extract_point_from_geometry(feature["geometry"])
        if point_coords:
            unmatched_assessments.append({
                "index": i,
                "folio": folio,
                "coords": point_coords,
                "taxable_value": TaxableValueCalculator.calculate_taxable_value(properties, tax_rates, description_mapping),
                "assessed_value": TotalAssessedValueCalculator.calculate_total_assessed_value(properties)
            })
    
    print(f"Processing {len(unmatched_assessments)} unmatched assessments...")
    
    # Process in smaller batches
    batch_size = 100
    for i in tqdm(range(0, len(unmatched_assessments), batch_size), desc="Matching by location"):
        batch = unmatched_assessments[i:i+batch_size]
        
        for assessment in batch:
            folio = assessment["folio"]
            point_coords = assessment["coords"]
            taxable_value = assessment["taxable_value"]
            assessed_value = assessment["assessed_value"]
            
            # Use a buffer for point search
            minx, miny = point_coords[0] - 0.01, point_coords[1] - 0.01
            maxx, maxy = point_coords[0] + 0.01, point_coords[1] + 0.01
            
            point = Point(point_coords)
            potential_parcels = list(spatial_idx.intersection((minx, miny, maxx, maxy)))
            found_match = False
            
            for idx in potential_parcels:
                parcel_id = idx_to_parcel_id.get(idx)
                if parcel_id and parcel_id in parcels_with_shapes:
                    if parcels_with_shapes[parcel_id]["shape"].contains(point):
                        # Found a containing parcel
                        unified_parcels[parcel_id]["taxable_value"] += taxable_value
                        unified_parcels[parcel_id]["assessed_value"] += assessed_value
                        unified_parcels[parcel_id]["combined_units"] += 1
                        unified_parcels[parcel_id]["combined_parcels"].append(folio)
                        second_pass_matches += 1
                        matched_count += 1
                        found_match = True
                        
                        if second_pass_matches <= 10:
                            print(f"Location match: Assessment FOLIO={folio} combined into Parcel FOLIO={parcel_id}")
                        
                        break
                    
            if not found_match and unmatched_count <= 20:
                print(f"Still unmatched after location check: Assessment FOLIO={folio}")
    
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
                filled_geometry = complex_developments[parcel_id]["filled"]["feature"]["geometry"]
                parcel_data["geometry"] = filled_geometry
        
        # Update the properties to include the combined unit count and taxable value
        parcel_data["properties"]["TaxableValue"] = parcel_data["taxable_value"]
        parcel_data["properties"]["AssessedValue"] = parcel_data["assessed_value"]
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
        "crs": parcels_geojson["crs"],  # Use the original CRS
        "features": output_features
    }
    
    # Write to file
    with open("unified_parcels.geojson", "w", encoding='utf-8') as outfile:
        json.dump(output_geojson, outfile)
    
    print(f"Process complete!")
    print(f"Total parcels processed: {len(parcels_with_shapes)}")
    print(f"Parcels in output file: {len(output_features)}")
    print(f"Parcels excluded (contained within others): {excluded_count}")
    print(f"Direct matches: {direct_match_count}")
    print(f"Geosource matches: {geosource_match_count}")
    print(f"Apartment/units combined: {apartment_count}")
    print(f"Second pass location matches: {second_pass_matches}")
    print(f"Unmatched assessment points: {unmatched_count}")
    print(f"Total matched: {matched_count}")
    print(f"Output written to unified_parcels.geojson")

if __name__ == "__main__":
    main()