from polygonholefiller import PolygonHoleFiller
from shapely.geometry import shape

class ComplexPolygonsIdentifier:

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