class ParcelsWithinComplexFinder:

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