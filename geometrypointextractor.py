from shapely.geometry import shape

class GeometryPointExtractor:

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