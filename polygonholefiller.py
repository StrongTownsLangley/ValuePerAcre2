from shapely.geometry import Polygon, shape, mapping
import shapely.ops

class PolygonHoleFiller:
    
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