class TotalAssessedValueCalculator:

    def calculate_total_assessed_value(properties):
        """Calculate the total assessed value (sum of all building and land values)."""
        # Check if properties has "Description" field, indicating it's from geosource
        if "Description" in properties:
            # Use GrossLand and GrossImprovements from geosource
            land_value = float(properties.get("GrossLand", 0))
            building_value = float(properties.get("GrossImprovements", 0))
            return round(land_value + building_value, 2)
        else:
            # Original logic for regular assessment properties
            total_value = 0
            
            # Look for any property keys that end with _Buildings or _Land
            for key, value in properties.items():
                if key.endswith('_Buildings') or key.endswith('_Land'):
                    try:
                        total_value += float(value)
                    except (ValueError, TypeError):
                        # Skip if the value can't be converted to float
                        pass
            
            return round(total_value, 2)