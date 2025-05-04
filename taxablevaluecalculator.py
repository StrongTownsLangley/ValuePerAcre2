class TaxableValueCalculator:
    
    def calculate_taxable_value(properties, tax_rates, description_mapping, tax_rate_divider=1000):
    #Calculate taxable value based on tax rates and property assessments.
        tax = 0
        
        # Check if properties has "Description" field, indicating it's from geosource
        if "Description" in properties:
            description = properties["Description"]
            # Find matching tax type for this description
            tax_type = None
            
            # Look up the tax rate by description
            for rate_entry in tax_rates:
                if rate_entry["description"] == description:
                    tax_type = rate_entry["code"]
                    rate = rate_entry["rate"]
                    
                    # Use GrossLand and GrossImprovements as values
                    land_value = float(properties.get("GrossLand", 0))
                    building_value = float(properties.get("GrossImprovements", 0))
                    
                    tax += building_value * (rate / tax_rate_divider)
                    tax += land_value * (rate / tax_rate_divider)
                    break
        else:
            # Handle regular assessment properties
            for rate_entry in tax_rates:
                tax_type = rate_entry["code"]
                rate = rate_entry["rate"]
                
                buildings_key = f"{tax_type}_Buildings"
                land_key = f"{tax_type}_Land"
                
                # Check if the keys exist in properties, otherwise use 0
                buildings_value = float(properties.get(buildings_key, 0))
                land_value = float(properties.get(land_key, 0))
                
                tax += buildings_value * (rate / tax_rate_divider)
                tax += land_value * (rate / tax_rate_divider)
        
        return round(tax, 2)