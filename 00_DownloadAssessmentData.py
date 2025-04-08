import requests
import json
import time
import os
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from tqdm import tqdm

def download_assessment_data(assessment_url, parcels_url, output_file, year=2024, max_records_per_request=2000, use_cache=False, cache_file=None):
    """
    Downloads assessment data records for a specific year from an ArcGIS REST service,
    looks up FOLIO information using ParcelKey, and saves the combined data to a JSON file.
    
    Args:
        assessment_url (str): Base URL of the ArcGIS REST service for assessment data
        parcels_url (str): Base URL of the ArcGIS REST service for parcels data
        output_file (str): Path to the output JSON file
        year (int): Year to filter records by (default: 2024)
        max_records_per_request (int): Maximum number of records per request (default: 2000)
        use_cache (bool): Whether to use cached assessment data if available
        cache_file (str): Path to the cache file for assessment data
    """
    # Ensure URLs end with a forward slash
    if not assessment_url.endswith('/'):
        assessment_url += '/'
    if not parcels_url.endswith('/'):
        parcels_url += '/'
    
    # Construct the query URLs
    assessment_query_url = assessment_url + "query"
    parcels_query_url = parcels_url + "query"
    
    # Get assessment data (either from cache or by downloading)
    all_records = []
    
    if use_cache and cache_file and os.path.exists(cache_file):
        print(f"Loading assessment data from cache file: {cache_file}")
        try:
            with open(cache_file, 'r') as f:
                cached_data = json.load(f)
                all_records = cached_data.get('features', [])
            print(f"Successfully loaded {len(all_records)} records from cache.")
        except Exception as e:
            print(f"Error loading cache file: {e}")
            print("Will download assessment data instead.")
            use_cache = False
    
    if not use_cache or not all_records:
        # Prepare query parameters for assessment data
        params = {
            'where': f'Year = {year}',
            'outFields': '*',  # Get all fields
            'returnGeometry': 'false',  # No geometry needed for table data
            'f': 'json',  # Request format
            'orderByFields': 'OBJECTID',  # Order by OBJECTID for pagination
        }
        
        # Download assessment data with a progress bar
        all_records = []
        offset = 0
        total_records = None
        
        print(f"Downloading assessment data for year {year}...")
        
        with tqdm(total=None, desc="Downloading assessment records", unit="records") as pbar:
            while total_records is None or offset < total_records:
                # Add pagination parameters
                params['resultOffset'] = offset
                params['resultRecordCount'] = max_records_per_request
                
                # Make the request
                try:
                    response = requests.get(assessment_query_url, params=params)
                    response.raise_for_status()  # Raise an exception for HTTP errors
                    
                    data = response.json()
                    
                    # Check for errors in the response
                    if 'error' in data:
                        print(f"Error: {data['error']['message']}")
                        return False
                    
                    # Get the records from the response
                    records = data.get('features', [])
                    
                    # Update total records count if not yet known
                    if total_records is None:
                        total_records = data.get('exceededTransferLimit', False)
                        if total_records:
                            # If transfer limit exceeded, we don't know the exact count yet
                            total_records = float('inf')
                        else:
                            # If not exceeded, total is the current count of records
                            total_records = len(records)
                            # Update progress bar total
                            pbar.total = total_records
                            pbar.refresh()
                    
                    # If no records, we're done
                    if not records:
                        break
                    
                    # Add records to our collection
                    all_records.extend(records)
                    
                    # Update offset for the next request
                    offset += len(records)
                    
                    # Update progress bar
                    pbar.update(len(records))
                    pbar.set_postfix({"Total": len(all_records)})
                    
                    # If there are no more records (exceededTransferLimit is False or not present), we're done
                    if not data.get('exceededTransferLimit', False):
                        break
                        
                    # Small delay to avoid overwhelming the server
                    time.sleep(0.5)
                    
                except requests.exceptions.RequestException as e:
                    print(f"Request failed: {e}")
                    return False
        
        print(f"Successfully downloaded {len(all_records)} assessment records for year {year}.")
        
        # Save to cache file if provided
        if cache_file:
            print(f"Saving assessment data to cache file: {cache_file}")
            with open(cache_file, 'w') as f:
                json.dump({'features': all_records}, f, indent=2)
    
    # Extract unique ParcelKeys from assessment records
    parcel_keys = set()
    for record in all_records:
        parcel_key = record.get('attributes', {}).get('ParcelKey')
        if parcel_key:
            parcel_keys.add(parcel_key)
    
    print(f"Found {len(parcel_keys)} unique ParcelKeys.")
    
    # Function to get FOLIO for a single ParcelKey (as a fallback)
    def get_folio_for_key(parcel_key):
        params = {
            'where': f"ParcelKey = '{parcel_key}'",
            'outFields': 'ParcelKey,FOLIO',
            'returnGeometry': 'false',
            'f': 'json'
        }
        
        try:
            response = requests.get(parcels_query_url, params=params)
            response.raise_for_status()
            data = response.json()
            
            if 'error' in data:
                return None
            
            # Extract FOLIO from the response
            features = data.get('features', [])
            if features and 'attributes' in features[0]:
                attrs = features[0]['attributes']
                if 'FOLIO' in attrs:
                    return attrs['FOLIO']
            
            return None
            
        except requests.exceptions.RequestException as e:
            if '404' in str(e):
                return get_folio_via_relationships(parcel_key)
            else:
                return None

    # Function to get FOLIOs for a batch of ParcelKeys
    def get_folios_for_batch(parcel_keys_batch):
        # Use OR to combine multiple ParcelKey conditions
        where_clause = " OR ".join([f"ParcelKey = '{key}'" for key in parcel_keys_batch])
        params = {
            'where': where_clause,
            'outFields': 'ParcelKey,FOLIO',
            'returnGeometry': 'false',
            'f': 'json'
        }
        
        try:
            response = requests.get(parcels_query_url, params=params)
            response.raise_for_status()
            data = response.json()
            
            if 'error' in data:
                # If error, fall back to individual queries
                result = {}
                for key in parcel_keys_batch:
                    folio = get_folio_for_key(key)
                    if folio:
                        result[key] = folio
                return result
            
            # Create a mapping of ParcelKey to FOLIO
            parcel_to_folio = {}
            for feature in data.get('features', []):
                attrs = feature.get('attributes', {})
                if 'ParcelKey' in attrs and 'FOLIO' in attrs:
                    parcel_to_folio[attrs['ParcelKey']] = attrs['FOLIO']
            
            return parcel_to_folio
            
        except requests.exceptions.RequestException as e:
            if '404' in str(e):
                # If 404, fall back to individual queries
                result = {}
                for key in parcel_keys_batch:
                    folio = get_folio_for_key(key)
                    if folio:
                        result[key] = folio
                return result
            else:
                return {}
    
    # Try to get FOLIO via relationships endpoint
    def get_folio_via_relationships(parcel_key):
        # First try to find the OBJECTID for this ParcelKey
        params = {
            'where': f"ParcelKey = '{parcel_key}'",
            'outFields': 'OBJECTID',
            'returnGeometry': 'false',
            'f': 'json'
        }
        
        try:
            response = requests.get(parcels_query_url, params=params)
            response.raise_for_status()
            data = response.json()
            
            features = data.get('features', [])
            if not features:
                return None
                
            object_id = features[0].get('attributes', {}).get('OBJECTID')
            if not object_id:
                return None
                
            # Now use the related records endpoint to get the FOLIO
            related_url = f"{parcels_url}queryRelatedRecords"
            related_params = {
                'relationshipId': 2,  # The Assessment relationship from the documentation
                'objectIds': object_id,
                'outFields': 'FOLIO',
                'f': 'json'
            }
            
            related_response = requests.get(related_url, params=related_params)
            related_response.raise_for_status()
            related_data = related_response.json()
            
            related_records = related_data.get('relatedRecordGroups', [])
            if related_records and 'relatedRecords' in related_records[0]:
                for record in related_records[0]['relatedRecords']:
                    if 'attributes' in record and 'FOLIO' in record['attributes']:
                        return record['attributes']['FOLIO']
            
            return None
            
        except requests.exceptions.RequestException:
            return None
    
    # Check if the user wants to look up FOLIO information
    lookup_folio = True
    if len(all_records) > 0:
        response = input("Do you want to look up FOLIO information for each ParcelKey? (y/n): ").lower()
        lookup_folio = response.startswith('y')
    
    parcel_to_folio = {}
    
    if lookup_folio:
        # Get FOLIO information for all ParcelKeys in batches
        print("Looking up FOLIO information for each ParcelKey (batch mode)...")
        
        # Convert set to list for batching
        parcel_keys_list = list(parcel_keys)
        total_keys = len(parcel_keys_list)
        
        # Define batch size - try 10 at a time instead of individual requests
        batch_size = 10
        
        # Create batches
        batches = [parcel_keys_list[i:i+batch_size] for i in range(0, total_keys, batch_size)]
        total_batches = len(batches)
        
        # Use ThreadPoolExecutor to parallelize batch requests
        parcel_to_folio = {}
        
        # Create a progress bar for batches
        with tqdm(total=total_batches, desc="Retrieving FOLIO IDs (batches)", unit="batch") as pbar:
            # Process batches with threading for better performance
            with ThreadPoolExecutor(max_workers=3) as executor:
                future_to_batch = {executor.submit(get_folios_for_batch, batch): i for i, batch in enumerate(batches)}
                
                for future in concurrent.futures.as_completed(future_to_batch):
                    batch_results = future.result()
                    parcel_to_folio.update(batch_results)
                    
                    # Update the progress bar
                    pbar.update(1)
                    pbar.set_postfix({"FOLIOs found": len(parcel_to_folio)})
                    
                    # Small delay to avoid overwhelming the server
                    time.sleep(0.05)
        
        print(f"Retrieved FOLIO information for {len(parcel_to_folio)} ParcelKeys.")
    
        # Add FOLIO information to assessment records
        print("Adding FOLIO information to assessment records...")
        
        missing_folio_count = 0
        with tqdm(total=len(all_records), desc="Adding FOLIO to records", unit="records") as pbar:
            for record in all_records:
                parcel_key = record.get('attributes', {}).get('ParcelKey')
                if parcel_key and parcel_key in parcel_to_folio:
                    record['attributes']['FOLIO'] = parcel_to_folio[parcel_key]
                else:
                    missing_folio_count += 1
                    record['attributes']['FOLIO'] = None
                pbar.update(1)
        
        if missing_folio_count > 0:
            print(f"Warning: Could not find FOLIO information for {missing_folio_count} records ({missing_folio_count/len(all_records)*100:.1f}%).")
    else:
        print("Skipping FOLIO lookup.")
    
    # Save all records to the output file
    print(f"Writing {len(all_records)} records to {output_file}...")
    
    with open(output_file, 'w') as f:
        json.dump({'features': all_records}, f, indent=2)
    
    print(f"Output saved to: {os.path.abspath(output_file)}")
    
    return True

if __name__ == "__main__":
    # URLs from the provided documentation
    assessment_url = "https://mapsvr.tol.ca/arcgisext02/rest/services/GeoSource/DynamicServices/MapServer/1129"
    parcels_url = "https://mapsvr.tol.ca/arcgisext02/rest/services/GeoSource/DynamicServices/MapServer/3"
    
    # Output files
    output_file = "geosource_assessments_2024.json"
    cache_file = "assessments_cache_2024.json"
    
    # Check if cache file exists
    use_cache = False
    if os.path.exists(cache_file):
        response = input(f"Cache file {cache_file} exists. Use cached assessment data? (y/n): ").lower()
        use_cache = response.startswith('y')
    
    # Download the data with FOLIO information
    download_assessment_data(assessment_url, parcels_url, output_file, use_cache=use_cache, cache_file=cache_file)