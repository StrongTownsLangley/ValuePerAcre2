python3 01_CreateUnifiedFile_Langley.py
python3 02_GeoJsonConvertToWSG.py unified_parcels.geojson unified_parcels_langley_4326.geojson 
python3 03_GenerateWebMap.py --input unified_parcels_langley_4326.geojson --output-folder langley_json --levels 50 --mode polygons