python3 01_CreateUnifiedFile_Moldovan.py
rem python3 02_GeoJsonConvertToWSG.py moldovan_unified_parcels.geojson moldovan_unified_parcels_4326.geojson 
python3 03_GenerateWebMap.py --input moldovan_unified_parcels.geojson --output-folder moldovan_json --levels 50 --mode polygons