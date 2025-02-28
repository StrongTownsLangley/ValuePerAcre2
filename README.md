# ValuePerAcre2

A modern Python implementation of the Value-per-Acre analysis tool, designed to create interactive maps for visualizing property tax revenue and land use efficiency.

![image](https://github.com/user-attachments/assets/35a498ae-e725-4642-a2df-676d9f3dbb92)

## What is Value-per-Acre and Why?

Value-per-Acre analysis is a method advocated by [Strong Towns](https://strongtowns.org) to evaluate the efficiency of land use by focusing on both its productivity and desirability relative to infrastructure investment.

Looking at a map with properties grouped and colorized based on their taxable value allows us to see which areas:
- Are the most and least productive
- Contribute the most and least revenue to city finances
- May be inefficiently using valuable land

While higher density development typically corresponds with higher taxable value, more desirable areas also attract higher home prices and, consequently, higher taxable value. This effect can compensate for areas of low taxable value like parks that increase overall desirability of neighboring properties.

## Demo

This tool was used to compile the map at https://strongtownslangley.org/maps?revenue-map

## Features

### Polygon Mode (Recommended)
Parcel data is used to assign colors to each individual parcel based on its value per acre.

![image](https://github.com/user-attachments/assets/f7571b1a-b58d-45d6-8aec-bf763ae56349)


### Point Mode
Each block (by default 100m²) is assigned a value based on the properties it contains.

## How to Use

This Python implementation consists of a pipeline of scripts that:

1. Create a unified parcels file with tax values from assessment data
2. Convert coordinate systems if needed
3. Generate the web map and associated files

### Quick Start

#### Prerequisites
1. Install Python 3.6 or higher from [python.org](https://python.org)
2. Install the required libraries using pip:

```bash
pip install shapely rtree tqdm pyproj geojson
```

Note: For rtree installation on Windows, you may need to install the [Microsoft C++ Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/).

#### Run the Pipeline
For Township of Langley data:
```bash
# Run the full pipeline
python3 01_CreateUnifiedFile_Langley.py
python3 02_GeoJsonConvertToWSG.py unified_parcels.geojson unified_parcels_langley_4326.geojson 
python3 03_GenerateWebMap.py --input unified_parcels_langley_4326.geojson --output-folder langley_json --levels 50 --mode polygons
```

For Moldovan data (or other data already in the right coordinate system):
```bash
python3 01_CreateUnifiedFile_Moldovan.py
python3 03_GenerateWebMap.py --input moldovan_unified_parcels.geojson --output-folder moldovan_json --levels 50 --mode polygons
```

### Understanding the Pipeline

#### Step 1: Create Unified Parcels File
The first script (`01_CreateUnifiedFile_*.py`) processes raw assessment and parcel data to create a unified GeoJSON file with tax values attached to each parcel.

```bash
python3 01_CreateUnifiedFile_Langley.py
```

#### Step 2: Coordinate System Conversion (if needed)
If your data isn't already in WGS84 (EPSG:4326) format needed for web maps, use the conversion script:

```bash
python3 02_GeoJsonConvertToWSG.py input.geojson output_4326.geojson
```

#### Step 3: Generate Web Map
The final script takes the unified parcels file and generates the web map and all supporting files:

```bash
python3 03_GenerateWebMap.py --input unified_parcels_langley_4326.geojson --output-folder output_json --levels 50 --mode polygons
```

Parameters:
- `--input` or `-i`: Path to your unified parcels GeoJSON file (required)
- `--output-folder` or `-o`: Where to save the output files (default: "json")
- `--levels` or `-l`: Number of levels to group the data into (default: 50)
- `--block-size` or `-b`: Size of blocks in meters when using points mode (default: 100.0)
- `--mode` or `-m`: Either "points" or "polygons" (default: "polygons")
  - "points" mode treats parcels as points (centroids) and groups them into blocks
  - "polygons" mode uses the actual parcel geometries for visualization

### Output Files

The tool generates several files in the output folder:

1. Level-specific GeoJSON files (`level_0.json`, `level_1.json`, etc.)
2. A level info summary file (`level_info.json`)
3. Two HTML visualization files:
   - `website.static.html`: Contains embedded data (larger file size but works offline)
   - `website.dynamic.html`: Loads data from the JSON files (smaller file size but requires a web server)

The `website.dynamic.html` is the preferred method for deployment on the internet as it allows you to customize the page and simply update the level files separately in the future.

## Data Format Examples

### Unified Parcels GeoJSON Format

This is the output of the first script and input to the third:

```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "id": "0020929006",
      "geometry": {
        "type": "Polygon",
        "coordinates": [...]
      },
      "properties": {
        "FOLIO": "0020929006",
        "PID": "000-561-291",
        "STREET": "204A ST",
        "HOUSE": 2144,
        "CITY": "LANGLEY",
        "TaxableValue": 7157.44,
        "CombinedUnits": 0,
        "CombinedParcels": ""
      }
    },
    ...
  ]
}
```

### Raw Tax Rates Format

If you're creating your own unified parcels file, the tax rates should be in this format:

```json
{
    "Residential": 3.80248,
    "Utilities": 42.53820,
    "SupportiveHousing": 2.37008,
    "MajorIndustry": 8.34595,
    "LightIndustry": 10.24375,
    "Business": 12.02044,
    "ManagedForest": 0,
    "Rec_NonProfit": 6.73828,
    "Farm": 15.97448
}
```

### Raw Assessments Format

Simplified example of the assessments file structure:

```json
{
"type": "FeatureCollection",
"features": [
  { 
    "type": "Feature", 
    "properties": { 
      "Residential_Buildings": 113700, 
      "Residential_Land": 0, 
      "Latitude": 49.004443779630002, 
      "Longitude": -122.64299023049,
      "PID": "000-561-123"
    } 
  },
  ...
]
}
```

## Other Similar Projects

- A commercial firm producing this kind of analysis is [Urban3](https://www.urbanthree.com/), who produce 3D renderings in addition to factoring in not just the Value but also the Cost-per-Acre.
- Fellow Strong Towns Local Conversation group [A Better Cobb](https://abettercobb.substack.com/) created their own set of Value per Acre tools, available at https://github.com/ABetterCobb/ValuePerAcre/

## Credits and Technical Details

This Python implementation is based on the original C# version (now deprecated) at https://github.com/StrongTownsLangley/ValuePerAcre/

Technologies used:
- Python 3 with libraries:
  - shapely: For geometry operations
  - json: For data processing
  - rtree: For spatial indexing (in the unified file creation)
- [Leaflet](https://github.com/Leaflet/Leaflet): For map visualization

## Help and Contributing

If you encounter any issues or have suggestions for improvement, please [open an issue](https://github.com/YourUsername/ValuePerAcre2/issues) on the GitHub repository. Pull requests are also welcome.

For additional help, join the Strong Towns Langley discord: https://discord.gg/MuAn3cFd8J and ask in the #🧮do-the-math channel.

## License

This program is released under the [Apache 2.0 License](LICENSE). If you use it for your website or project, please provide credit to **Strong Towns Langley** and preferably link to this GitHub.
