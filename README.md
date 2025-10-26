# 🌾 Precision Agriculture Analysis Toolkit

A comprehensive Python toolkit for precision agriculture analysis using satellite imagery and remote sensing data. This project provides both **analysis libraries** and **educational notebooks** to help farmers, researchers, and agronomists make data-driven decisions.

## 🎯 Project Vision

This toolkit aims to provide accessible, open-source tools for:
- 📊 **Agricultural Monitoring**: Track crop health, growth stages, and field conditions
- 🛰️ **Remote Sensing Analysis**: Process and analyze satellite imagery at scale
- 💧 **Resource Management**: Optimize irrigation, fertilization, and other inputs
- 📈 **Yield Prediction**: Forecast crop yields using multi-temporal analysis
- 🌍 **Environmental Impact**: Monitor soil health, carbon sequestration, and sustainability metrics

## 🚀 Current Capabilities

### Vegetation Index Analysis (Available Now)
Compute multiple vegetation indices from free satellite data (Sentinel-2, Landsat) to monitor crop health, vegetation vigor, and plant stress.

**Available Indices:**
- **NDVI** - Normalized Difference Vegetation Index (general vegetation health)
- **EVI** - Enhanced Vegetation Index (dense vegetation, atmospheric correction)
- **SAVI** - Soil Adjusted Vegetation Index (sparse vegetation, early season)
- **NDRE** - Normalized Difference Red Edge (chlorophyll, nitrogen status) *Sentinel-2 only*
- **GNDVI** - Green NDVI (photosynthetic activity, nitrogen)

**Key Features:**
- Automatic data retrieval via STAC APIs
- Support for multiple AOI formats (GeoJSON, bounding boxes, coordinates)
- Efficient multi-index computation (loads each band once)
- Built-in side-by-side visualizations and statistics
- Both CLI and Python library interfaces
- Interactive Jupyter tutorials for learning

## 🔮 Coming Soon

- **Additional Indices**: NDWI, NDMI, MSI, and more
- **Time Series Analysis**: Track changes over growing seasons
- **Yield Modeling**: Predictive analytics for crop production
- **Soil Moisture Estimation**: Remote sensing-based soil monitoring
- **Field Boundary Detection**: Automated field delineation
- **Crop Classification**: Machine learning-based crop type identification

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/chris/precision-ag.git
cd precision-ag

# Create and activate a virtual environment (recommended)
# Option 1: Using uv (fastest)
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Option 2: Using virtualenv
virtualenv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install the package
pip install -e .

# Or with optional dependencies for development and notebooks
pip install -e ".[dev,notebook]"
```

### Usage Examples

#### Simple NDVI Analysis

**Command Line:**
```bash
# Quick NDVI computation using a bounding box
ndvi-compute \
    --aoi "-121.0,37.0,-120.5,37.5" \
    --start-date 2024-06-01 \
    --end-date 2024-06-30
```

**Python Library:**
```python
from precision_ag.ndvi_from_aoi import compute_ndvi_for_aoi

results = compute_ndvi_for_aoi(
    aoi_input=[-121.0, 37.0, -120.5, 37.5],
    start_date="2024-06-01",
    end_date="2024-06-30"
)
```

#### Multi-Index Vegetation Analysis

**Command Line:**
```bash
# Compute multiple vegetation indices at once
vegetation-index-compute \
    --aoi my_field.geojson \
    --start-date 2024-06-01 \
    --end-date 2024-06-30 \
    --indices ndvi evi savi ndre gndvi

# List available indices for your satellite collection
vegetation-index-compute --list-indices
```

**Python Library:**
```python
from precision_ag.vegetation_indices import compute_vegetation_indices_for_aoi

# Compute multiple indices efficiently (loads each band once)
results = compute_vegetation_indices_for_aoi(
    aoi_input="field.geojson",
    start_date="2024-06-01",
    end_date="2024-06-30",
    indices=['ndvi', 'evi', 'savi', 'ndre', 'gndvi'],
    output_dir="my_analysis",
    visualize=True  # Creates side-by-side comparison plots
)
```

### Interactive Tutorials

Launch the Jupyter notebooks to learn interactively:

```bash
# Tutorial 1: Introduction to NDVI
jupyter notebook notebooks/NDVI_Tutorial.ipynb

# Tutorial 2: Comprehensive Vegetation Health (coming soon)
# Multiple indices, comparison, and decision making
```

**Documentation:**
- [NDVI Module](precision-ag/ndvi_from_aoi.py) - Simple NDVI computation
- [Vegetation Indices Module](precision-ag/vegetation_indices.py) - Multi-index analysis
- Run `ndvi-compute --help` or `vegetation-index-compute --help` for CLI options

## 📁 Project Structure

```
precision-ag/
├── precision-ag/                      # Analysis libraries
│   ├── ndvi_from_aoi.py              # Simple NDVI computation
│   ├── vegetation_indices.py         # Multi-index analysis (NDVI, EVI, SAVI, NDRE, GNDVI)
│   └── ...                            # More modules coming soon
├── notebooks/                         # Educational tutorials
│   ├── NDVI_Tutorial.ipynb           # Tutorial 1: NDVI basics
│   └── Vegetation_Health_Tutorial.ipynb  # Tutorial 2: Multi-index analysis (coming soon)
├── pyproject.toml                     # Project configuration
├── README.md                          # This file
└── output/                            # Generated outputs (auto-created)
```

## 🔧 Requirements

- Python >= 3.8
- numpy >= 1.20.0
- rasterio >= 1.3.0
- matplotlib >= 3.5.0
- pystac-client >= 0.7.0
- planetary-computer >= 1.0.0
- requests >= 2.28.0
- shapely >= 2.0.0

## 🌍 Data Sources

This tool uses free, public satellite data via STAC APIs:

| Satellite | Resolution | Revisit | Coverage | STAC Catalog |
|-----------|-----------|---------|----------|--------------|
| **Sentinel-2** | 10m | 5 days | Global | [Microsoft Planetary Computer](https://planetarycomputer.microsoft.com/) |
| **Landsat 8/9** | 30m | 16 days | Global | [Earth Search AWS](https://earth-search.aws.element84.com/) |

## 💡 Use Cases

This toolkit supports various precision agriculture applications:

- 🌾 **Crop Health Monitoring**: Detect stress, disease, and nutrient deficiencies early
- 💧 **Irrigation Management**: Optimize water usage based on vegetation and soil data
- 🌱 **Growth Stage Tracking**: Monitor crop development throughout the season
- 📊 **Yield Forecasting**: Predict harvest outcomes using multi-temporal analysis
- 🗺️ **Field Mapping**: Delineate management zones for variable rate applications
- 🌍 **Sustainability Reporting**: Track environmental metrics and carbon footprint
- 🔬 **Research & Development**: Support agricultural research with reproducible analysis

## 🛠️ Development

### Running Tests
```bash
pytest
```

### Code Formatting
```bash
black precision-ag/
```

### Type Checking
```bash
mypy precision-ag/
```

## 📝 Outputs

Analysis tools generate georeferenced raster files (GeoTIFF), visualizations (PNG/PDF), and statistical summaries. All outputs are compatible with standard GIS software (QGIS, ArcGIS) and can be used for further analysis or reporting.

## 📚 Learning & Documentation

This project emphasizes both **practical tools** and **educational resources**:

- **Analysis Libraries**: Production-ready Python modules for data processing
- **Tutorial Notebooks**: Interactive Jupyter notebooks explaining concepts, methods, and best practices
- **Documentation**: Inline code documentation and detailed docstrings
- **Examples**: Real-world use cases demonstrating various agricultural scenarios

Whether you're a researcher, agronomist, or developer, you'll find resources suited to your needs.

## 🤝 Contributing

Contributions are welcome! This project is in active development. Areas where you can help:
- Adding new analysis modules (vegetation indices, soil metrics, yield models)
- Creating educational notebooks and tutorials
- Improving documentation and examples
- Bug fixes and performance improvements

Please feel free to submit a Pull Request or open an issue to discuss new features.

## 📄 License

MIT License - see LICENSE file for details

## 🙏 Acknowledgments

- [Microsoft Planetary Computer](https://planetarycomputer.microsoft.com/) - Free satellite data access
- [STAC](https://stacspec.org/) - Standardized geospatial data catalog
- [ESA Copernicus](https://www.copernicus.eu/) - Sentinel missions
- [NASA/USGS](https://www.usgs.gov/landsat-missions) - Landsat program

## 📧 Contact

For questions or issues, please open an issue on GitHub.

---

**Building the future of precision agriculture, one pixel at a time 🌾🛰️**

