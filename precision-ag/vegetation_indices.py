#!/usr/bin/env python3
"""
Vegetation Index Computation from Satellite Data

This module provides a flexible framework for computing multiple vegetation indices
from satellite imagery using STAC APIs. It supports:
- NDVI (Normalized Difference Vegetation Index)
- EVI (Enhanced Vegetation Index)
- SAVI (Soil Adjusted Vegetation Index)
- NDRE (Normalized Difference Red Edge)
- GNDVI (Green Normalized Difference Vegetation Index)

Usage:
    from vegetation_indices import VegetationIndexComputer
    
    computer = VegetationIndexComputer(collection="sentinel-2-l2a")
    indices = computer.compute_indices(item, aoi, ['ndvi', 'evi', 'savi'])
"""

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import warnings

import numpy as np
import rasterio
from rasterio.mask import mask
from rasterio.warp import transform_geom
import matplotlib.pyplot as plt
from pystac_client import Client
import planetary_computer
from shapely.geometry import shape, box


class VegetationIndexComputer:
    """
    Compute various vegetation indices from satellite imagery.
    
    This class handles satellite data retrieval and computation of multiple
    vegetation indices optimized for different agricultural scenarios.
    """
    
    # STAC Catalog URLs
    PLANETARY_COMPUTER_CATALOG = "https://planetarycomputer.microsoft.com/api/stac/v1"
    EARTH_SEARCH_CATALOG = "https://earth-search.aws.element84.com/v1"
    
    # Vegetation Index Definitions
    INDICES = {
        'ndvi': {
            'name': 'Normalized Difference Vegetation Index',
            'bands': ['red', 'nir'],
            'formula': '(NIR - Red) / (NIR + Red)',
            'range': (-1, 1),
            'description': 'General vegetation health and vigor',
            'best_for': 'General purpose vegetation monitoring'
        },
        'evi': {
            'name': 'Enhanced Vegetation Index',
            'bands': ['blue', 'red', 'nir'],
            'formula': '2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1)',
            'range': (-1, 1),
            'description': 'Improved sensitivity in high biomass regions',
            'best_for': 'Dense vegetation, reduced atmospheric effects'
        },
        'savi': {
            'name': 'Soil Adjusted Vegetation Index',
            'bands': ['red', 'nir'],
            'formula': '((NIR - Red) / (NIR + Red + L)) * (1 + L), L=0.5',
            'range': (-1, 1),
            'description': 'Minimizes soil brightness effects',
            'best_for': 'Sparse vegetation, early season crops, exposed soil'
        },
        'ndre': {
            'name': 'Normalized Difference Red Edge',
            'bands': ['rededge', 'nir'],
            'formula': '(NIR - RedEdge) / (NIR + RedEdge)',
            'range': (-1, 1),
            'description': 'Sensitive to chlorophyll content',
            'best_for': 'Precision nitrogen management, disease detection',
            'note': 'Sentinel-2 only (requires Red Edge bands)'
        },
        'gndvi': {
            'name': 'Green Normalized Difference Vegetation Index',
            'bands': ['green', 'nir'],
            'formula': '(NIR - Green) / (NIR + Green)',
            'range': (-1, 1),
            'description': 'Sensitive to chlorophyll concentration',
            'best_for': 'Nitrogen status, photosynthetic activity'
        }
    }
    
    # Band mapping for different satellites
    BAND_MAPPING = {
        'sentinel-2-l2a': {
            'blue': 'B02',
            'green': 'B03',
            'red': 'B04',
            'rededge': 'B05',  # Can also use B06 or B07
            'nir': 'B08',
            'swir1': 'B11',
            'swir2': 'B12'
        },
        'landsat-c2-l2': {
            'blue': 'blue',
            'green': 'green',
            'red': 'red',
            'nir': 'nir08',
            'swir1': 'swir16',
            'swir2': 'swir22'
        }
    }
    
    def __init__(self, catalog_url: str = None, collection: str = "sentinel-2-l2a"):
        """
        Initialize the Vegetation Index Computer.
        
        Args:
            catalog_url: URL of the STAC catalog. Defaults to Planetary Computer.
            collection: Satellite collection to use. Options:
                - 'sentinel-2-l2a' (default, 10m resolution, has Red Edge)
                - 'landsat-c2-l2' (30m resolution, no Red Edge)
        """
        self.catalog_url = catalog_url or self.PLANETARY_COMPUTER_CATALOG
        self.collection = collection
        self.catalog = Client.open(self.catalog_url)
        
        # Check which indices are available for this collection
        self._validate_collection()
        
    def _validate_collection(self):
        """Validate that the collection supports the required bands."""
        if self.collection not in self.BAND_MAPPING:
            available = ', '.join(self.BAND_MAPPING.keys())
            raise ValueError(
                f"Collection '{self.collection}' not supported. "
                f"Available collections: {available}"
            )
        
        # Warn if Red Edge is not available (Landsat)
        if 'rededge' not in self.BAND_MAPPING[self.collection]:
            warnings.warn(
                f"Collection '{self.collection}' does not support Red Edge bands. "
                f"NDRE index will not be available.",
                UserWarning
            )
    
    def get_available_indices(self) -> List[str]:
        """
        Get list of indices available for the current collection.
        
        Returns:
            List of index names (e.g., ['ndvi', 'evi', 'savi'])
        """
        available = []
        band_map = self.BAND_MAPPING[self.collection]
        
        for idx_name, idx_info in self.INDICES.items():
            # Check if all required bands are available
            bands_available = all(
                band in band_map for band in idx_info['bands']
            )
            if bands_available:
                available.append(idx_name)
        
        return available
    
    def print_index_info(self, index_name: str = None):
        """
        Print information about vegetation indices.
        
        Args:
            index_name: Specific index to show info for, or None for all
        """
        if index_name:
            if index_name not in self.INDICES:
                print(f"Unknown index: {index_name}")
                return
            indices = {index_name: self.INDICES[index_name]}
        else:
            indices = self.INDICES
        
        available = self.get_available_indices()
        
        print(f"\n{'='*70}")
        print(f"VEGETATION INDICES - Collection: {self.collection}")
        print(f"{'='*70}\n")
        
        for name, info in indices.items():
            status = "✅ Available" if name in available else "❌ Not available"
            print(f"{name.upper()} - {info['name']}")
            print(f"Status: {status}")
            print(f"Formula: {info['formula']}")
            print(f"Range: {info['range']}")
            print(f"Best for: {info['best_for']}")
            if 'note' in info:
                print(f"Note: {info['note']}")
            print()
    
    def load_aoi(self, aoi_input: Union[str, Path, Dict, List]) -> Dict:
        """
        Load Area of Interest from various input formats.
        
        Args:
            aoi_input: Can be:
                - Path to GeoJSON file
                - GeoJSON dict
                - Bounding box as [min_lon, min_lat, max_lon, max_lat]
                - List of coordinate tuples [(lon1, lat1), (lon2, lat2), ...]
                
        Returns:
            GeoJSON geometry dict
        """
        if isinstance(aoi_input, (str, Path)):
            with open(aoi_input, 'r') as f:
                geojson = json.load(f)
            
            if geojson.get('type') == 'Feature':
                return geojson['geometry']
            elif geojson.get('type') == 'FeatureCollection':
                return geojson['features'][0]['geometry']
            else:
                return geojson
                
        elif isinstance(aoi_input, dict):
            return aoi_input
            
        elif isinstance(aoi_input, list):
            if len(aoi_input) == 4 and all(isinstance(x, (int, float)) for x in aoi_input):
                min_lon, min_lat, max_lon, max_lat = aoi_input
                return {
                    "type": "Polygon",
                    "coordinates": [[
                        [min_lon, min_lat],
                        [max_lon, min_lat],
                        [max_lon, max_lat],
                        [min_lon, max_lat],
                        [min_lon, min_lat]
                    ]]
                }
            elif len(aoi_input) >= 3 and all(
                isinstance(coord, (tuple, list)) and len(coord) >= 2 
                for coord in aoi_input
            ):
                coords = [[float(coord[0]), float(coord[1])] for coord in aoi_input]
                if coords[0] != coords[-1]:
                    coords.append(coords[0])
                return {
                    "type": "Polygon",
                    "coordinates": [coords]
                }
            else:
                raise ValueError(
                    "List input must be either a bounding box [min_lon, min_lat, max_lon, max_lat] "
                    "or a list of coordinate tuples [(lon1, lat1), (lon2, lat2), ...]"
                )
        else:
            raise ValueError(
                "AOI input must be a file path, GeoJSON dict, bounding box list, "
                "or list of coordinate tuples"
            )
    
    def search_satellite_data(
        self,
        aoi_geometry: Dict,
        start_date: str,
        end_date: str,
        max_cloud_cover: float = 20
    ) -> List:
        """
        Search for satellite imagery covering the AOI.
        
        Args:
            aoi_geometry: GeoJSON geometry dict
            start_date: Start date in YYYY-MM-DD format
            end_date: End date in YYYY-MM-DD format
            max_cloud_cover: Maximum cloud cover percentage (0-100)
            
        Returns:
            List of STAC items
        """
        search = self.catalog.search(
            collections=[self.collection],
            intersects=aoi_geometry,
            datetime=f"{start_date}/{end_date}",
            query={"eo:cloud_cover": {"lt": max_cloud_cover}}
        )
        
        items = list(search.items())
        print(f"Found {len(items)} scenes matching criteria")
        
        return items
    
    def _load_band(
        self,
        item,
        band_name: str,
        aoi_geometry: Dict
    ) -> Tuple[np.ndarray, Dict, object]:
        """
        Load and clip a single band.
        
        Args:
            item: STAC item
            band_name: Logical band name (e.g., 'red', 'nir')
            aoi_geometry: GeoJSON geometry to clip to
            
        Returns:
            Tuple of (clipped_array, metadata, transform)
        """
        # Get the actual band key for this collection
        band_map = self.BAND_MAPPING[self.collection]
        if band_name not in band_map:
            raise ValueError(
                f"Band '{band_name}' not available in collection '{self.collection}'"
            )
        
        band_key = band_map[band_name]
        
        # Sign assets if using Planetary Computer
        if 'planetarycomputer' in self.catalog_url:
            item = planetary_computer.sign(item)
        
        # Get band URL
        band_href = item.assets[band_key].href
        
        # Load and clip
        with rasterio.open(band_href) as src:
            # Transform AOI to match raster CRS
            aoi_transformed = aoi_geometry
            if aoi_geometry.get('crs') != src.crs or 'crs' not in aoi_geometry:
                aoi_transformed = transform_geom(
                    'EPSG:4326',
                    src.crs,
                    aoi_geometry
                )
            
            # Check intersection
            raster_bounds = box(*src.bounds)
            aoi_shapely = shape(aoi_transformed)
            
            if not raster_bounds.intersects(aoi_shapely):
                raise ValueError(
                    f"AOI does not intersect with this scene's coverage area."
                )
            
            # Clip to AOI
            aoi_shape = [aoi_transformed]
            band_clipped, band_transform = mask(src, aoi_shape, crop=True)
            band_clipped = band_clipped[0].astype(float)
            
            # Get metadata
            meta = src.meta.copy()
            meta.update({
                "driver": "GTiff",
                "height": band_clipped.shape[0],
                "width": band_clipped.shape[1],
                "transform": band_transform,
                "count": 1,
                "dtype": 'float32'
            })
        
        return band_clipped, meta, band_transform
    
    def _compute_index(
        self,
        bands: Dict[str, np.ndarray],
        index_name: str
    ) -> np.ndarray:
        """
        Compute a vegetation index from loaded bands.
        
        Args:
            bands: Dictionary of band arrays (e.g., {'red': array, 'nir': array})
            index_name: Name of index to compute
            
        Returns:
            Computed index array
        """
        index_name = index_name.lower()
        
        if index_name == 'ndvi':
            nir = bands['nir']
            red = bands['red']
            denominator = nir + red
            return np.where(
                denominator != 0,
                (nir - red) / denominator,
                np.nan
            )
        
        elif index_name == 'evi':
            nir = bands['nir']
            red = bands['red']
            blue = bands['blue']
            denominator = nir + 6 * red - 7.5 * blue + 1
            return np.where(
                denominator != 0,
                2.5 * (nir - red) / denominator,
                np.nan
            )
        
        elif index_name == 'savi':
            nir = bands['nir']
            red = bands['red']
            L = 0.5  # Soil brightness correction factor
            denominator = nir + red + L
            return np.where(
                denominator != 0,
                ((nir - red) / denominator) * (1 + L),
                np.nan
            )
        
        elif index_name == 'ndre':
            nir = bands['nir']
            rededge = bands['rededge']
            denominator = nir + rededge
            return np.where(
                denominator != 0,
                (nir - rededge) / denominator,
                np.nan
            )
        
        elif index_name == 'gndvi':
            nir = bands['nir']
            green = bands['green']
            denominator = nir + green
            return np.where(
                denominator != 0,
                (nir - green) / denominator,
                np.nan
            )
        
        else:
            raise ValueError(f"Unknown index: {index_name}")
    
    def compute_indices(
        self,
        item,
        aoi_geometry: Dict,
        indices: List[str] = ['ndvi'],
        output_dir: Optional[str] = None
    ) -> Dict[str, Tuple[np.ndarray, Dict]]:
        """
        Compute multiple vegetation indices efficiently.
        
        This method loads each required band only once, even when computing
        multiple indices.
        
        Args:
            item: STAC item
            aoi_geometry: GeoJSON geometry to clip to
            indices: List of index names to compute (e.g., ['ndvi', 'evi'])
            output_dir: Optional directory to save rasters
            
        Returns:
            Dictionary mapping index names to (array, metadata) tuples
        """
        # Validate indices
        available = self.get_available_indices()
        for idx in indices:
            if idx not in available:
                raise ValueError(
                    f"Index '{idx}' not available for collection '{self.collection}'. "
                    f"Available indices: {available}"
                )
        
        print(f"\nProcessing scene: {item.id}")
        print(f"Date: {item.datetime}")
        print(f"Computing indices: {', '.join(indices)}")
        
        # Determine which bands we need to load
        required_bands = set()
        for idx in indices:
            required_bands.update(self.INDICES[idx]['bands'])
        
        print(f"Loading bands: {', '.join(sorted(required_bands))}")
        
        # Load all required bands (only once each)
        bands = {}
        band_metadata = {}
        metadata = None
        
        for band_name in required_bands:
            try:
                band_array, meta, transform = self._load_band(item, band_name, aoi_geometry)
                bands[band_name] = band_array
                band_metadata[band_name] = (meta, transform)
                if metadata is None:
                    metadata = meta
            except Exception as e:
                raise ValueError(f"Error loading {band_name} band: {str(e)}")
        
        # Resample bands to match the highest resolution (largest array)
        # This is needed because Sentinel-2 bands have different resolutions
        # (e.g., Red Edge at 20m vs NIR at 10m)
        target_shape = max(bands.values(), key=lambda x: x.size).shape
        
        for band_name, band_array in bands.items():
            if band_array.shape != target_shape:
                print(f"    Resampling {band_name} from {band_array.shape} to {target_shape}")
                # Use scipy zoom for resampling
                from scipy.ndimage import zoom
                zoom_factors = (target_shape[0] / band_array.shape[0],
                              target_shape[1] / band_array.shape[1])
                bands[band_name] = zoom(band_array, zoom_factors, order=1)
        
        # Update metadata to reflect the target shape
        metadata.update({
            "height": target_shape[0],
            "width": target_shape[1]
        })
        
        # Compute each index
        results = {}
        
        for index_name in indices:
            print(f"  Computing {index_name.upper()}...")
            
            index_array = self._compute_index(bands, index_name)
            
            # Clip to valid range
            valid_range = self.INDICES[index_name]['range']
            index_array = np.clip(index_array, valid_range[0], valid_range[1])
            
            # Save if output directory provided
            if output_dir:
                output_path = Path(output_dir)
                output_path.mkdir(exist_ok=True, parents=True)
                
                output_file = output_path / f"{index_name}_{item.id}.tif"
                with rasterio.open(output_file, 'w', **metadata) as dst:
                    dst.write(index_array.astype('float32'), 1)
                print(f"    Saved to: {output_file}")
            
            results[index_name] = (index_array, metadata)
        
        print("✅ All indices computed successfully!\n")
        
        return results
    
    def visualize_indices(
        self,
        results: Dict[str, Tuple[np.ndarray, Dict]],
        scene_id: str = "",
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = None
    ):
        """
        Create side-by-side visualization of multiple indices.
        
        Args:
            results: Dictionary from compute_indices()
            scene_id: Scene identifier for title
            save_path: Optional path to save figure
            figsize: Figure size (width, height). Auto-calculated if None.
        """
        n_indices = len(results)
        
        if n_indices == 0:
            print("No indices to visualize")
            return
        
        # Calculate layout
        if n_indices <= 3:
            nrows, ncols = 1, n_indices
        else:
            ncols = 3
            nrows = (n_indices + ncols - 1) // ncols
        
        # Calculate figure size
        if figsize is None:
            figsize = (6 * ncols, 5 * nrows)
        
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        
        # Handle single subplot case
        if n_indices == 1:
            axes = [axes]
        elif nrows > 1:
            axes = axes.flatten()
        
        for idx, (index_name, (index_array, _)) in enumerate(results.items()):
            ax = axes[idx] if n_indices > 1 else axes[0]
            
            # Plot
            im = ax.imshow(index_array, cmap='RdYlGn', vmin=-1, vmax=1)
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=ax, shrink=0.8)
            cbar.set_label(index_name.upper(), rotation=270, labelpad=15)
            
            # Title and labels
            info = self.INDICES[index_name]
            ax.set_title(f"{index_name.upper()} - {info['name']}", 
                        fontsize=11, fontweight='bold')
            ax.set_xlabel('Column (pixels)')
            ax.set_ylabel('Row (pixels)')
            
            # Add statistics
            valid_data = index_array[~np.isnan(index_array)]
            if len(valid_data) > 0:
                stats_text = (
                    f"Mean: {valid_data.mean():.3f}\n"
                    f"Std: {valid_data.std():.3f}\n"
                    f"Range: [{valid_data.min():.3f}, {valid_data.max():.3f}]"
                )
                ax.text(0.02, 0.98, stats_text,
                       transform=ax.transAxes,
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                       fontsize=8)
        
        # Hide extra subplots
        if nrows > 1:
            for idx in range(n_indices, len(axes)):
                axes[idx].axis('off')
        
        # Main title
        if scene_id:
            fig.suptitle(f'Vegetation Indices Comparison - {scene_id}', 
                        fontsize=14, fontweight='bold', y=0.995)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")
        
        plt.show()
    
    def compute_statistics(
        self,
        results: Dict[str, Tuple[np.ndarray, Dict]]
    ) -> Dict[str, Dict]:
        """
        Compute statistics for all indices.
        
        Args:
            results: Dictionary from compute_indices()
            
        Returns:
            Dictionary mapping index names to statistics dicts
        """
        stats = {}
        
        for index_name, (index_array, _) in results.items():
            valid_data = index_array[~np.isnan(index_array)]
            
            if len(valid_data) > 0:
                stats[index_name] = {
                    'mean': float(valid_data.mean()),
                    'median': float(np.median(valid_data)),
                    'std': float(valid_data.std()),
                    'min': float(valid_data.min()),
                    'max': float(valid_data.max()),
                    'valid_pixels': len(valid_data),
                    'total_pixels': index_array.size
                }
            else:
                stats[index_name] = {
                    'error': 'No valid data'
                }
        
        return stats


def compute_vegetation_indices_for_aoi(
    aoi_input: Union[str, Path, Dict, List],
    start_date: str,
    end_date: str,
    indices: List[str] = ['ndvi', 'evi', 'savi'],
    output_dir: str = "output",
    max_cloud_cover: float = 20,
    collection: str = "sentinel-2-l2a",
    visualize: bool = True
) -> List[Tuple[str, Dict]]:
    """
    High-level function to compute multiple vegetation indices for an AOI.
    
    Args:
        aoi_input: AOI as file path, GeoJSON dict, or bounding box
        start_date: Start date (YYYY-MM-DD)
        end_date: End date (YYYY-MM-DD)
        indices: List of index names to compute
        output_dir: Directory to save outputs
        max_cloud_cover: Maximum cloud cover percentage
        collection: Satellite collection to use
        visualize: Whether to create visualizations
        
    Returns:
        List of tuples (scene_id, results_dict)
    """
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Initialize computer
    computer = VegetationIndexComputer(collection=collection)
    
    # Print available indices
    print(f"\n{'='*70}")
    print(f"Vegetation Index Computer - {collection}")
    print(f"{'='*70}")
    print(f"Requested indices: {', '.join(indices)}")
    print(f"Available indices: {', '.join(computer.get_available_indices())}")
    print(f"{'='*70}\n")
    
    # Load AOI
    print("Loading AOI...")
    aoi_geometry = computer.load_aoi(aoi_input)
    
    # Search for data
    print(f"\nSearching for {collection} data from {start_date} to {end_date}...")
    items = computer.search_satellite_data(
        aoi_geometry,
        start_date,
        end_date,
        max_cloud_cover
    )
    
    if not items:
        print("No scenes found matching criteria!")
        return []
    
    # Process each scene
    all_results = []
    skipped = 0
    
    for i, item in enumerate(items):
        scene_id = item.id
        
        print(f"\n{'='*70}")
        print(f"Processing scene {i+1}/{len(items)}")
        print(f"{'='*70}")
        
        try:
            # Compute indices
            results = computer.compute_indices(
                item,
                aoi_geometry,
                indices=indices,
                output_dir=str(output_path)
            )
            
            all_results.append((scene_id, results))
            
            # Visualize
            if visualize:
                fig_path = output_path / f"comparison_{scene_id}.png"
                computer.visualize_indices(
                    results,
                    scene_id=f"{scene_id} ({item.datetime.strftime('%Y-%m-%d')})",
                    save_path=str(fig_path)
                )
            
            # Print statistics
            stats = computer.compute_statistics(results)
            print("\n📊 Statistics Summary:")
            for idx_name, idx_stats in stats.items():
                if 'error' not in idx_stats:
                    print(f"\n  {idx_name.upper()}:")
                    print(f"    Mean: {idx_stats['mean']:.3f}")
                    print(f"    Range: [{idx_stats['min']:.3f}, {idx_stats['max']:.3f}]")
                    print(f"    Valid pixels: {idx_stats['valid_pixels']:,}")
        
        except ValueError as e:
            if "does not intersect" in str(e):
                print(f"⚠️  Scene doesn't cover AOI, skipping...")
                skipped += 1
                continue
            else:
                raise
        except Exception as e:
            print(f"\n⚠️  Error processing scene: {str(e)}")
            skipped += 1
            continue
    
    print(f"\n{'='*70}")
    print(f"PROCESSING COMPLETE")
    print(f"{'='*70}")
    print(f"✅ Successfully processed: {len(all_results)} scenes")
    if skipped > 0:
        print(f"⚠️  Skipped: {skipped} scenes")
    print(f"📁 Outputs saved to: {output_path.absolute()}")
    print(f"{'='*70}\n")
    
    return all_results


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description="Compute vegetation indices from satellite data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compute NDVI, EVI, and SAVI for a field
  vegetation-index-compute --aoi field.geojson \\
      --start-date 2024-06-01 --end-date 2024-06-30 \\
      --indices ndvi evi savi
  
  # Compute all available indices using bounding box
  vegetation-index-compute --aoi "-121.0,37.0,-120.5,37.5" \\
      --start-date 2024-06-01 --end-date 2024-06-30 \\
      --indices ndvi evi savi ndre gndvi
  
  # List available indices
  vegetation-index-compute --list-indices
        """
    )
    parser.add_argument(
        '--aoi',
        help='Area of Interest: path to GeoJSON file or bounding box as "min_lon,min_lat,max_lon,max_lat"'
    )
    parser.add_argument(
        '--start-date',
        help='Start date in YYYY-MM-DD format'
    )
    parser.add_argument(
        '--end-date',
        help='End date in YYYY-MM-DD format'
    )
    parser.add_argument(
        '--indices',
        nargs='+',
        default=['ndvi'],
        help='Vegetation indices to compute (e.g., ndvi evi savi ndre gndvi)'
    )
    parser.add_argument(
        '--output-dir',
        default='output',
        help='Output directory for rasters and figures (default: output)'
    )
    parser.add_argument(
        '--max-cloud-cover',
        type=float,
        default=20,
        help='Maximum cloud cover percentage (default: 20)'
    )
    parser.add_argument(
        '--collection',
        default='sentinel-2-l2a',
        choices=['sentinel-2-l2a', 'landsat-c2-l2'],
        help='Satellite collection to use (default: sentinel-2-l2a)'
    )
    parser.add_argument(
        '--no-visualize',
        action='store_true',
        help='Skip visualization'
    )
    parser.add_argument(
        '--list-indices',
        action='store_true',
        help='List available indices and exit'
    )
    
    args = parser.parse_args()
    
    # Handle list indices
    if args.list_indices:
        computer = VegetationIndexComputer(collection=args.collection)
        computer.print_index_info()
        return
    
    # Validate required arguments
    if not args.aoi or not args.start_date or not args.end_date:
        parser.error("--aoi, --start-date, and --end-date are required (or use --list-indices)")
    
    # Parse AOI
    if args.aoi.endswith('.geojson') or args.aoi.endswith('.json'):
        aoi_input = args.aoi
    else:
        try:
            bbox = [float(x.strip()) for x in args.aoi.split(',')]
            if len(bbox) != 4:
                raise ValueError
            aoi_input = bbox
        except:
            parser.error("AOI must be a GeoJSON file path or bounding box as 'min_lon,min_lat,max_lon,max_lat'")
    
    # Run computation
    compute_vegetation_indices_for_aoi(
        aoi_input=aoi_input,
        start_date=args.start_date,
        end_date=args.end_date,
        indices=args.indices,
        output_dir=args.output_dir,
        max_cloud_cover=args.max_cloud_cover,
        collection=args.collection,
        visualize=not args.no_visualize
    )


if __name__ == "__main__":
    main()

