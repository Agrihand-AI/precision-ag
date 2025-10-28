#!/usr/bin/env python3
"""
Management Zone Analysis for Precision Agriculture

This module provides tools for creating management zones based on satellite-derived
indices. Inspired by precision agriculture research, particularly the work of
Ali Mirzakhani Nafchi (SDSU Extension) on soil variability and management zones.

Management zones are areas within a field that exhibit similar characteristics
and can be managed uniformly. This approach enables:
- Variable rate application of inputs
- Targeted soil sampling
- Precision irrigation management
- Optimized resource allocation

Usage:
    from precision_ag.zone_analysis import ManagementZoneAnalyzer
    
    analyzer = ManagementZoneAnalyzer()
    zones = analyzer.create_zones(indices_data, n_zones=3)
    analyzer.visualize_zones(zones, field_data)
"""

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import rasterio
from rasterio.transform import from_bounds
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA


class ManagementZoneAnalyzer:
    """
    Create and analyze management zones for precision agriculture.
    
    This class implements clustering-based zone delineation using
    multiple satellite-derived indices. The approach complements
    ground-based methods (e.g., EC mapping) by providing:
    - Broad spatial coverage
    - Temporal analysis capabilities
    - Cost-effective preliminary assessment
    """
    
    def __init__(self):
        """Initialize the Management Zone Analyzer."""
        self.scaler = StandardScaler()
        self.zone_colors = ['#d73027', '#fc8d59', '#fee090', '#91bfdb', '#4575b4']
        
    def prepare_data(
        self,
        indices_results: Dict[str, Tuple[np.ndarray, Dict]],
        mask_invalid: bool = True
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Prepare index data for clustering analysis.
        
        Args:
            indices_results: Dictionary from compute_indices() containing
                           (index_array, metadata) tuples
            mask_invalid: Whether to mask NaN/invalid values
            
        Returns:
            Tuple of (stacked_data, valid_mask) where:
                - stacked_data: (n_pixels, n_indices) array
                - valid_mask: boolean array of shape matching original rasters
        """
        index_names = list(indices_results.keys())
        arrays = []
        shapes = []
        
        print(f"\nPreparing data for {len(index_names)} indices...")
        
        for idx_name in index_names:
            index_array, metadata = indices_results[idx_name]
            shapes.append(index_array.shape)
            
            # Flatten for clustering
            flat = index_array.flatten()
            arrays.append(flat)
            
            # Count valid pixels
            valid_count = np.sum(~np.isnan(flat))
            print(f"  {idx_name.upper()}: {valid_count:,} valid pixels")
        
        # Verify all arrays have same shape
        if len(set(shapes)) > 1:
            raise ValueError(f"Index arrays have different shapes: {shapes}")
        
        # Stack indices (each column is one index)
        stacked = np.column_stack(arrays)
        
        # Create mask for valid data (all indices must be valid)
        if mask_invalid:
            valid_mask = ~np.any(np.isnan(stacked), axis=1)
            print(f"\nValid pixels (all indices): {np.sum(valid_mask):,} / {len(valid_mask):,}")
        else:
            valid_mask = np.ones(stacked.shape[0], dtype=bool)
        
        return stacked, valid_mask.reshape(shapes[0])
    
    def tune_zones(
        self,
        indices_results: Dict[str, Tuple[np.ndarray, Dict]],
        min_clusters: int = 2,
        max_clusters: int = 10,
        n_init: int = 10,
        max_iter: int = 300
    ) -> Dict:
        """
        Tunes the number of clusters for K-means clustering to find the optimal
        number of management zones based on silhouette score.
        
        Args:
            indices_results: Dictionary from compute_indices()
            min_clusters: Minimum number of clusters to try
            max_clusters: Maximum number of clusters to try
            n_init: Number of initializations for K-means
            max_iter: Maximum number of iterations for K-means
            
        Returns:
            Dictionary containing:
                - 'optimal_n_clusters': The number of clusters with the highest silhouette score
                - 'silhouette_scores': Silhouette scores for each number of clusters
        """
        print(f"\n{'='*70}")
        print("TUNING ZONE CLUSTERS")
        print(f"{'='*70}")
        
        stacked_data, valid_mask = self.prepare_data(indices_results)
        valid_data = stacked_data[valid_mask.flatten()]
        
        # Standardize features (important for clustering)
        scaled_data = self.scaler.fit_transform(valid_data)
        
        silhouette_scores = {}
        for n_clusters in range(min_clusters, max_clusters + 1):
            print(f"  Testing {n_clusters} clusters...")
            kmeans = KMeans(
                n_clusters=n_clusters,
                random_state=42,
                n_init=n_init,
                max_iter=max_iter
            )
            labels = kmeans.fit_predict(scaled_data)
            
            # Calculate silhouette score
            if len(set(labels)) > 1: # Silhouette score is not defined for a single cluster
                silhouette_scores[n_clusters] = silhouette_score(scaled_data, labels)
            else:
                silhouette_scores[n_clusters] = 0.0 # Assign a low score for single cluster
            
            print(f"    Silhouette score for {n_clusters} clusters: {silhouette_scores[n_clusters]:.3f}")
        
        optimal_n_clusters = max(silhouette_scores, key=silhouette_scores.get)
        print(f"\nOptimal number of clusters: {optimal_n_clusters} (Silhouette score: {silhouette_scores[optimal_n_clusters]:.3f})")
        
        return {
            'optimal_n_clusters': optimal_n_clusters,
            'silhouette_scores': silhouette_scores
        }
    
    def create_zones(
        self,
        indices_results: Dict[str, Tuple[np.ndarray, Dict]],
        n_zones: int = 3,
        method: str = 'kmeans',
        random_state: int = 42
    ) -> Dict:
        """
        Create management zones using clustering analysis.
        
        This method groups pixels with similar spectral characteristics,
        creating zones that can be managed uniformly. The approach is
        similar to EC-based zone delineation but uses satellite data.
        
        Args:
            indices_results: Dictionary from compute_indices()
            n_zones: Number of management zones to create (typically 3-5)
            method: Clustering method ('kmeans', 'hierarchical')
            random_state: Random seed for reproducibility
            
        Returns:
            Dictionary containing:
                - 'zones': 2D array with zone labels (0 to n_zones-1, -1 for invalid)
                - 'centers': Cluster centers in original data space
                - 'metadata': Original metadata from indices
                - 'index_names': List of indices used
                - 'statistics': Statistics for each zone
        """
        if method not in ['kmeans']:
            raise ValueError(f"Method '{method}' not supported. Use 'kmeans'.")
        
        print(f"\n{'='*70}")
        print(f"CREATING {n_zones} MANAGEMENT ZONES")
        print(f"{'='*70}")
        
        # Prepare data
        stacked_data, valid_mask = self.prepare_data(indices_results)
        original_shape = valid_mask.shape
        
        # Extract valid data for clustering
        valid_data = stacked_data[valid_mask.flatten()]
        
        # Standardize features (important for clustering)
        print(f"\nStandardizing features...")
        scaled_data = self.scaler.fit_transform(valid_data)
        
        # Perform clustering
        print(f"Performing K-means clustering with {n_zones} clusters...")
        kmeans = KMeans(
            n_clusters=n_zones,
            random_state=random_state,
            n_init=10,
            max_iter=300
        )
        labels = kmeans.fit_predict(scaled_data)
        
        # Create zone array (full size, with -1 for invalid pixels)
        zones = np.full(valid_mask.size, -1, dtype=np.int32)
        zones[valid_mask.flatten()] = labels
        zones = zones.reshape(original_shape)
        
        # Transform cluster centers back to original scale
        centers_scaled = kmeans.cluster_centers_
        centers_original = self.scaler.inverse_transform(centers_scaled)
        
        # Calculate statistics for each zone
        statistics = self._calculate_zone_statistics(
            zones, indices_results, n_zones
        )
        
        # Get metadata from first index
        first_index = list(indices_results.values())[0]
        metadata = first_index[1].copy()
        
        result = {
            'zones': zones,
            'centers': centers_original,
            'metadata': metadata,
            'index_names': list(indices_results.keys()),
            'n_zones': n_zones,
            'statistics': statistics,
            'valid_mask': valid_mask
        }
        
        print(f"\n✅ Zone creation complete!")
        self._print_zone_summary(result)
        
        return result
    
    def _calculate_zone_statistics(
        self,
        zones: np.ndarray,
        indices_results: Dict[str, Tuple[np.ndarray, Dict]],
        n_zones: int
    ) -> Dict:
        """Calculate statistics for each management zone."""
        stats = {}
        
        for zone_id in range(n_zones):
            zone_mask = zones == zone_id
            zone_pixels = np.sum(zone_mask)
            
            zone_stats = {
                'pixel_count': int(zone_pixels),
                'area_percentage': float(zone_pixels / np.sum(zones >= 0) * 100),
                'indices': {}
            }
            
            # Calculate statistics for each index
            for idx_name, (index_array, _) in indices_results.items():
                zone_values = index_array[zone_mask]
                valid_values = zone_values[~np.isnan(zone_values)]
                
                if len(valid_values) > 0:
                    zone_stats['indices'][idx_name] = {
                        'mean': float(np.mean(valid_values)),
                        'std': float(np.std(valid_values)),
                        'min': float(np.min(valid_values)),
                        'max': float(np.max(valid_values))
                    }
            
            stats[f'zone_{zone_id}'] = zone_stats
        
        return stats
    
    def _print_zone_summary(self, result: Dict):
        """Print summary of created zones."""
        print(f"\n{'='*70}")
        print("ZONE SUMMARY")
        print(f"{'='*70}\n")
        
        stats = result['statistics']
        index_names = result['index_names']
        
        for zone_key in sorted(stats.keys()):
            zone_num = int(zone_key.split('_')[1])
            zone_data = stats[zone_key]
            
            print(f"Zone {zone_num}:")
            print(f"  Pixels: {zone_data['pixel_count']:,} ({zone_data['area_percentage']:.1f}%)")
            print(f"  Characteristics:")
            
            for idx_name in index_names:
                if idx_name in zone_data['indices']:
                    idx_stats = zone_data['indices'][idx_name]
                    print(f"    {idx_name.upper()}: {idx_stats['mean']:.3f} "
                          f"(±{idx_stats['std']:.3f})")
            print()
    
    def visualize_zones(
        self,
        result: Dict,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (14, 10),
        show_indices: bool = True
    ):
        """
        Visualize management zones and underlying indices.
        
        Args:
            result: Dictionary from create_zones()
            save_path: Optional path to save figure
            figsize: Figure size (width, height)
            show_indices: Whether to show individual indices alongside zones
        """
        zones = result['zones']
        n_zones = result['n_zones']
        index_names = result['index_names']
        
        # Determine layout
        if show_indices and len(index_names) > 0:
            n_plots = 1 + len(index_names)
            ncols = min(3, n_plots)
            nrows = (n_plots + ncols - 1) // ncols
        else:
            nrows, ncols = 1, 1
        
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        if n_plots == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        
        # Plot zones
        ax = axes[0]
        
        # Create colormap for zones
        n_colors = min(n_zones, len(self.zone_colors))
        zone_cmap = ListedColormap(self.zone_colors[:n_colors])
        
        # Plot with invalid areas shown in gray
        zones_display = np.ma.masked_where(zones < 0, zones)
        im = ax.imshow(zones_display, cmap=zone_cmap, vmin=0, vmax=n_zones-1)
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Management Zone', rotation=270, labelpad=20)
        cbar.set_ticks(range(n_zones))
        cbar.set_ticklabels([f'Zone {i}' for i in range(n_zones)])
        
        ax.set_title('Management Zones', fontsize=14, fontweight='bold')
        ax.set_xlabel('Column (pixels)')
        ax.set_ylabel('Row (pixels)')
        
        # Add area statistics
        stats_text = []
        for zone_id in range(n_zones):
            zone_mask = zones == zone_id
            pct = np.sum(zone_mask) / np.sum(zones >= 0) * 100
            stats_text.append(f"Zone {zone_id}: {pct:.1f}%")
        
        ax.text(0.02, 0.98, '\n'.join(stats_text),
               transform=ax.transAxes,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
               fontsize=9)
        
        # Plot individual indices if requested
        if show_indices and 'indices_data' in result:
            for idx, idx_name in enumerate(index_names, start=1):
                if idx < len(axes):
                    ax = axes[idx]
                    index_array = result['indices_data'][idx_name]
                    
                    # Handle different index ranges
                    if idx_name == 'si':
                        # Normalize SI to [0, 1] using min-max scaling
                        valid_data = index_array[~np.isnan(index_array)]
                        if len(valid_data) > 0:
                            data_min = valid_data.min()
                            data_max = valid_data.max()
                            if data_max > data_min:
                                display_array = (index_array - data_min) / (data_max - data_min)
                            else:
                                display_array = np.zeros_like(index_array)
                            vmin, vmax = 0, 1
                        else:
                            display_array = index_array
                            vmin, vmax = 0, 1
                    else:
                        # Use [-1, 1] for other indices
                        display_array = index_array
                        vmin, vmax = -1, 1
                    
                    im = ax.imshow(display_array, cmap='RdYlGn', vmin=vmin, vmax=vmax)
                    plt.colorbar(im, ax=ax, shrink=0.8)
                    
                    ax.set_title(f'{idx_name.upper()}', fontsize=12, fontweight='bold')
                    ax.set_xlabel('Column (pixels)')
                    ax.set_ylabel('Row (pixels)')
        
        # Hide unused subplots
        for idx in range(n_plots, len(axes)):
            axes[idx].axis('off')
        
        fig.suptitle('Management Zone Analysis', fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"\nFigure saved to: {save_path}")
        
        plt.show()
    
    def export_zones(
        self,
        result: Dict,
        output_path: Union[str, Path],
        zone_interpretation: Optional[Dict[int, str]] = None
    ):
        """
        Export management zones as a GeoTIFF file.
        
        This creates a georeferenced raster file that can be imported into
        farm management software, GIS systems, or used to create prescription maps.
        
        Args:
            result: Dictionary from create_zones()
            output_path: Path for output GeoTIFF file
            zone_interpretation: Optional dict mapping zone IDs to interpretations
                                (e.g., {0: 'Low productivity', 1: 'Medium', 2: 'High'})
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(exist_ok=True, parents=True)
        
        zones = result['zones']
        metadata = result['metadata'].copy()
        
        # Update metadata for integer zones
        metadata.update({
            'dtype': 'int16',
            'nodata': -1
        })
        
        # Write zones to file
        with rasterio.open(output_path, 'w', **metadata) as dst:
            dst.write(zones.astype('int16'), 1)
            
            # Add descriptions if provided
            if zone_interpretation:
                dst.set_band_description(1, "Management Zones")
        
        print(f"\n✅ Zones exported to: {output_path}")
        
        if zone_interpretation:
            print("\nZone Interpretation:")
            for zone_id, description in zone_interpretation.items():
                print(f"  Zone {zone_id}: {description}")
    
    def recommend_actions(
        self,
        result: Dict,
        index_thresholds: Optional[Dict[str, Dict]] = None
    ) -> Dict[int, List[str]]:
        """
        Generate management recommendations for each zone.
        
        Based on zone characteristics and index values, suggests actions
        such as variable rate applications, soil testing priorities, etc.
        
        Args:
            result: Dictionary from create_zones()
            index_thresholds: Optional custom thresholds for decision making
            
        Returns:
            Dictionary mapping zone IDs to lists of recommendations
        """
        recommendations = {}
        stats = result['statistics']
        
        # Default thresholds (can be customized based on research)
        if index_thresholds is None:
            index_thresholds = {
                'ndvi': {'low': 0.3, 'high': 0.7},
                'ndmi': {'low': 0.2, 'high': 0.5},
                'bsi': {'low': 0.0, 'high': 0.3}
            }
        
        for zone_key in sorted(stats.keys()):
            zone_id = int(zone_key.split('_')[1])
            zone_data = stats[zone_key]
            actions = []
            
            # Analyze each index
            for idx_name, idx_stats in zone_data['indices'].items():
                mean_value = idx_stats['mean']
                
                # NDVI-based recommendations
                if idx_name == 'ndvi':
                    if mean_value < index_thresholds['ndvi']['low']:
                        actions.append("❗ Low vegetation vigor - investigate cause (soil, nutrients, water)")
                        actions.append("  → Consider soil testing and nutrient analysis")
                    elif mean_value > index_thresholds['ndvi']['high']:
                        actions.append("✅ Healthy vegetation - standard management")
                
                # NDMI-based recommendations  
                elif idx_name == 'ndmi':
                    if mean_value < index_thresholds['ndmi']['low']:
                        actions.append("💧 Low moisture - increase irrigation or monitor closely")
                    elif mean_value > index_thresholds['ndmi']['high']:
                        actions.append("💧 Good moisture retention - reduce irrigation if applicable")
                
                # BSI-based recommendations
                elif idx_name == 'bsi':
                    if mean_value > index_thresholds['bsi']['high']:
                        actions.append("🏜️  High bare soil exposure - assess erosion risk")
                        actions.append("  → Consider cover crops or mulching")
            
            if not actions:
                actions.append("✅ Normal conditions - standard management practices apply")
            
            recommendations[zone_id] = actions
        
        # Print recommendations
        print(f"\n{'='*70}")
        print("MANAGEMENT RECOMMENDATIONS")
        print(f"{'='*70}\n")
        
        for zone_id in sorted(recommendations.keys()):
            print(f"Zone {zone_id} ({stats[f'zone_{zone_id}']['area_percentage']:.1f}% of field):")
            for action in recommendations[zone_id]:
                print(f"  {action}")
            print()
        
        return recommendations


def create_management_zones_for_field(
    indices_results: Dict[str, Tuple[np.ndarray, Dict]],
    n_zones: int = 3,
    output_dir: str = "zones",
    visualize: bool = True,
    export: bool = True
) -> Dict:
    """
    High-level function to create management zones from index data.
    
    Args:
        indices_results: Dictionary from compute_indices() or compute_agricultural_indices_for_aoi()
        n_zones: Number of management zones (typically 3-5)
        output_dir: Directory for outputs
        visualize: Whether to create visualization
        export: Whether to export zones as GeoTIFF
        
    Returns:
        Dictionary with zone analysis results
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Create analyzer
    analyzer = ManagementZoneAnalyzer()
    
    # Create zones
    result = analyzer.create_zones(indices_results, n_zones=n_zones)
    
    # Store indices data for visualization
    result['indices_data'] = {name: data[0] for name, data in indices_results.items()}
    
    # Visualize
    if visualize:
        viz_path = output_path / "management_zones.png"
        analyzer.visualize_zones(result, save_path=str(viz_path))
    
    # Export
    if export:
        geotiff_path = output_path / "management_zones.tif"
        analyzer.export_zones(result, geotiff_path)
    
    # Generate recommendations
    analyzer.recommend_actions(result)
    
    print(f"\n{'='*70}")
    print(f"Zone analysis complete! Outputs saved to: {output_path.absolute()}")
    print(f"{'='*70}\n")
    
    return result


if __name__ == "__main__":
    print("Management Zone Analysis Module")
    print("=" * 70)
    print(__doc__)

