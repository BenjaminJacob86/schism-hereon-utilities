"""
Convert NetCDF files to Zarr format for cloud applications.

This script converts structured grid NetCDF files (created by interpolation_hybrid.py)
to Zarr format, which is optimized for cloud storage and parallel access.

Features:
    - Converts one or multiple NetCDF files to Zarr
    - Option to delete original NetCDF files after successful conversion
    - Preserves all metadata and attributes
    - Supports chunking optimization for cloud storage
"""

import os
import sys
import argparse
from glob import glob
import xarray as xr
import zarr

# Try to import Blosc from numcodecs (required for compression)
try:
    from numcodecs import Blosc
except ImportError:
    try:
        # Fallback: some zarr versions may have Blosc directly
        from zarr import Blosc
    except ImportError:
        # If Blosc is not available, we'll use a different compressor or None
        Blosc = None
        print("Warning: Blosc compressor not available. Compression may be limited.")


def merge_and_convert_netcdf_to_zarr(nc_files, zarr_file, delete_nc=False, chunk_sizes=None, use_default_chunks=True):
    """
    Merge multiple NetCDF files along the time dimension and convert to a single Zarr file.
    
    Parameters:
    -----------
    nc_files : list
        List of paths to input NetCDF files (will be sorted and concatenated along time)
    zarr_file : str
        Path to output Zarr file/directory
    delete_nc : bool
        If True, delete original NetCDF files after successful conversion
    chunk_sizes : dict, optional
        Dictionary of chunk sizes for dimensions
    use_default_chunks : bool
        If True and chunk_sizes is None, applies default chunking
    """
    print(f"Merging {len(nc_files)} NetCDF files into single Zarr: {zarr_file}")
    print(f"  Input files: {len(nc_files)} files")
    
    try:
        # Open all datasets lazily
        datasets = [xr.open_dataset(f) for f in nc_files]
        
        # Concatenate along time dimension
        # xarray will automatically align coordinates
        print("  Concatenating datasets along time dimension...")
        ds_merged = xr.concat(datasets, dim='time', data_vars='minimal', coords='minimal', compat='override')
        
        # Close individual datasets
        for ds in datasets:
            ds.close()
        
        print(f"  Merged dataset: {len(ds_merged.time)} time steps, {len(ds_merged.data_vars)} variables")
        
        # Apply default chunking if not specified and defaults are enabled
        if chunk_sizes is None and use_default_chunks:
            default_chunks = {}
            if 'time' in ds_merged.dims:
                default_chunks['time'] = min(24, len(ds_merged.time))
            if 'latitude' in ds_merged.dims:
                default_chunks['latitude'] = min(100, len(ds_merged.latitude))
            if 'longitude' in ds_merged.dims:
                default_chunks['longitude'] = min(100, len(ds_merged.longitude))
            if 'depth' in ds_merged.dims:
                default_chunks['depth'] = len(ds_merged.depth)
            
            if default_chunks:
                chunk_sizes = default_chunks
                print(f"  Using default chunk sizes: {chunk_sizes}")
        
        # Set chunk sizes if provided
        if chunk_sizes:
            # Build chunks dict with only dimensions that exist in the dataset
            # xarray's chunk() expects dimension names as keys
            chunks_dict = {}
            for dim in ds_merged.dims:
                if dim in chunk_sizes:
                    chunks_dict[dim] = chunk_sizes[dim]
                else:
                    # Use full dimension size if not specified
                    chunks_dict[dim] = len(ds_merged[dim])
            # Apply chunking to the dataset
            ds_merged = ds_merged.chunk(chunks_dict)
        
        # Write merged dataset to Zarr
        print("  Writing merged dataset to Zarr...")
        ds_merged.to_zarr(zarr_file, mode='w')
        ds_merged.close()
        
        print(f"  ✓ Successfully merged and converted to: {zarr_file}")
        
        # Delete NetCDF files if requested
        if delete_nc:
            for nc_file in nc_files:
                os.remove(nc_file)
                print(f"  ✓ Deleted original NetCDF file: {nc_file}")
        
        return True
        
    except Exception as e:
        print(f"  ✗ Error merging and converting: {e}")
        import traceback
        traceback.print_exc()
        return False


def convert_netcdf_to_zarr(nc_file, zarr_file, delete_nc=False, chunk_sizes=None, use_default_chunks=True):
    """
    Convert a NetCDF file to Zarr format.
    
    Parameters:
    -----------
    nc_file : str
        Path to input NetCDF file
    zarr_file : str
        Path to output Zarr file/directory
    delete_nc : bool
        If True, delete the NetCDF file after successful conversion
    chunk_sizes : dict, optional
        Dictionary of chunk sizes for dimensions (e.g., {'time': 24, 'lat': 100, 'lon': 100})
        If None and use_default_chunks=True, uses sensible defaults for cloud storage
        If None and use_default_chunks=False, lets zarr determine optimal chunks
    use_default_chunks : bool
        If True and chunk_sizes is None, applies default chunking optimized for cloud access
    """
    print(f"Converting: {nc_file} -> {zarr_file}")
    
    try:
        # Open NetCDF file
        ds = xr.open_dataset(nc_file)
        
        # Apply default chunking if not specified and defaults are enabled
        if chunk_sizes is None and use_default_chunks:
            # Default chunk sizes optimized for cloud storage and common access patterns
            # These are reasonable defaults but can be overridden
            default_chunks = {}
            if 'time' in ds.dims:
                default_chunks['time'] = min(24, len(ds.time))  # Daily chunks if available
            if 'latitude' in ds.dims:
                default_chunks['latitude'] = min(100, len(ds.latitude))  # Reasonable spatial chunk
            if 'longitude' in ds.dims:
                default_chunks['longitude'] = min(100, len(ds.longitude))  # Reasonable spatial chunk
            if 'depth' in ds.dims:
                default_chunks['depth'] = len(ds.depth)  # Full depth dimension
            
            if default_chunks:
                chunk_sizes = default_chunks
                print(f"  Using default chunk sizes: {chunk_sizes}")
        
        # Set chunk sizes if provided
        if chunk_sizes:
            # Build chunks dict with only dimensions that exist in the dataset
            # xarray's chunk() expects dimension names as keys
            chunks_dict = {}
            for dim in ds.dims:
                if dim in chunk_sizes:
                    chunks_dict[dim] = chunk_sizes[dim]
                else:
                    # Use full dimension size if not specified
                    chunks_dict[dim] = len(ds[dim])
            # Apply chunking to the dataset
            ds = ds.chunk(chunks_dict)
        
        # Convert to Zarr
        # xarray's to_zarr method will automatically use appropriate compression
        # Zarr will use Blosc if available (via numcodecs), otherwise fall back to other compressors
        # We don't need to specify compressor explicitly - let zarr/xarray choose defaults
        # This avoids the "Expected a BytesBytesCodec" error
        
        # Write to Zarr without explicit encoding/compression settings
        # xarray will preserve all metadata and use sensible defaults for compression and chunking
        ds.to_zarr(zarr_file, mode='w')
        ds.close()
        
        print(f"  ✓ Successfully converted to: {zarr_file}")
        
        # Delete NetCDF file if requested
        if delete_nc:
            os.remove(nc_file)
            print(f"  ✓ Deleted original NetCDF file: {nc_file}")
        
        return True
        
    except Exception as e:
        print(f"  ✗ Error converting {nc_file}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Convert NetCDF files to Zarr format for cloud applications',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert single file
  python convert_to_zarr.py output.nc output.zarr
  
  # Convert all NetCDF files in directory (individual files)
  python convert_to_zarr.py --input-dir ./outputs/ --pattern "*.nc"
  
  # Merge all files into one Zarr file
  python convert_to_zarr.py --input-dir ./outputs/ --pattern "*.nc" --merge
  
  # Merge with custom output filename
  python convert_to_zarr.py --input-dir ./outputs/ --pattern "*.nc" --merge --merge-output merged_data.zarr
  
  # Convert and delete originals
  python convert_to_zarr.py --input-dir ./outputs/ --delete-nc
  
  # Convert with custom chunking
  python convert_to_zarr.py output.nc output.zarr --chunks time=24 lat=100 lon=100
        """
    )
    
    parser.add_argument('input', nargs='?', help='Input NetCDF file (or use --input-dir)')
    parser.add_argument('output', nargs='?', help='Output Zarr file/directory (or use --input-dir)')
    
    parser.add_argument('--input-dir', type=str, help='Directory containing NetCDF files to convert')
    parser.add_argument('--pattern', type=str, default='*.nc', 
                       help='File pattern to match (default: *.nc)')
    parser.add_argument('--output-dir', type=str, 
                       help='Output directory for Zarr files (default: same as input-dir)')
    parser.add_argument('--delete-nc', action='store_true',
                       help='Delete original NetCDF files after successful conversion')
    parser.add_argument('--chunks', nargs='+', metavar='DIM=SIZE',
                       help='Chunk sizes for dimensions (e.g., time=24 lat=100 lon=100)')
    parser.add_argument('--no-default-chunks', action='store_true',
                       help='Disable default chunking, let zarr determine optimal chunks automatically')
    parser.add_argument('--merge', action='store_true',
                       help='Merge all input files into a single Zarr file (concatenate along time dimension)')
    parser.add_argument('--merge-output', type=str,
                       help='Output filename for merged Zarr file (only used with --merge)')
    
    args = parser.parse_args()
    
    # Parse chunk sizes if provided
    chunk_sizes = None
    if args.chunks:
        chunk_sizes = {}
        for chunk_arg in args.chunks:
            if '=' in chunk_arg:
                dim, size = chunk_arg.split('=', 1)
                try:
                    chunk_sizes[dim] = int(size)
                except ValueError:
                    print(f"Warning: Invalid chunk size '{chunk_arg}', ignoring")
            else:
                print(f"Warning: Invalid chunk format '{chunk_arg}', use DIM=SIZE")
    
    # Process files
    if args.input_dir:
        # Batch processing mode
        if not os.path.isdir(args.input_dir):
            print(f"Error: Input directory does not exist: {args.input_dir}")
            sys.exit(1)
        
        output_dir = args.output_dir if args.output_dir else args.input_dir
        
        # Find all matching files
        pattern = os.path.join(args.input_dir, args.pattern)
        nc_files = sorted(glob(pattern))
        
        if not nc_files:
            print(f"No files found matching pattern: {pattern}")
            sys.exit(1)
        
        if args.merge:
            # Merge mode: combine all files into one Zarr
            if args.merge_output:
                zarr_file = os.path.join(output_dir, args.merge_output)
            else:
                # Generate default merged filename from pattern
                # Extract base name from pattern (remove wildcards)
                base_pattern = args.pattern.replace('*', '').replace('.nc', '')
                if base_pattern:
                    zarr_file = os.path.join(output_dir, f"{base_pattern}_merged.zarr")
                else:
                    zarr_file = os.path.join(output_dir, "merged.zarr")
            
            print(f"Found {len(nc_files)} file(s) to merge")
            print(f"Output file: {zarr_file}\n")
            
            use_defaults = not args.no_default_chunks
            success = merge_and_convert_netcdf_to_zarr(nc_files, zarr_file, args.delete_nc, chunk_sizes, use_defaults)
            
            print(f"\n{'='*70}")
            if success:
                print(f"Merge and conversion complete: 1 merged Zarr file created")
            else:
                print(f"Merge and conversion failed")
            print(f"{'='*70}\n")
            
            sys.exit(0 if success else 1)
        else:
            # Individual file mode: convert each file separately
            print(f"Found {len(nc_files)} file(s) to convert")
            print(f"Output directory: {output_dir}\n")
            
            success_count = 0
            for nc_file in nc_files:
                # Generate output filename
                base_name = os.path.splitext(os.path.basename(nc_file))[0]
                zarr_file = os.path.join(output_dir, base_name + '.zarr')
                
                use_defaults = not args.no_default_chunks
                if convert_netcdf_to_zarr(nc_file, zarr_file, args.delete_nc, chunk_sizes, use_defaults):
                    success_count += 1
            
            print(f"\n{'='*70}")
            print(f"Conversion complete: {success_count}/{len(nc_files)} files converted successfully")
            print(f"{'='*70}\n")
        
    elif args.input and args.output:
        # Single file mode
        if not os.path.exists(args.input):
            print(f"Error: Input file does not exist: {args.input}")
            sys.exit(1)
        
        use_defaults = not args.no_default_chunks
        success = convert_netcdf_to_zarr(args.input, args.output, args.delete_nc, chunk_sizes, use_defaults)
        sys.exit(0 if success else 1)
        
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()

