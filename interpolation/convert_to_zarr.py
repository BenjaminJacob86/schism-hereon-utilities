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
            # Apply chunking to all data variables
            for var in ds.data_vars:
                chunks = {}
                for dim in ds[var].dims:
                    if dim in chunk_sizes:
                        chunks[dim] = chunk_sizes[dim]
                    else:
                        chunks[dim] = len(ds[dim])  # Use full dimension if not specified
                ds[var] = ds[var].chunk(chunks)
        
        # Convert to Zarr
        # Use zarr backend with appropriate encoding
        encoding = {}
        for var in ds.data_vars:
            encoding[var] = {
                'compressor': zarr.Blosc(cname='lz4', clevel=5, shuffle=1),
                'chunks': None  # Let zarr determine optimal chunks
            }
        
        # Write to Zarr
        ds.to_zarr(zarr_file, mode='w', encoding=encoding)
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
  
  # Convert all NetCDF files in directory
  python convert_to_zarr.py --input-dir ./outputs/ --pattern "*.nc"
  
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

