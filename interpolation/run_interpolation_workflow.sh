#!/bin/bash

###############################################################################
# SCHISM Interpolation Workflow Script
#
# This script runs the complete interpolation workflow:
#   1. Interpolates SCHISM unstructured grid output to structured grid
#   2. Optionally converts NetCDF output to Zarr format
#   3. Optionally deletes NetCDF files after Zarr conversion
#
# Usage:
#   ./run_interpolation_workflow.sh [OPTIONS]
#
# Options:
#   --interp-script PATH      Path to interpolation script (default: interpolation_hybrid.py)
#   --convert-script PATH     Path to conversion script (default: convert_to_zarr.py)
#   --convert                 Enable Zarr conversion after interpolation
#   --delete-nc               Delete NetCDF files after successful Zarr conversion
#   --chunks DIM=SIZE         Chunk sizes for Zarr (e.g., time=24 lat=100 lon=100)
#   --help                    Show this help message
#
# Examples:
#   # Run interpolation only
#   ./run_interpolation_workflow.sh
#
#   # Run interpolation and convert to Zarr
#   ./run_interpolation_workflow.sh --convert
#
#   # Run full workflow with conversion and cleanup
#   ./run_interpolation_workflow.sh --convert --delete-nc
#
#   # Run with custom chunking
#   ./run_interpolation_workflow.sh --convert --chunks time=24 lat=100 lon=100
###############################################################################

set -e  # Exit on error

# Default values
INTERP_SCRIPT="interpolation_hybrid.py"
CONVERT_SCRIPT="convert_to_zarr.py"
DO_CONVERT=false
DELETE_NC=false
CHUNKS=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show help
show_help() {
    head -n 30 "$0" | grep -E "^# " | sed 's/^# //' | sed 's/^#//'
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --interp-script)
            INTERP_SCRIPT="$2"
            shift 2
            ;;
        --convert-script)
            CONVERT_SCRIPT="$2"
            shift 2
            ;;
        --convert)
            DO_CONVERT=true
            shift
            ;;
        --delete-nc)
            DELETE_NC=true
            shift
            ;;
        --chunks)
            CHUNKS="$2"
            shift 2
            ;;
        --help|-h)
            show_help
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate that scripts exist
if [ ! -f "$INTERP_SCRIPT" ]; then
    print_error "Interpolation script not found: $INTERP_SCRIPT"
    exit 1
fi

if [ "$DO_CONVERT" = true ] && [ ! -f "$CONVERT_SCRIPT" ]; then
    print_error "Conversion script not found: $CONVERT_SCRIPT"
    exit 1
fi

# Check if Python is available
if ! command -v python &> /dev/null && ! command -v python3 &> /dev/null; then
    print_error "Python not found. Please install Python to run this script."
    exit 1
fi

# Use python3 if available, otherwise python
if command -v python3 &> /dev/null; then
    PYTHON_CMD="python3"
else
    PYTHON_CMD="python"
fi

print_info "Starting SCHISM interpolation workflow"
echo "=========================================="
print_info "Interpolation script: $INTERP_SCRIPT"
if [ "$DO_CONVERT" = true ]; then
    print_info "Conversion script: $CONVERT_SCRIPT"
    print_info "Zarr conversion: ENABLED"
    if [ "$DELETE_NC" = true ]; then
        print_info "Delete NetCDF after conversion: ENABLED"
    fi
    if [ -n "$CHUNKS" ]; then
        print_info "Custom chunking: $CHUNKS"
    fi
else
    print_info "Zarr conversion: DISABLED"
fi
echo "=========================================="
echo ""

# Step 1: Run interpolation
print_info "Step 1: Running interpolation..."
echo ""

if $PYTHON_CMD "$INTERP_SCRIPT"; then
    print_success "Interpolation completed successfully"
else
    print_error "Interpolation failed. Exiting."
    exit 1
fi

echo ""

# Step 2: Convert to Zarr (if requested)
if [ "$DO_CONVERT" = true ]; then
    print_info "Step 2: Converting NetCDF files to Zarr format..."
    echo ""
    
    # Build conversion command
    CONVERT_CMD="$PYTHON_CMD $CONVERT_SCRIPT --input-dir . --pattern 'schism-wwm*.nc'"
    
    if [ "$DELETE_NC" = true ]; then
        CONVERT_CMD="$CONVERT_CMD --delete-nc"
    fi
    
    if [ -n "$CHUNKS" ]; then
        CONVERT_CMD="$CONVERT_CMD --chunks $CHUNKS"
    fi
    
    if eval $CONVERT_CMD; then
        print_success "Zarr conversion completed successfully"
    else
        print_error "Zarr conversion failed."
        exit 1
    fi
    
    echo ""
fi

# Summary
echo "=========================================="
print_success "Workflow completed successfully!"
echo "=========================================="

if [ "$DO_CONVERT" = true ]; then
    print_info "Output files:"
    echo "  - NetCDF files: *.nc"
    echo "  - Zarr files: *.zarr"
    if [ "$DELETE_NC" = true ]; then
        echo "  (NetCDF files have been deleted)"
    fi
else
    print_info "Output files: *.nc"
fi

