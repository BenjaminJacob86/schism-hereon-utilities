# SCHISM Docker Container - Local Testing Guide

This guide explains how to build, run, and test the SCHISM Docker container interactively on your Windows computer for debugging before deploying to the Helm chart.

## Prerequisites

- **Docker Desktop** installed and running on Windows
- **PowerShell** or **Command Prompt**
- Access to S3/MinIO credentials (for downloading grid and forcing files)
- Grid and forcing files already uploaded to S3 (from the Grid+Forcing generator)

## Quick Start

### 1. Build the Docker Image

Open PowerShell in the SCHISM2 directory and run:

```powershell
docker build -t schism:local .
```

**Note:** The build process compiles SCHISM twice:
- First build: Creates `pschism_PREC_EVAP_TVD-SB` (with atmospheric forcing)
- Second build: Creates `pschism_TVD-SB` (without atmospheric forcing)

Both executables will be in `/SCHISM/GBconfig/bin/`.

### 2. Test Run (Interactive Mode)

Run the container interactively to debug:

```powershell
docker run -it --rm `
  -v ${PWD}/test_output:/output-data `
  -e GRID_FOLDER="NWS_NEW" `
  -e FORCING_FOLDER="NWS_NEW/forcing_output" `
  -e AWS_S3_ENDPOINT="minio.dive.edito.eu" `
  -e AWS_ACCESS_KEY_ID="YOUR_ACCESS_KEY" `
  -e AWS_SECRET_ACCESS_KEY="YOUR_SECRET_KEY" `
  -e AWS_SESSION_TOKEN="YOUR_SESSION_TOKEN" `
  -e User_Name="user-jacobb" `
  -e EDITO_INFRA_OUTPUT="/output-data" `
  schism:local /bin/bash
```

**Note:** Replace the S3 credentials with your actual values.

### 3. Manual Step-by-Step Testing

Once inside the container, you can run each step manually:

```bash
# Navigate to setup directory
cd /app/setup

# Run the download script manually
/opt/miniforge/envs/forcing/bin/python /app/dowload_forcings.py

# Check what was downloaded
ls -la

# Check if sflux folder exists
ls -la sflux/  # (if present)

# Check param.nml was updated correctly
cat param.nml | grep nws

# Run exec.sh manually
bash exec.sh

# Check outputs
ls -la outputs/
```

## Environment Variables

### Required for Production (EDITO Platform)

These are automatically set by the Helm chart:

- `GRID_FOLDER` - S3 folder path for grid files (e.g., `EFWS`)
- `FORCING_FOLDER` - S3 folder path for forcing files (e.g., `EFWS/forcing_output`)
- `AWS_S3_ENDPOINT` - S3 endpoint URL (e.g., `minio.dive.edito.eu`)
- `AWS_ACCESS_KEY_ID` - S3 access key
- `AWS_SECRET_ACCESS_KEY` - S3 secret key
- `AWS_SESSION_TOKEN` - S3 session token
- `User_Name` - User identifier (e.g., `user-jacobb`)
- `EDITO_INFRA_OUTPUT` - Output directory path (e.g., `/output-data`)

### Local Testing Mode

If `EDITO_INFRA_OUTPUT` is not set, the script runs in local test mode:
- Uses hardcoded S3 credentials (in `dowload_forcings.py`)
- Outputs to `./outputs/` directory
- Uses hardcoded `grid_selection` and `forcing_selection` values

## Testing Scenarios

### Scenario 1: Test WITHOUT Atmospheric Forcing (no sflux)

```powershell
docker run -it --rm `
  -v ${PWD}/test_output:/output-data `
  -e GRID_FOLDER="YOUR_GRID_FOLDER" `
  -e FORCING_FOLDER="YOUR_FORCING_FOLDER" `
  -e AWS_S3_ENDPOINT="minio.dive.edito.eu" `
  -e AWS_ACCESS_KEY_ID="YOUR_KEY" `
  -e AWS_SECRET_ACCESS_KEY="YOUR_SECRET" `
  -e AWS_SESSION_TOKEN="YOUR_TOKEN" `
  -e User_Name="user-jacobb" `
  -e EDITO_INFRA_OUTPUT="/output-data" `
  schism:local
```

**Expected behavior:**
- Downloads grid files
- Downloads forcing files (no `sflux/` folder)
- Sets `nws=0` in `param.nml`
- Uses `pschism_TVD-SB` binary
- Creates `exec.sh` that selects the non-atmospheric binary

### Scenario 2: Test WITH Atmospheric Forcing (with sflux)

```powershell
# Same command as above, but ensure your FORCING_FOLDER contains sflux/ subdirectory
docker run -it --rm `
  -v ${PWD}/test_output:/output-data `
  -e GRID_FOLDER="YOUR_GRID_FOLDER" `
  -e FORCING_FOLDER="YOUR_FORCING_FOLDER_WITH_SFLUX" `
  -e AWS_S3_ENDPOINT="minio.dive.edito.eu" `
  -e AWS_ACCESS_KEY_ID="YOUR_KEY" `
  -e AWS_SECRET_ACCESS_KEY="YOUR_SECRET" `
  -e AWS_SESSION_TOKEN="YOUR_TOKEN" `
  -e User_Name="user-jacobb" `
  -e EDITO_INFRA_OUTPUT="/output-data" `
  schism:local
```

**Expected behavior:**
- Downloads grid files
- Downloads forcing files including `sflux/` folder
- Sets `nws=2` in `param.nml`
- Uses `pschism_PREC_EVAP_TVD-SB` binary
- Creates `exec.sh` that selects the atmospheric binary

## Debugging Tips

### 1. Check Downloaded Files

```bash
# Inside container
cd /app/setup
ls -la
ls -la sflux/  # Check if sflux exists
```

### 2. Verify param.nml Configuration

```bash
# Check nws parameter
cat param.nml | grep nws

# Should show:
# nws = 0  (without sflux)
# nws = 2  (with sflux)
```

### 3. Check exec.sh Script

```bash
cat exec.sh

# Should show:
# - pschism_PREC_EVAP_TVD-SB (if sflux exists)
# - pschism_TVD-SB (if sflux doesn't exist)
```

### 4. Test Binary Selection Manually

```bash
# Check which binary exists
ls -la /SCHISM/GBconfig/bin/

# Test binary selection logic
if [ -d "sflux" ] && [ "$(ls -A sflux 2>/dev/null)" ]; then
    echo "sflux found - use PREC_EVAP binary"
else
    echo "sflux not found - use standard binary"
fi
```

### 5. Run Python Script with Debug Output

```bash
# Run with verbose output
cd /app/setup
/opt/miniforge/envs/forcing/bin/python /app/dowload_forcings.py

# Check for errors in download
# Look for "sflux available" or "sflux not available" messages
```

### 6. Check S3 Download Structure

The download function preserves directory structure:
- `EFWS/forcing_output/sflux/file.nc` → downloads to `./sflux/file.nc`
- `EFWS_forcing_output/sflux/file.nc` → downloads to `./sflux/file.nc`

Verify files are in the correct location:

```bash
find . -name "*.nc" -type f
```

### 7. Monitor Container Logs

```powershell
# Run container and see all output
docker run --rm `
  -v ${PWD}/test_output:/output-data `
  -e GRID_FOLDER="YOUR_GRID" `
  -e FORCING_FOLDER="YOUR_FORCING" `
  # ... other env vars ...
  schism:local
```

## Common Issues and Solutions

### Issue: "sflux directory not found" but files exist

**Solution:** Check that files were actually downloaded:
```bash
ls -la sflux/
# If empty or missing, check S3 path in FORCING_FOLDER
```

### Issue: Wrong binary selected

**Solution:** The `exec.sh` checks for `sflux/` directory. Verify:
```bash
# Check if sflux exists and is not empty
[ -d "sflux" ] && [ "$(ls -A sflux 2>/dev/null)" ] && echo "sflux OK" || echo "sflux missing"
```

### Issue: S3 download fails

**Solution:** Verify credentials and paths:
```bash
# Test S3 connection (inside container)
python3 -c "import boto3; s3 = boto3.client('s3', endpoint_url='https://minio.dive.edito.eu', ...); print(s3.list_buckets())"
```

### Issue: Files downloaded but in wrong location

**Solution:** The download function preserves relative paths. Check:
```bash
# Files should be relative to /app/setup/
pwd
ls -R
```

## Testing Before Helm Chart Deployment

### Step 1: Test Locally with Real Data

1. Ensure grid and forcing files are uploaded to S3
2. Test with actual S3 paths you'll use in production
3. Verify both scenarios (with/without sflux)

### Step 2: Verify Output Structure

```bash
# Check outputs directory structure
ls -la /output-data/
# Should contain:
# - *.nc files (SCHISM outputs)
# - description.txt
# - schism_edito.py
# - plot_outputs.ipynb
```

### Step 3: Test Environment Variable Handling

Test that the container correctly reads all environment variables:

```powershell
docker run -it --rm `
  -e GRID_FOLDER="TEST_GRID" `
  -e FORCING_FOLDER="TEST_FORCING" `
  # ... other vars ...
  schism:local /bin/bash

# Inside container, check env vars
env | grep -E "GRID|FORCING|AWS|User|EDITO"
```

### Step 4: Validate exec.sh Creation

The `exec.sh` script is created dynamically. Verify it's correct:

```bash
# After running download script
cat exec.sh
chmod +x exec.sh
bash -n exec.sh  # Syntax check
```

## File Structure

```
/app/
├── dowload_forcings.py    # Main download and setup script
├── entrypoint.sh          # Container entrypoint
├── schism_edito.py        # Post-processing script
├── plot_outputs.ipynb     # Visualization notebook
└── setup/
    └── param.nml          # SCHISM parameter file (updated by script)

/SCHISM/
└── GBconfig/
    └── bin/
        ├── pschism_TVD-SB              # Standard binary
        └── pschism_PREC_EVAP_TVD-SB   # With atmospheric forcing
```

## Workflow Summary

1. **Download Phase** (`dowload_forcings.py`):
   - Downloads grid files from S3
   - Downloads forcing files from S3 (preserves `sflux/` if present)
   - Updates `param.nml` (sets `nws=2` if sflux exists, else `nws=0`)
   - Creates `exec.sh` with appropriate binary selection

2. **Execution Phase** (`exec.sh`):
   - Checks for `sflux/` directory
   - Selects correct SCHISM binary
   - Runs SCHISM simulation

3. **Output Phase** (`entrypoint.sh`):
   - Copies NetCDF outputs to output directory
   - Copies post-processing scripts

## Next Steps

After successful local testing:

1. Update Helm chart with correct environment variables
2. Test on EDITO platform with a small test case
3. Monitor logs for any issues
4. Scale up to production runs

## Additional Resources

- SCHISM documentation: https://github.com/schism-dev/schism
- Docker documentation: https://docs.docker.com/
- EDITO platform documentation (internal)

