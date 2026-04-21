#!/bin/bash
cd /app/setup/

# Run the Python script
/opt/miniforge/envs/forcing/bin/python /app/dowload_forcings.py

# Wait for exec.sh to be created
while [ ! -f /app/setup/exec.sh ]; do
    echo "Waiting for exec.sh to be created..."
    sleep 1
done

# Make it executable just in case
chmod +x /app/setup/exec.sh

# Run the script
bash /app/setup/exec.sh

# Final step — copy NetCDF outputs to the output path
cp /app/setup/outputs/*.nc ${EDITO_INFRA_OUTPUT}/
cp /app/setup/schism_edito.py ${EDITO_INFRA_OUTPUT}/
cp /app/setup/plot_outputs.ipynb ${EDITO_INFRA_OUTPUT}/
