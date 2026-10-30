#!/usr/bin/env python3
"""
Scriptized version of validation toolbox in-situ validation workflow.

Usage:
    prepare the config.yaml with file and validation settings
    python run_validation.py --config /full/path/to/config.yaml

Notes:
- Run this script from within (or with access to) the toolbox repository.
- Ensure the validation-toolbox environment is activated.
"""

import os
import sys
import argparse
import logging
import yaml
from pathlib import Path

## go to validation-toolbox repo
toolbox_dir='/gpfs/work/ksddata/wp3-validation-toolbox/'
os.chdir(toolbox_dir)

# run command
#python run_validation.py --config /gpfs/work/ksddata/wp3-validation-toolbox/configs/ugrid_ssh_seyed.yaml

# -----------------------------------------------------------------------------
# Utility: Find repository root
# -----------------------------------------------------------------------------
def find_repo_root(start: Path) -> Path:
    """Search upwards for repository root (identified by toolbox/config.py)."""
    for candidate in [start, *start.parents]:
        if (candidate / "toolbox" / "config.py").exists():
            return candidate
    return start


# -----------------------------------------------------------------------------
# Setup environment and imports
# -----------------------------------------------------------------------------
def setup_imports():
    """Ensure local toolbox is imported instead of site-packages version."""
    repo_root = find_repo_root(Path.cwd().resolve())

    # Prioritize local repo
    if sys.path:
        sys.path = [p for p in sys.path if Path(p).resolve() != repo_root]
    sys.path.insert(0, str(repo_root))

    # Clear cached toolbox modules
    for mod_name in list(sys.modules.keys()):
        if mod_name == "toolbox" or mod_name.startswith("toolbox."):
            del sys.modules[mod_name]

    return repo_root


# -----------------------------------------------------------------------------
# Configure logging
# -----------------------------------------------------------------------------
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
    )
    return logging.getLogger("validation_script")


# -----------------------------------------------------------------------------
# Load configuration
# -----------------------------------------------------------------------------
def load_config(config_path, logger):
    from toolbox.config import ValidationConfig

    try:
        with open(config_path, "r") as f:
            config_dict = yaml.safe_load(f)

        logger.info(f"Loaded configuration from: {config_path}")

        config = ValidationConfig.load_file(str(config_path))

        # Force non-interactive execution
        config.skip_confirmation = True
        config.yes = True

        logger.info(f"Model URL: {config.model_url}")
        logger.info(f"Observations URI: {config.observations_uri}")
        logger.info(f"Output Directory: {config.output_dir}")
        logger.info(f"insitu_type: {config.insitu_type}")
        logger.info(f"obs_var: {config.obs_var}")

        return config

    except FileNotFoundError:
        logger.error(f"Configuration file not found: {config_path}")
    except Exception as e:
        logger.error(f"Error loading configuration: {e}", exc_info=True)

    return None


# -----------------------------------------------------------------------------
# Validation workflow
# -----------------------------------------------------------------------------
def run_validation_workflow(config, logger):
    from toolbox.insitu_mooring_validation import (
        setup_environment,
        load_datasets,
        process_data,
        compute_validation_metrics,
        save_validation_results,
    )

    validation_successful = False
    validation_results_dir = None

    try:
        logger.info(f"Starting {config.dataset_name.upper()} validation workflow")

        # Step 1: Setup environment
        config.output_dir = setup_environment(config.output_dir)
        validation_results_dir = config.output_dir
        logger.info(f"Output directory: {validation_results_dir}")

        # Step 2-3: Load datasets
        logger.info("Loading datasets...")
        model_dataset, obs_data, model_domain, variable_units = load_datasets(config)
        logger.info("Datasets loaded successfully.")

        # Step 6-7: Process data
        logger.info("Processing data...")
        model_insitu, obs_data_arr, obs_data_resample = process_data(
            model_dataset, obs_data, config
        )
        logger.info("Processing complete.")

        # Step 8: Compute metrics
        logger.info("Computing validation metrics...")
        validation_data = compute_validation_metrics(
            model_insitu, obs_data_arr, obs_data_resample, variable_units
        )
        logger.info("Metrics computed.")

        # Step 9: Save results
        logger.info("Saving results...")
        metrics_file, insitu_file, summary_file = save_validation_results(
            validation_data, config, variable_units
        )

        logger.info(f"Metrics: {metrics_file}")
        logger.info(f"Datasets: {insitu_file}")
        logger.info(f"Summary: {summary_file}")

        logger.info("Validation completed successfully.")
        validation_successful = True

    except Exception as e:
        logger.error(f"Validation failed: {e}", exc_info=True)

    return validation_successful, validation_results_dir


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def run_plotting_step(config, results_dir, logger):
    from toolbox.insitu_plots import process_validation_results as run_plotting

    input_dir = results_dir
    plot_output_dir = os.path.join(results_dir, "plots")
    dataset_name = config.dataset_name
    plot_format = "png"

    try:
        logger.info("Starting plotting...")
        logger.info(f"Input: {input_dir}")
        logger.info(f"Output: {plot_output_dir}")

        run_plotting(
            input_dir=input_dir,
            output_dir=plot_output_dir,
            dataset_name=dataset_name,
            save_format=plot_format,
        )

        logger.info("Plotting complete.")

    except FileNotFoundError as e:
        logger.error(f"Plotting failed (missing files): {e}")
    except Exception as e:
        logger.error(f"Plotting failed: {e}", exc_info=True)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Run validation workflow")
    parser.add_argument(
        "--config",
        required=True,
        help="Full path to configuration YAML file",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plotting step",
    )

    args = parser.parse_args()

    logger = setup_logging()
    repo_root = setup_imports()

    logger.info(f"Repo root: {repo_root}")

    # Import toolbox AFTER path setup
    try:
        import toolbox
        logger.info(f"toolbox imported from: {toolbox.__file__}")
    except ImportError as e:
        logger.error(f"Failed to import toolbox: {e}")
        sys.exit(1)

    # Load config
    config = load_config(args.config, logger)
    if not config:
        logger.error("Exiting due to configuration error.")
        sys.exit(1)

    # Run validation
    success, results_dir = run_validation_workflow(config, logger)

    # Plotting
    if success and results_dir and not args.no_plots:
        run_plotting_step(config, results_dir, logger)
    elif not success:
        logger.warning("Skipping plotting due to failed validation.")


if __name__ == "__main__":
    main()