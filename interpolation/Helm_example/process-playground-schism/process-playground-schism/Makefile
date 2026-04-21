PROJECT_NAME = process-playground

.ONESHELL:
.SHELLFLAGS = -ec
SHELL := /bin/bash

ENVIRONMENT_FILE_NAME = environment.yaml

define conda-command
	micromamba $1 || mamba $1 || conda $1
endef

create-environment:
	$(call conda-command, env create --file ${ENVIRONMENT_FILE_NAME} --name ${PROJECT_NAME})

check-format:
	pre-commit run --all-files --show-diff-on-failure

lint-charts:
	$${PWD}/charts.sh lint_all_charts

publish-charts:
	$${PWD}/charts.sh publish_all_charts
